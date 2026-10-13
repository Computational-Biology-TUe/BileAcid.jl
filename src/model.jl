# Rate constants in the form used by the right-hand side, derived once from `Parameters`
# and `Constants`. This is where the structure of the physiology is configured: which gut
# segments carry the bile acid transporter (ASBT) or bacteria, which biotransformations
# occur, and how hepatic extraction depends on species and conjugation.

"Small-intestinal segments with the proximal bacterial activity `k_bact_psi`."
const PROXIMAL_SI = 1:5
"Small-intestinal segments with the slower transit `k_si_2` (terminal ileum)."
const TERMINAL_ILEUM = 9:10
"Small-intestinal segments with active uptake via the apical sodium-dependent BA transporter."
const ASBT_SEGMENTS = 9:10

"""
    GutSegment

Properties of one segment of the intestinal tract.
"""
Base.@kwdef struct GutSegment{T}
    k_transit::T                    # basal transit to the next segment (1/min), scaled by the reflex
    bacteria::T                     # bacterial activity: scales deconjugation and biotransformation
    k_passive::SVector{N_FORMS, T}  # passive uptake per conjugation form (1/min)
    asbt::Bool                      # active uptake by ASBT
    biotransformation::Bool         # bacterial biotransformation of unconjugated BA
end

"""
    BileAcidModel(p = Parameters(), c = Constants(); meal_times = Float64[])

The model in the form used by `rhs!`: rate constants derived from parameters `p`
and constants `c`, and the start times of meals (min), which trigger the postprandial
increase of gallbladder emptying and intestinal transit.
"""
struct BileAcidModel{T, P <: Parameters, C <: Constants}
    # Liver
    synthesis::BAMatrix{T}                    # de novo synthesis per species and form (amount/min)
    hepatic_extraction::BAMatrix{T}           # fraction extracted from blood per pass
    frac_tauro::T                             # fraction conjugated with taurine (rest: glycine)
    k_secretion::T                            # secretion into bile (1/min)
    k_spill::T                                # spill-over into plasma (1/min)
    liver_biotransformation::SMatrix{N_SPECIES, N_SPECIES, T, N_SPECIES^2}
    # Gallbladder
    frac_gb::T                                # fraction of secreted bile stored in the gallbladder
    k_gb_emptying::T                          # basal emptying (1/min), scaled by the reflex
    # Intestinal tract
    segments::SVector{N_GUT, GutSegment{T}}
    k_deconjugation::SVector{N_FORMS, T}      # per conjugation form and unit bacterial activity (1/min)
    vmax_asbt::SVector{N_SPECIES, T}          # maximal ASBT uptake per species (amount/min)
    km_asbt::T                                # ASBT Michaelis constant (amount)
    gut_biotransformation::SMatrix{N_SPECIES, N_SPECIES, T, N_SPECIES^2}  # per unit bacterial activity
    # Plasma
    V_P::T                                    # plasma volume (L)
    Q_L::T                                    # hepatic plasma flow (L/min)
    # Postprandial reflex (β, δ), see `gi_reflex`
    gb_reflex::NTuple{2, T}
    gut_reflex::NTuple{2, T}
    meal_times::Vector{Float64}
    # Where it all came from
    parameters::P
    constants::C
end

function BileAcidModel(
        p::Parameters = Parameters(), c::Constants = Constants();
        meal_times::AbstractVector{<:Real} = Float64[]
    )
    T = promote_type(eltype(p), eltype(c))
    issorted(meal_times) || throw(ArgumentError("meal_times must be sorted"))

    # De novo synthesis in the liver: CA, other BA as given, CDCA makes up the remainder.
    # Newly synthesised BAs are conjugated with taurine or glycine.
    composition = by_species(
        T;
        CA = p.frac_CA, DCA = c.frac_DCA, UDCA = c.frac_UDCA, LCA = c.frac_LCA,
        LCAs = c.frac_LCAs, o = c.frac_o,
        CDCA = 1 - (p.frac_CA + c.frac_DCA + c.frac_UDCA + c.frac_LCA + c.frac_LCAs + c.frac_o),
    )
    conjugation = SVector{N_FORMS, T}(0, p.frac_tu, 1 - p.frac_tu)
    synthesis = composition .* conjugation' .* p.k_u

    # Hepatic extraction depends on species; unconjugated BAs are extracted less efficiently
    extraction = by_species(
        T;
        CA = p.frac_li_c_1, CDCA = p.frac_li_c_2, DCA = p.frac_li_c_2, UDCA = p.frac_li_c_2,
        LCA = p.frac_li_c_3, LCAs = p.frac_li_s, o = p.frac_li_c_3,
    )
    hepatic_extraction = hcat(p.frac_li_u * extraction, extraction, extraction)

    liver_biotransformation = conversion_matrix(T, ((:LCA, :LCAs, p.k_tr_6),))
    gut_biotransformation = conversion_matrix(
        T, (
            (:CA, :DCA, p.k_tr_1),
            (:CDCA, :UDCA, p.k_tr_2),
            (:CDCA, :LCA, p.k_tr_3),
            (:UDCA, :LCA, p.k_tr_4),
            (:LCAs, :LCA, p.k_tr_5),
            (:CA, :o, p.k_tr_7),
            (:CDCA, :o, p.k_tr_7),
            (:DCA, :o, p.k_tr_7),
            (:UDCA, :o, p.k_tr_7),
            (:LCA, :o, p.k_tr_7),
            (:LCAs, :o, p.k_tr_7),
        )
    )

    # Intestinal tract: small-intestinal segments followed by colon segments.
    # Bacterial biotransformation is only modelled in the colon.
    small_intestine = ntuple(N_SI) do k
        GutSegment{T}(;
            k_transit = k in TERMINAL_ILEUM ? p.k_si_2 : p.k_si_1,
            bacteria = k in PROXIMAL_SI ? c.k_bact_psi : p.k_bact_dsi,
            k_passive = SVector(p.k_xi_up_si, c.k_xi_tp_si, c.k_xi_gp_si),
            asbt = k in ASBT_SEGMENTS,
            biotransformation = false,
        )
    end
    colon = ntuple(N_CO) do _
        GutSegment{T}(;
            k_transit = p.k_co_1,
            bacteria = c.k_bact_co,
            k_passive = SVector(p.k_xi_up_co, c.k_xi_tp_co, c.k_xi_gp_co),
            asbt = false,
            biotransformation = true,
        )
    end

    vmax_asbt = by_species(
        T;
        CA = p.vmax_asbt, CDCA = p.vmax_asbt, DCA = p.vmax_asbt, UDCA = p.vmax_asbt,
        LCA = p.vmax_asbt, LCAs = p.vmax_asbt_s, o = p.vmax_asbt,
    )

    return from_fields(
        BileAcidModel{T, typeof(p), typeof(c)};
        synthesis, hepatic_extraction, liver_biotransformation,
        frac_tauro = p.frac_tu, k_secretion = p.k_xl, k_spill = c.k_pl,
        frac_gb = p.frac_gb, k_gb_emptying = p.k_xg,
        segments = SVector(small_intestine..., colon...),
        k_deconjugation = SVector(0, p.k_ut, p.k_ug),
        vmax_asbt, km_asbt = p.km_asbt, gut_biotransformation,
        V_P = c.V_P, Q_L = c.Q_L,
        gb_reflex = (p.GIr_beta_GB, p.GIr_delta_GB),
        gut_reflex = (p.GIr_beta_SI, p.GIr_delta_SI),
        meal_times = collect(Float64, meal_times),
        parameters = p, constants = c,
    )
end

"""
    from_fields(S; fields...)

Construct an `S` from keyword arguments named after its fields, given in any order.
Missing or unknown fields are an error.
"""
function from_fields(::Type{S}; fields...) where {S}
    names = fieldnames(S)
    missing_fields, unknown = setdiff(names, keys(fields)), setdiff(keys(fields), names)
    isempty(missing_fields) && isempty(unknown) ||
        throw(ArgumentError("$S: missing fields $missing_fields, unknown fields $unknown"))
    return S(map(name -> fields[name], names)...)
end

"Copy of model `m` with different meal times."
with_meals(m::BileAcidModel, meal_times) = BileAcidModel(m.parameters, m.constants; meal_times)

"Species-indexed vector with the given entries by name; unnamed species are zero."
function by_species(::Type{T}; values...) where {T}
    v = zeros(MVector{N_SPECIES, T})
    for (name, value) in values
        v[species_index(name)] = value
    end
    return SVector(v)
end

"""
    conversion_matrix(T, reactions)

Matrix `M` such that `M * x` is the rate of change of species amounts `x` due to the
first-order conversions `reactions`, given as `(from, to, rate)`.
"""
function conversion_matrix(::Type{T}, reactions) where {T}
    M = zeros(MMatrix{N_SPECIES, N_SPECIES, T})
    for (from, to, rate) in reactions
        i, j = species_index(from), species_index(to)
        M[j, i] += rate
        M[i, i] -= rate
    end
    return SMatrix(M)
end
