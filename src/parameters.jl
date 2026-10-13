# Model parameters (calibrated) and physiological constants (fixed).
#
# Names and ordering follow the MATLAB implementation (`BA_human_model_v9_to_fit_par.m` and
# `BA_human_model_v9_to_fit_con.m`); defaults are the calibrated values `result.p_opt` and
# `result.model_info.c_vec` from `Results/refined_result_MST.mat`.
#
# Units: time in minutes, volumes in litres. Plasma states are concentrations (amount per
# litre of plasma), all other states are amounts.

"""
    Parameters(; kwargs...)
    Parameters(p::Parameters; kwargs...)   # copy of `p` with some values replaced
    Parameters(v::AbstractVector)          # from the 33-element MATLAB parameter vector

Calibrated model parameters. Defaults are the calibrated values of Sips et al.
"""
Base.@kwdef struct Parameters{T <: Real}
    # Postprandial (gastro-intestinal) reflex, see `gi_reflex`
    GIr_beta_GB::T = 70.274395774668108       # time to peak of gallbladder response (min)
    GIr_beta_SI::T = 10.267628634698122       # time to peak of intestinal transit response (min)
    GIr_delta_GB::T = 44.346484115570696      # peak relative increase of gallbladder emptying (-)
    GIr_delta_SI::T = 2.9317984990954673      # peak relative increase of intestinal transit (-)
    # Bile secretion and gallbladder
    k_xl::T = 0.11753367389353529             # hepatic secretion into bile (1/min)
    frac_gb::T = 0.94811260281437848          # fraction of secreted bile stored in the gallbladder (-)
    k_xg::T = 0.0028652013849729966           # basal gallbladder emptying (1/min)
    # Intestinal uptake
    k_xi_up_si::T = 9.0395651931279392e-5     # passive uptake of unconjugated BA, small intestine (1/min)
    k_xi_up_co::T = 0.00019177103357983748    # passive uptake of unconjugated BA, colon (1/min)
    vmax_asbt_s::T = 0.010022096708157457     # maximal ASBT uptake rate of sulfated LCA (amount/min)
    vmax_asbt::T = 439.04707353910908         # maximal ASBT uptake rate of other BA (amount/min)
    km_asbt::T = 9661.5151069162675           # ASBT Michaelis constant (total BA amount in segment)
    # First-pass hepatic extraction (fraction extracted per pass)
    frac_li_u::T = 0.62500066325501069        # unconjugated relative to conjugated (-)
    frac_li_c_1::T = 0.95208405739114754      # conjugated CA (-)
    frac_li_c_2::T = 0.83890649438764708      # conjugated CDCA, DCA, UDCA (-)
    frac_li_c_3::T = 0.77518344092211489      # conjugated LCA and other BA (-)
    frac_li_s::T = 0.043525697636848319       # conjugated sulfated LCA (-)
    # Hepatic conjugation
    frac_tu::T = 0.3679586755580963           # fraction conjugated with taurine; the rest with glycine (-)
    # Bacterial deconjugation (rates per unit bacterial activity)
    k_ut::T = 0.019896305814298521            # tauro-conjugated BA (1/min)
    k_ug::T = 0.0025174644227249307           # glyco-conjugated BA (1/min)
    k_bact_dsi::T = 0.0768606985390793        # bacterial activity in the distal small intestine (-)
    # Intestinal transit
    k_si_1::T = 0.15217525579356136           # basal transit, small-intestinal segments 1-8 (1/min)
    k_si_2::T = 0.010429159954816367          # basal transit, small-intestinal segments 9-10 (1/min)
    k_co_1::T = 0.0021237149095303305         # basal transit, colon segments (1/min)
    # De novo synthesis
    k_u::T = 0.82363442312373258              # total hepatic BA synthesis (amount/min)
    frac_CA::T = 0.53005475636372823          # fraction of synthesis that is CA (-)
    # Biotransformation rates (colon: per unit bacterial activity; k_tr_6: liver)
    k_tr_1::T = 0.0026788415925414811         # CA   -> DCA   (7α-dehydroxylation) (1/min)
    k_tr_2::T = 0.0020790155306778467         # CDCA -> UDCA  (7-epimerisation) (1/min)
    k_tr_3::T = 1.7277798801154969e-7         # CDCA -> LCA   (7α-dehydroxylation) (1/min)
    k_tr_4::T = 0.0040848035348395465         # UDCA -> LCA   (7β-dehydroxylation) (1/min)
    k_tr_5::T = 0.00072634610015056265        # LCAs -> LCA   (desulfation) (1/min)
    k_tr_6::T = 0.066956822470394528          # LCA  -> LCAs  (hepatic sulfation) (1/min)
    k_tr_7::T = 8.3122533183443829e-5         # BA   -> other BA (1/min)
end

"""
    Constants(; kwargs...)
    Constants(c::Constants; kwargs...)     # copy of `c` with some values replaced
    Constants(c_vec::AbstractVector)       # from the 128-element MATLAB constants vector

Physiological constants that are not calibrated. Only the constants that enter the model
equations are included; the MATLAB constants vector contains many more entries that the
equations never use.
"""
Base.@kwdef struct Constants{T <: Real}
    V_P::T = 3.15                     # plasma volume (L)
    Q_L::T = 0.82500000000000007      # hepatic plasma flow (L/min)
    k_pl::T = 0.0                     # spill-over from liver into plasma (1/min)
    k_xi_tp_si::T = 0.0               # passive uptake of tauro-conjugated BA, small intestine (1/min)
    k_xi_tp_co::T = 0.0               # passive uptake of tauro-conjugated BA, colon (1/min)
    k_xi_gp_si::T = 0.0               # passive uptake of glyco-conjugated BA, small intestine (1/min)
    k_xi_gp_co::T = 0.0               # passive uptake of glyco-conjugated BA, colon (1/min)
    k_bact_psi::T = 0.0               # bacterial activity, proximal small intestine (-)
    k_bact_co::T = 1.0                # bacterial activity, colon (-)
    frac_DCA::T = 0.0                 # fraction of synthesis that is DCA (-)
    frac_UDCA::T = 0.0                # fraction of synthesis that is UDCA (-)
    frac_LCA::T = 0.0                 # fraction of synthesis that is LCA (-)
    frac_LCAs::T = 0.0                # fraction of synthesis that is sulfated LCA (-)
    frac_o::T = 0.0                   # fraction of synthesis that is other BA (-)
end

"Positions of the `Constants` fields in the MATLAB constants vector `c_vec`."
const MATLAB_CONSTANT_INDEX = (
    V_P = 8, Q_L = 9, k_pl = 32,
    k_xi_tp_si = 2, k_xi_tp_co = 3, k_xi_gp_si = 4, k_xi_gp_co = 5,
    k_bact_psi = 14, k_bact_co = 15,
    frac_DCA = 20, frac_UDCA = 24, frac_LCA = 25, frac_LCAs = 26, frac_o = 30,
)

for S in (:Parameters, :Constants)
    @eval begin
        # Promote mixed argument types, e.g. `Parameters(frac_gb = 0)`
        $S(args::Real...) = $S(promote(args...)...)
        # Copy with replaced fields
        $S(x::$S; kwargs...) = $S(; _fields(x)..., kwargs...)
        Base.eltype(::Type{<:$S{T}}) where {T} = T
        # Values in field order (for Parameters: the MATLAB parameter vector)
        Base.Vector(x::$S) = [getfield(x, i) for i in 1:fieldcount(typeof(x))]
        Base.show(io::IO, x::$S) = print(io, $(string(S)), _fields(x))
    end
end

_fields(x) = NamedTuple{fieldnames(typeof(x))}(ntuple(i -> getfield(x, i), fieldcount(typeof(x))))

function Parameters(v::AbstractVector{<:Real})
    length(v) == fieldcount(Parameters) ||
        throw(DimensionMismatch("expected $(fieldcount(Parameters)) parameters, got $(length(v))"))
    return Parameters(v...)
end

function Constants(c_vec::AbstractVector{<:Real})
    length(c_vec) == 128 ||
        throw(DimensionMismatch("expected the 128-element MATLAB constants vector, got $(length(c_vec))"))
    return Constants(; map(i -> c_vec[i], MATLAB_CONSTANT_INDEX)...)
end
