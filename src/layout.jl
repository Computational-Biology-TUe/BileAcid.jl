# State layout: bile acid species, conjugation forms and compartments, and how they map
# onto the state vector.
#
# Every compartment holds one amount per (species, form) pair, stored as a contiguous block
# of N_SPECIES × N_FORMS states: species vary fastest, then form. Blocks are ordered as in
# COMPARTMENTS and followed by the cumulative flux counters. This is exactly the ordering
# of the MATLAB model, so state vectors can be exchanged between the two implementations.

"Bile acid species (`o` collects all other bile acids)."
const SPECIES = (:CA, :CDCA, :DCA, :UDCA, :LCA, :LCAs, :o)

"Conjugation forms: unconjugated, taurine-conjugated, glycine-conjugated."
const FORMS = (:u, :t, :g)

const N_SPECIES = length(SPECIES)
const N_FORMS = length(FORMS)
const BLOCK = N_SPECIES * N_FORMS

"Number of small-intestinal segments (`si1` = duodenum, …, `si$N_SI` = terminal ileum)."
const N_SI = 10
"Number of colon segments."
const N_CO = 5
"The gut is a chain of `N_SI` small-intestinal segments followed by `N_CO` colon segments."
const N_GUT = N_SI + N_CO

"Names of the gut segments, from duodenum to distal colon."
const GUT_SEGMENTS = (ntuple(k -> Symbol(:si, k), N_SI)..., ntuple(k -> Symbol(:co, k), N_CO)...)

"Compartments, in state-vector order: liver, gallbladder, gut segments, faeces, plasma."
const COMPARTMENTS = (:li, :gb, GUT_SEGMENTS..., :fa, :pl)
const N_COMPARTMENTS = length(COMPARTMENTS)

"Cumulative fluxes, stored after the compartment blocks."
const CUMULATIVE = (:cum_secreted_BA, :cum_colon_BA, :cum_excreted_BA)

"Length of the state vector."
const N_STATES = N_COMPARTMENTS * BLOCK + length(CUMULATIVE)

"Amounts of each species (rows) in each conjugation form (columns) within one compartment."
const BAMatrix{T} = SMatrix{N_SPECIES, N_FORMS, T, BLOCK}

# --- Names and indices -------------------------------------------------------------------

function _position(name::Symbol, list, what)
    i = findfirst(==(name), list)
    i === nothing && throw(ArgumentError("unknown $what :$name; expected one of $list"))
    return i
end

compartment_index(name::Symbol) = _position(name, COMPARTMENTS, "compartment")
species_index(name::Symbol) = _position(name, SPECIES, "species")
form_index(name::Symbol) = _position(name, FORMS, "conjugation form")

"""
    state_index(compartment, form, species)
    state_index(cumulative)

Position in the state vector of e.g. `state_index(:pl, :t, :CA)` (taurocholate in plasma)
or `state_index(:cum_excreted_BA)`.
"""
function state_index(c::Symbol, f::Symbol, s::Symbol)
    return block_offset(compartment_index(c)) +
        (form_index(f) - 1) * N_SPECIES + species_index(s)
end
state_index(cumulative::Symbol) =
    N_COMPARTMENTS * BLOCK + _position(cumulative, CUMULATIVE, "cumulative flux")

"""
    state_names()

Names of all states in state-vector order, following the MATLAB convention
`<compartment>_<form>_<species>` (e.g. `"si3_t_CA"`).
"""
function state_names()
    names = String[]
    for c in COMPARTMENTS, f in FORMS, s in SPECIES
        push!(names, "$(c)_$(f)_$(s)")
    end
    append!(names, string.(CUMULATIVE))
    return names
end

# --- Positions used by the right-hand side -----------------------------------------------

const UNCONJ = form_index(:u)
const TCONJ = form_index(:t)
const GCONJ = form_index(:g)

const LIVER = compartment_index(:li)
const GALLBLADDER = compartment_index(:gb)
const GUT = map(compartment_index, GUT_SEGMENTS)
const FAECES = compartment_index(:fa)
const PLASMA = compartment_index(:pl)

const CUM_SECRETED = state_index(:cum_secreted_BA)   # bile entering the duodenum
const CUM_COLON = state_index(:cum_colon_BA)         # chyme entering the colon
const CUM_EXCRETED = state_index(:cum_excreted_BA)   # faecal excretion

# --- Reading and writing compartment blocks ----------------------------------------------

@inline block_offset(c::Integer) = (c - 1) * BLOCK

"""
    compartment(u, c)

Contents of compartment `c` (index or name, e.g. `:pl`) in state vector `u`, as a
species × form `BAMatrix`.
"""
@inline function compartment(u::AbstractVector{T}, c::Integer) where {T}
    o = block_offset(c)
    return BAMatrix{T}(ntuple(i -> @inbounds(u[o + i]), Val(BLOCK)))
end
compartment(u::AbstractVector, name::Symbol) = compartment(u, compartment_index(name))

@inline function set_compartment!(du::AbstractVector, c::Integer, x::BAMatrix)
    o = block_offset(c)
    @inbounds for i in 1:BLOCK
        du[o + i] = x[i]
    end
    return du
end
