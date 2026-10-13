# Compartment equations.
#
# Each compartment holds a species × form matrix of bile acid amounts (`BAMatrix`). The
# functions below give the rate of change of one compartment from its contents and the
# fluxes it exchanges with other compartments; `rhs!` connects them.

"""
    liver(x, extracted, secretion, spill, m) -> dx

Liver: de novo synthesis, uptake of the bile acids `extracted` from portal and systemic
blood, and hepatic biotransformation, minus `secretion` into bile and `spill` into plasma.
Unconjugated bile acids are conjugated with taurine or glycine on uptake.
"""
@inline function liver(x::BAMatrix, extracted, secretion, spill, m::BileAcidModel)
    dx = m.synthesis + conjugate(extracted, m.frac_tauro) +
        m.liver_biotransformation * x - secretion - spill
    return conjugated_only(dx)
end

"""
    gallbladder(x, inflow, m, reflex) -> (dx, emptying)

Gallbladder: stores the `inflow` of secreted bile and empties into the duodenum, much
faster after a meal (`reflex`).
"""
@inline function gallbladder(x::BAMatrix, inflow, m::BileAcidModel, reflex)
    emptying = (m.k_gb_emptying * reflex) * x
    return conjugated_only(inflow - emptying), emptying
end

"""
    gut_segment(x, inflow, seg::GutSegment, m, reflex) -> (dx, outflow, uptake)

One segment of the intestinal tract with contents `x`, receiving `inflow` from upstream
(bile, for the duodenum). Bile acids leave by transit to the next segment (`outflow`,
accelerated after a meal by `reflex`) and by uptake into portal blood (`uptake`); bacteria
deconjugate and biotransform them within the segment.
"""
@inline function gut_segment(x::BAMatrix, inflow, seg::GutSegment, m::BileAcidModel, reflex)
    outflow = (seg.k_transit * reflex) * x
    uptake = x .* seg.k_passive'
    if seg.asbt
        uptake += active_uptake(x, m)
    end
    dx = inflow - outflow - uptake + deconjugation(x, seg.bacteria * m.k_deconjugation)
    if seg.biotransformation
        dx += unconjugated(seg.bacteria * (m.gut_biotransformation * x[:, UNCONJ]))
    end
    return dx, outflow, uptake
end

"""
    plasma(absorbed, extracted, spill, m) -> dx

Plasma (concentrations): bile acids `absorbed` from the gut that are not `extracted` by
the liver on first pass, minus hepatic clearance of plasma bile acids (also part of
`extracted`), plus `spill` from the liver.
"""
@inline plasma(absorbed, extracted, spill, m::BileAcidModel) = (absorbed - extracted + spill) / m.V_P

# --- Processes ---------------------------------------------------------------------------

"""
Saturable uptake by the apical sodium-dependent bile acid transporter, competitively
inhibited by all bile acids in the segment.
"""
@inline active_uptake(x::BAMatrix, m::BileAcidModel) = m.vmax_asbt .* x ./ (sum(x) + m.km_asbt)

"""
Bacterial deconjugation, at rate `k` per conjugation form: conjugated bile acids are
converted to the unconjugated form of the same species.
"""
@inline function deconjugation(x::BAMatrix, k::SVector{N_FORMS})
    flux = x .* k'
    return hcat(flux[:, TCONJ] + flux[:, GCONJ], -flux[:, TCONJ], -flux[:, GCONJ])
end

"""
Hepatic conjugation: unconjugated bile acids are conjugated with taurine (fraction
`frac_tauro`) or glycine.
"""
@inline function conjugate(x::BAMatrix, frac_tauro)
    free = x[:, UNCONJ]
    return hcat(zero(free), x[:, TCONJ] + frac_tauro * free, x[:, GCONJ] + (1 - frac_tauro) * free)
end

"Matrix with species amounts `v` in the unconjugated form only."
@inline unconjugated(v::SVector{N_SPECIES}) = hcat(v, zero(v), zero(v))

"""
Liver and gallbladder only contain conjugated bile acids: their unconjugated states are
kept constant (at zero), as in the MATLAB model.
"""
@inline conjugated_only(dx::BAMatrix) = hcat(zero(dx[:, UNCONJ]), dx[:, TCONJ], dx[:, GCONJ])
