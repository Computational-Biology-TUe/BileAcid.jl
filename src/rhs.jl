# Right-hand side of the ODE system: connects the compartments of the enterohepatic
# circulation.
#
#   synthesis
#       │
#       ▼    secretion                 emptying
#     liver ─────┬────> gallbladder ──────────┐
#     ▲   │      └────────────────────────────┴──> si1 → … → si10 → co1 → … → co5 → faeces
#     │   │ spill                                   └─────── uptake ───────┘
#     │   ▼                                                    │
#     │ plasma ◄─── not extracted ─── portal blood ◄───────────┘
#     │   │                               │
#     └───┴──────── hepatic extraction ◄──┘

"""
    rhs!(du, u, m::BileAcidModel, t)

Time derivative `du` of state `u` (see `state_names()` for the layout) at time `t` (min).
"""
function rhs!(du, u, m::BileAcidModel, t)
    gb_reflex, gut_reflex = postprandial_reflex(m, t)

    li = compartment(u, LIVER)
    gb = compartment(u, GALLBLADDER)
    pl = compartment(u, PLASMA)

    # Bile: hepatic secretion is partly stored in the gallbladder, partly flows directly
    # into the duodenum.
    secretion = m.k_secretion * li
    dgb, emptying = gallbladder(gb, m.frac_gb * secretion, m, gb_reflex)
    bile = (1 - m.frac_gb) * secretion + emptying

    # Intestinal tract, from duodenum to faeces
    absorbed, colon_entry, faecal = gut!(du, u, bile, m, gut_reflex)

    # The liver extracts bile acids from portal blood (absorbed from the gut) and from
    # systemic plasma; what escapes first-pass extraction enters plasma.
    extracted = m.hepatic_extraction .* (absorbed + m.Q_L * pl)
    spill = m.k_spill * li

    set_compartment!(du, LIVER, liver(li, extracted, secretion, spill, m))
    set_compartment!(du, GALLBLADDER, dgb)
    set_compartment!(du, FAECES, faecal)
    set_compartment!(du, PLASMA, plasma(absorbed, extracted, spill, m))

    du[CUM_SECRETED] = sum(bile)
    du[CUM_COLON] = sum(colon_entry)
    du[CUM_EXCRETED] = sum(faecal)
    return nothing
end

"""
    gut!(du, u, bile, m, reflex) -> (absorbed, colon_entry, faecal)

Intestinal tract: contents move segment by segment, starting with `bile` entering the
duodenum. Writes the rate of change of each segment into `du`; returns the total uptake
into portal blood, the flux from the small intestine into the colon, and the flux into
the faeces.
"""
@inline function gut!(du, u, bile, m::BileAcidModel, reflex)
    inflow = bile
    absorbed = zero(bile)
    colon_entry = zero(bile)
    for k in 1:N_GUT
        dx, outflow, uptake = gut_segment(compartment(u, GUT[k]), inflow, m.segments[k], m, reflex)
        set_compartment!(du, GUT[k], dx)
        absorbed += uptake
        if k == N_SI
            colon_entry = outflow
        end
        inflow = outflow
    end
    return absorbed, colon_entry, inflow
end
