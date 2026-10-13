# Julia version of BA_model_package/BA_model_EXAMPLE.m.
#
# Run from the repository root with:  julia --project=. examples/example.jl

using BileAcids
using Printf

# Total bile acid concentration in plasma, and the bile acid pool (everything except faeces)
plasma_total(u) = sum(compartment(u, :pl))
pool(u) = sum(sum(compartment(u, c)) for c in COMPARTMENTS if c ∉ (:fa, :pl))

function describe_last_day(label, sol, schedule)
    day = last_day(schedule)
    on_day = findall(t -> day[1] <= t <= day[2], sol.t)
    plasma = plasma_total.(sol.u[on_day])
    @printf("%-32s pool %6.0f   plasma BA %.2f-%.2f (mean %.2f)\n",
        label, pool(sol.u[end]), minimum(plasma), maximum(plasma), sum(plasma) / length(plasma))
end

# ------------------------------------------------------------------------------------------
# Example 1: long-term change of a parameter value
# ------------------------------------------------------------------------------------------
# Gallbladder dysfunction (all bile flows directly to the duodenum), with meals 4 h, 4 h
# and 16 h apart. Each simulation covers 20 days without meals and 35 days with meals.

standard = MealSchedule()                    # gaps of 6, 6 and 12 hours
schedule1 = MealSchedule(gaps = (4, 4, 16))

healthy = simulate(Parameters(), Constants(), standard; saveat = 10.0)
example1 = simulate(Parameters(frac_gb = 0), Constants(), schedule1; saveat = 10.0)

println("Example 1: final day")
describe_last_day("healthy", healthy, standard)
describe_last_day("no gallbladder storage", example1, schedule1)

# ------------------------------------------------------------------------------------------
# Example 2: short-term change of a parameter value
# ------------------------------------------------------------------------------------------
# Start from the end of the healthy simulation and simulate one more day, without the
# postprandial acceleration of intestinal transit.

one_day = MealSchedule(fasting_days = 0, meal_days = 1)
first_meal = 0.0:1.0:360.0
reference = simulate(Parameters(), Constants(), one_day; u0 = healthy.u[end], saveat = first_meal)
example2 = simulate(Parameters(GIr_delta_SI = 0), Constants(), one_day; u0 = healthy.u[end], saveat = first_meal)

println("\nExample 2: plasma bile acids after the first meal")
for (label, sol) in (("healthy", reference), ("no transit response", example2))
    plasma = plasma_total.(sol.u)
    peak = argmax(plasma)
    @printf("%-32s peak %.2f at %3.0f min, fasting %.2f\n", label, plasma[peak], sol.t[peak], plasma[1])
end

# Time courses of individual states, e.g. for plotting:
#   sol[state_index(:pl, :g, :CDCA), :]     glyco-CDCA in plasma
#   compartment(sol.u[end], :li)            liver contents (species × form matrix)
