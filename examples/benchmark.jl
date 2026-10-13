# Timing of the right-hand side and of the standard 55-day simulation.
#
# Run from the repository root with:  julia --project=. examples/benchmark.jl

using BileAcids
using Printf

best_time(f; repeats) = minimum(@elapsed(f()) for _ in 1:repeats)

function time_rhs(m; n = 100_000)
    u, du = rand(N_STATES), zeros(N_STATES)
    rhs!(du, u, m, 30.0)
    return best_time(() -> (for _ in 1:n; rhs!(du, u, m, 30.0); end); repeats = 5) / n
end
t_rhs = time_rhs(BileAcidModel(; meal_times = MealSchedule().meal_times))
@printf("rhs!:                                  %6.0f ns\n", 1.0e9 * t_rhs)

for (label, tolerances) in (
        ("simulate(), default tolerances:", (;)),
        ("simulate(), reltol = 1e-4:", (; reltol = 1.0e-4, abstol = 1.0e-6)),
    )
    simulate(; saveat = 60.0, tolerances...)
    @printf("%-38s %6.2f s\n", label, best_time(() -> simulate(; saveat = 60.0, tolerances...); repeats = 3))
end

matlab = only(readlines(joinpath(@__DIR__, "..", "test", "reference", "sim_matlab_time.csv")))
@printf("MATLAB simulate_human_day (ode15s):    %6.2f s\n", parse(Float64, matlab))
