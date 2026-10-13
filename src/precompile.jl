# Compile the common simulation path when the package is precompiled, so that the first
# `simulate` call in a session is fast.
@setup_workload begin
    schedule = MealSchedule(fasting_days = 1, meal_days = 1)
    @compile_workload begin
        sol = simulate(Parameters(), Constants(), schedule; saveat = 60.0)
        compartment(sol.u[end], :pl)
    end
end
