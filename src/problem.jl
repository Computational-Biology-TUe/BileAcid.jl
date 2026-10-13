# Setting up and solving the ODE problem.

"""
    jacobian_prototype(m::BileAcidModel)

Sparse matrix with the structural non-zeros of the Jacobian of `rhs!`. It is detected by
tracing `rhs!`, so it follows any change to the model equations automatically, and it
does not depend on parameter values.
"""
function jacobian_prototype(m::BileAcidModel)
    pattern = SparseConnectivityTracer.jacobian_sparsity(
        (du, u) -> rhs!(du, u, m, 0.0), zeros(N_STATES), zeros(N_STATES),
        TracerSparsityDetector()
    )
    return SparseMatrixCSC{Float64, Int}(pattern)
end

"""
    meal_callback(m::BileAcidModel)

Restarts the integration at every meal, as the MATLAB implementation does by starting a
new solve per meal. Without it, a solver arriving at a meal with the large step size of
the preceding quiet period can step straight over the postprandial response; multistep
methods (FBDF, QNDF), which only evaluate the right-hand side at the end of a step, do
so without noticing. The callback discards the multistep history and reduces the next
step to a tenth of the fastest reflex time-to-peak.
"""
function meal_callback(m::BileAcidModel)
    dt_meal = min(m.gb_reflex[1], m.gut_reflex[1]) / 10
    is_meal(_, t, _) = insorted(t, m.meal_times)
    function restart!(integrator)
        set_proposed_dt!(integrator, min(SciMLBase.get_proposed_dt(integrator), dt_meal))
        derivative_discontinuity!(integrator, true)
        return nothing
    end
    return DiscreteCallback(is_meal, restart!; save_positions = (false, false))
end

"""
    ODEProblem(m::BileAcidModel, u0, tspan; kwargs...)

ODE problem for model `m`. Steps end exactly at meal times, where `meal_callback`
restarts the integration. Implicit solvers exploit the sparsity of the Jacobian (sparse
automatic differentiation and sparse LU factorisation).

Passing `callback` to `solve` replaces the meal callback; use
`CallbackSet(meal_callback(m), your_callback)` to keep it.
"""
function SciMLBase.ODEProblem(m::BileAcidModel, u0, tspan; kwargs...)
    f = ODEFunction{true, SciMLBase.FullSpecialize}(rhs!; jac_prototype = jacobian_prototype(m))
    tstops = filter(t -> tspan[1] < t < tspan[2], m.meal_times)
    return ODEProblem(f, u0, tspan, m; tstops, callback = meal_callback(m), kwargs...)
end

"""
    simulate(p = Parameters(), c = Constants(), schedule = MealSchedule(); u0, solver, kwargs...)

Simulate the model with parameters `p` and constants `c` under the meal `schedule`,
starting from `u0` (default: no bile acids anywhere). Remaining keyword arguments are
passed to `solve`, e.g. `saveat`. Returns the `ODESolution`; see `compartment` and
`state_index` for extracting results.

The default solver and tolerances reproduce a tight-tolerance MATLAB solution to about
1e-5 (relative to the range of each state).
"""
function simulate(
        p::Parameters = Parameters(), c::Constants = Constants(), schedule::MealSchedule = MealSchedule();
        u0 = zeros(promote_type(eltype(p), eltype(c)), N_STATES),
        solver = QNDF(), reltol = 1.0e-6, abstol = 1.0e-8, kwargs...
    )
    m = BileAcidModel(p, c; schedule.meal_times)
    return solve(ODEProblem(m, u0, schedule.tspan), solver; reltol, abstol, kwargs...)
end
