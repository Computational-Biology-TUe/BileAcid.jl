using BileAcids
using BileAcids: BLOCK, PLASMA, block_offset
using DelimitedFiles
using ForwardDiff
using Test

# Reference data from the MATLAB implementation, see reference/generate_reference.m
readref(name) = readdlm(joinpath(@__DIR__, "reference", name), ',')

"Total amount of bile acids in the body and faeces (plasma states are concentrations)."
function total_amount(u, V_P)
    plasma = block_offset(PLASMA) .+ (1:BLOCK)
    return sum(u[1:(first(plasma) - 1)]) + V_P * sum(u[plasma])
end

"Largest error of each state relative to its own range, over all states."
function max_relative_error(X, Xref)
    scale = maximum(abs, Xref; dims = 1)
    keep = vec(scale) .> 1.0e-12
    return maximum((maximum(abs, X .- Xref; dims = 1) ./ scale)[keep])
end

@testset "BileAcids" begin
    # skipped
    # @testset "state layout matches MATLAB" begin
    #     ode = joinpath(@__DIR__, "..", "BA_model_package", "Human_MST", "BA_human_model_v9_to_fit_full_ode.m")
    #     matlab_names = [m[1] for m in eachmatch(r"^(\w+) = x\(\d+\);"m, read(ode, String))]
    #     @test state_names() == matlab_names
    #     @test state_index(:si3, :t, :CA) == findfirst(==("si3_t_CA"), matlab_names)
    #     @test state_index(:cum_excreted_BA) == N_STATES
    #     u = collect(1.0:N_STATES)
    #     @test compartment(u, :pl)[2, 3] == u[state_index(:pl, :g, :CDCA)]
    # end

    @testset "parameters and constants" begin
        p = Parameters()
        @test Parameters(Vector(p)) == p
        @test Parameters(p; frac_gb = 0).frac_gb === 0.0
        @test Parameters(p; frac_gb = 0).k_xl === p.k_xl
        # The first reference case uses the calibrated constants vector
        @test Constants(readref("rhs_c.csv")[1, :]) == Constants()
        @test Parameters(readref("rhs_p.csv")[1, :]) == p
    end

    @testset "postprandial reflex and meal schedule" begin
        @test gi_reflex(10.0, 2.0, 10.0) ≈ 3.0      # peak 1 + δ at t = β
        @test gi_reflex(10.0, 2.0, 0.0) == 1.0
        @test gi_reflex(10.0, 2.0, -5.0) == 1.0     # before the meal
        s = MealSchedule(fasting_days = 1, meal_days = 2, gaps = (6, 6, 12))
        @test s.meal_times == [0, 1440, 1800, 2160, 2880, 3240, 3600]
        @test s.tspan == (0.0, 4320.0)
        @test last_day(s) == (2880.0, 4320.0)
        @test MealSchedule(fasting_days = 0, meal_days = 1).meal_times == [0, 360, 720]
        @test MealSchedule(meal_days = 0).meal_times == [0]
    end

    # @testset "right-hand side matches MATLAB" begin
    #     X, T, P, C, DX = readref.(("rhs_x.csv", "rhs_t.csv", "rhs_p.csv", "rhs_c.csv", "rhs_dx.csv"))
    #     du = zeros(N_STATES)
    #     # Cases 1-20: calibrated parameters. Cases 21-40: perturbed parameters with all
    #     # constants that are zero in the calibrated set switched on.
    #     for i in axes(X, 1)
    #         m = BileAcidModel(Parameters(P[i, :]), Constants(C[i, :]); meal_times = [T[i, 2]])
    #         rhs!(du, X[i, :], m, T[i, 1])
    #         tol = 1.0e-12 .* (abs.(DX[i, :]) .+ maximum(abs, DX[i, :]))
    #         @test all(abs.(du .- DX[i, :]) .<= tol)
    #     end
    # end

    @testset "rhs! does not allocate" begin
        m = BileAcidModel(; meal_times = [0.0])
        u, du = rand(N_STATES), zeros(N_STATES)
        rhs!(du, u, m, 30.0)
        @test (@allocated rhs!(du, u, m, 30.0)) == 0
    end

    @testset "Jacobian sparsity pattern" begin
        # All terms switched on, so that the true Jacobian has every structural non-zero
        c = Constants(k_pl = 0.01, k_bact_psi = 0.05, k_xi_tp_si = 1.0e-4, k_xi_gp_si = 1.0e-4,
            k_xi_tp_co = 1.0e-4, k_xi_gp_co = 1.0e-4, frac_DCA = 0.01)
        m = BileAcidModel(Parameters(), c; meal_times = [0.0])
        pattern = jacobian_prototype(m)
        @test pattern == jacobian_prototype(BileAcidModel())   # independent of values
        J = ForwardDiff.jacobian((du, u) -> rhs!(du, u, m, 30.0), zeros(N_STATES), rand(N_STATES))
        @test all(!iszero(pattern[i]) for i in findall(!iszero, J))
    end

    @testset "derivatives with respect to parameters" begin
        u = readref("rhs_x.csv")[1, :]
        function f(k_xl)
            du = zeros(typeof(k_xl), N_STATES)
            rhs!(du, u, BileAcidModel(Parameters(Parameters(); k_xl); meal_times = [0.0]), 30.0)
            return du
        end
        h = 1.0e-6
        @test ForwardDiff.derivative(f, 0.1) ≈ (f(0.1 + h) - f(0.1 - h)) / 2h rtol = 1.0e-6
    end

    @testset "simulation matches MATLAB" begin
        t, X = vec(readref("sim_tight_t.csv")), readref("sim_tight_x.csv")
        sol = simulate(; reltol = 1.0e-8, abstol = 1.0e-10, saveat = t)
        @test sol.t == t
        @test max_relative_error(stack(sol.u; dims = 1), X) < 1.0e-5
        sol = simulate(; saveat = t)   # default tolerances
        @test max_relative_error(stack(sol.u; dims = 1), X) < 1.0e-4
    end

    @testset "mass balance" begin
        # Bile acids are only created by synthesis (k_u in total) and never destroyed, so
        # the total amount grows at rate k_u (given that liver and gallbladder contain no
        # unconjugated bile acids, as in any simulation).
        p = Parameters(Parameters(); frac_gb = 0.5, k_u = 1.3)
        c = Constants(k_pl = 0.01, k_bact_psi = 0.05, k_xi_tp_si = 1.0e-4, frac_DCA = 0.05)
        m = BileAcidModel(p, c; meal_times = [0.0])
        u, du = 100 .* rand(N_STATES), zeros(N_STATES)
        for s in SPECIES, organ in (:li, :gb)
            u[state_index(organ, :u, s)] = 0
        end
        rhs!(du, u, m, 20.0)
        @test total_amount(du, c.V_P) ≈ p.k_u rtol = 1.0e-12
        # ... and in a simulation from zero, up to the solver tolerance
        sol = simulate(p, c, MealSchedule(fasting_days = 2, meal_days = 3))
        @test total_amount(sol.u[end], c.V_P) ≈ p.k_u * sol.t[end] rtol = 1.0e-5
        # Everything that reaches the faeces has been counted as excreted
        @test sum(compartment(sol.u[end], :fa)) ≈ sol.u[end][state_index(:cum_excreted_BA)] rtol = 1.0e-12
    end
end
