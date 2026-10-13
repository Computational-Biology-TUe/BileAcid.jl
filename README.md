# BileAcids.jl

A Julia implementation of the human bile acid model of Sips et al., *In silico analysis
identifies intestinal transit as a key determinant of systemic bile acid metabolism*,
ported from the MATLAB code in [`BA_model_package`](BA_model_package).

Seven bile acid species (CA, CDCA, DCA, UDCA, LCA, sulfated LCA, other), each
unconjugated or conjugated with taurine or glycine, are followed through liver,
gallbladder, ten small-intestinal segments, five colon segments, faeces and plasma.

The right-hand side reproduces the MATLAB model to rounding error (relative error
< 1e-13), and the standard 55-day simulation runs in 0.73 s instead of 43 s.

## Quick start

```julia
using BileAcids

sol = simulate()        # 20 days without meals, then 35 days with meals at 6 h, 6 h, 12 h

# Change parameters, constants or the meal pattern
sol = simulate(Parameters(frac_gb = 0), Constants(), MealSchedule(gaps = (4, 4, 16)))

# Results: amounts per species × form for a compartment, or the trajectory of one state
compartment(sol.u[end], :li)              # liver contents (7 × 3 matrix)
sol[state_index(:pl, :g, :CDCA), :]       # glyco-CDCA in plasma over time
```

`simulate` also accepts `u0` and `solver`; other keyword arguments go to `solve` (`saveat`, `reltol`, ...).
Time is in minutes. Plasma states are concentrations; all other states are amounts.

[`examples/example.jl`](examples/example.jl) is a port of `BA_model_EXAMPLE.m`.

## Code structure

The state vector has the same layout as in MATLAB: 19 compartments of 21 states each
(7 species × 3 forms, species varying fastest), followed by 3 cumulative fluxes. Inside
the model, each compartment is handled as a 7 × 3 static matrix of species × forms.

| File | Contents | Edit this to … |
|---|---|---|
| [`src/layout.jl`](src/layout.jl) | species, conjugation forms, compartments, state indexing | add a compartment or species, change the number of gut segments |
| [`src/parameters.jl`](src/parameters.jl) | `Parameters` (calibrated) and `Constants` (fixed), with defaults and units | add or change a parameter |
| [`src/model.jl`](src/model.jl) | `BileAcidModel`: rates derived from the parameters, properties of each gut segment, biotransformation reactions | change where a process occurs or add a reaction |
| [`src/meals.jl`](src/meals.jl) | postprandial reflex, `MealSchedule` | change the meal response or protocol |
| [`src/compartments.jl`](src/compartments.jl) | equations of each compartment (liver, gallbladder, gut segment, plasma) | change the equations of a compartment |
| [`src/rhs.jl`](src/rhs.jl) | `rhs!`: connects the compartments | add or reroute fluxes between compartments |
| [`src/problem.jl`](src/problem.jl) | `ODEProblem` setup, `simulate` | change solver settings |

`rhs!` reads like a flow diagram of the enterohepatic circulation. The gut is a chain of
segments that share one equation, `gut_segment`. Differences between segments, such as
ASBT uptake in the terminal ileum or biotransformation in the colon, are data in
`GutSegment` rather than separate code.

## Modifying the model

**Parameters.** `Parameters(p; k_xl = 0.1)` copies `p` with some values replaced.
`Parameters(v)` and `Vector(p)` convert from and to the MATLAB parameter vector, and
`Constants(c_vec)` reads the MATLAB constants vector.

**Where processes occur.** See the top of `src/model.jl`, e.g. `ASBT_SEGMENTS = 9:10`,
and the `GutSegment` definitions. Biotransformations are listed by name as
`(from, to, rate)`, e.g. `(:CA, :DCA, p.k_tr_1)`.

**Adding a compartment.** For example, urinary excretion from plasma:

1. `src/layout.jl`: append `:ur` to `COMPARTMENTS` and define `const URINE = compartment_index(:ur)`.
2. `src/parameters.jl`: add `k_renal::T = 0.0  # renal clearance (L/min)` to `Constants`.
3. `src/model.jl`: add a field `k_renal::T` to `BileAcidModel`, and `k_renal = c.k_renal`
   in the constructor call at the end of `BileAcidModel(p, c)`.
4. `src/rhs.jl`:
   ```julia
   renal = m.k_renal * pl
   set_compartment!(du, PLASMA, plasma(absorbed, extracted, spill, m) - renal / m.V_P)
   set_compartment!(du, URINE, renal)
   ```

The state vector length, state names and Jacobian sparsity pattern follow automatically.
This example has been run as described. Appending at the end keeps the original states
in their MATLAB positions. The tests that compare against MATLAB will fail after a
structural change, as they should.

## Solving

`simulate` and `ODEProblem(m::BileAcidModel, u0, tspan)` set up three things:

- **Sparse Jacobian.** The sparsity pattern is detected by tracing `rhs!`
  (`jacobian_prototype`). Implicit solvers use it for sparse automatic differentiation
  and sparse LU, which is 5× faster than dense.
- **Meal restarts.** Steps end exactly at meal times (`tstops`), and `meal_callback`
  restarts the integrator there with a small step. This matters: a BDF solver arriving
  at a meal with the large step of a quiet period can step over the whole postprandial
  response and still report success. Without the callback, FBDF and QNDF gave errors
  of order 100%. The MATLAB code avoids this by starting a new `ode15s` call per meal.
  If you pass your own `callback` to `solve`, use `CallbackSet(meal_callback(m), cb)`.
- **Solver.** QNDF with `reltol = 1e-6, abstol = 1e-8`. FBDF converged poorly on this
  problem; Rodas5P and KenCarp47 work, but are slower.

For repeated solves (e.g. fitting), build the problem yourself; setting it up costs well
under a millisecond:

```julia
using BileAcids, OrdinaryDiffEq
schedule = MealSchedule()
m = BileAcidModel(Parameters(k_xl = 0.1), Constants(); schedule.meal_times)
sol = solve(ODEProblem(m, zeros(N_STATES), schedule.tspan), Rodas5P(); reltol = 1e-6, abstol = 1e-8)
```

The model is generic in its number type, so `ForwardDiff` works through `BileAcidModel`
and `rhs!`.

### Performance

Measured on the development machine (Windows, Julia 1.13, MATLAB R2022b) for the
standard protocol (20 + 35 days). `simulate()` takes 0.73 s with its default settings;
MATLAB's `simulate_human_day` as shipped took 43 s. `rhs!` takes about 0.4 µs and does
not allocate. Thanks to a precompile workload, the first `simulate` in a session is as
fast as later ones (`using BileAcids` takes about 4 s);
[`examples/benchmark.jl`](examples/benchmark.jl) reproduces these timings.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="benchmark/tolerances-dark.png">
  <img alt="Two panels against relative tolerance from 1e-3 to 1e-9: error, and run time, for Julia QNDF, FBDF, Rodas5P and KenCarp47 and for MATLAB ode15s with and without JPattern." src="benchmark/tolerances.png">
</picture>

Each solver integrated the protocol at seven tolerances, with `abstol = reltol / 100`.

- **Error.** The largest deviation of any state from a reference solution, relative to
  that state's largest value. It is checked at the end of every meal interval and every
  10 minutes on the last day.
- **Reference.** Rodas5P at `reltol = 1e-12`. It is accurate to about 1e-10 and agrees
  with ode15s at `1e-10` to within 1.5e-9.
- **MATLAB.** One `ode15s` call per meal interval, as in `simulate_human_day`. It runs
  once as shipped, and once with the Jacobian sparsity pattern passed as `JPattern`.
- **Julia timings.** The fastest of three runs for runs under 5 s.

What the sweep shows:

- **The same tolerance does not give the same accuracy.** At equal tolerance, ode15s
  is 3–6× more accurate than QNDF. At `reltol = 1e-6`, for example, the error is 2.6e-6
  for ode15s and 1.6e-5 for QNDF (its default). Compared at equal error, an error of
  about 2e-6 takes:
  - 1.0 s with QNDF (`reltol = 1e-7`);
  - 13 s with ode15s + JPattern (`reltol = 1e-6`);
  - 45 s with ode15s as shipped (`reltol = 1e-6`).
- **Much of MATLAB's cost is the dense finite-difference Jacobian.** `JPattern` makes
  ode15s 2.2–3.4× faster. QNDF is still about 13× faster than that at equal error.
- **FBDF's error stalls at about 4e-3** until `reltol = 1e-7`; avoid it for this model.
  KenCarp47 converges slowly. Rodas5P is accurate, but 3–10× slower than QNDF.

At `reltol = 1e-3` the sweep's ode15s run takes 19 s rather than 43 s. The 43 s run of
`simulate_human_day` used ode15s's default `AbsTol = 1e-6`, and also includes evaluating
the diagnostic variables at every output time of the last day.

<details>
<summary>Results as a table (run time, error)</summary>

| reltol | Julia QNDF (default) | MATLAB ode15s | MATLAB ode15s + JPattern | Julia Rodas5P | Julia KenCarp47 | Julia FBDF |
|---|---|---|---|---|---|---|
| 1e-03 | 0.45 s, 4.5e-03 | 19 s, 1.1e-03 | 7 s, 1.1e-03 | 1.2 s, 3.0e-04 | 0.63 s, 3.1e-03 | 0.42 s, 4.4e-03 |
| 1e-04 | 0.4 s, 8.9e-04 | 25 s, 1.5e-04 | 9 s, 1.5e-04 | 1.8 s, 3.3e-04 | 0.89 s, 6.0e-04 | 0.54 s, 4.0e-03 |
| 1e-05 | 0.53 s, 8.7e-05 | 34 s, 2.2e-05 | 11 s, 2.2e-05 | 2.7 s, 6.8e-05 | 1.3 s, 2.3e-04 | 0.65 s, 4.0e-03 |
| 1e-06 | 0.73 s, 1.6e-05 | 45 s, 2.6e-06 | 13 s, 2.6e-06 | 4.2 s, 1.1e-05 | 2 s, 4.0e-05 | 0.91 s, 4.4e-03 |
| 1e-07 | 1 s, 1.7e-06 | 54 s, 4.0e-07 | 17 s, 4.0e-07 | 6.5 s, 1.0e-06 | 2.9 s, 2.0e-05 | 1.2 s, 1.6e-03 |
| 1e-08 | 1.3 s, 2.3e-07 | 48 s, 6.2e-08 | 22 s, 6.2e-08 | 10 s, 1.2e-07 | 4.6 s, 5.9e-06 | 1.7 s, 1.3e-05 |
| 1e-09 | 1.6 s, 3.1e-08 | 65 s, 9.5e-09 | 29 s, 9.5e-09 | 16 s, 9.9e-09 | 7 s, 1.7e-06 | 2.6 s, 2.0e-06 |

</details>

To reproduce, run from the repository root. The scripts and results are in
[`benchmark/`](benchmark), and the commands take about 2, 7 and 3 minutes:

```
julia --project=benchmark -e "using Pkg; Pkg.instantiate()"
julia --project=benchmark benchmark/reference.jl     # reference solution
matlab -batch "run('benchmark/sweep_matlab.m')"      # MATLAB sweep
julia --project=benchmark benchmark/sweep_julia.jl   # Julia sweep
julia --project=benchmark benchmark/plot.jl          # figures and table
```

## Differences from the MATLAB implementation

The equations are the same. The implementation differs in a few places:

- One solve over the whole protocol, with restarts at meals, instead of one solve per
  meal. As in MATLAB, only the most recent meal drives the postprandial reflex, and the
  fasting period starts with a postprandial response at `t = 0` (MATLAB's `tau = 0`).
- Hepatic extraction is written as extraction from portal plus systemic blood. This is
  algebraically identical to MATLAB's separate first-pass and plasma-clearance terms.
- As in MATLAB, unconjugated bile acids in liver and gallbladder stay at their initial
  value (zero), because they are conjugated on hepatic uptake.
- `Constants` only holds the 15 constants that enter the equations. Most of the 128
  entries in MATLAB's constants vector are unused.
- Not ported: the transit tracer model and `simulate_transit_exp`, the diagnostic
  variables of `*_var_full.m`, the cost function and calibration code, and plotting.

## Tests

```
julia --project=. -e "using Pkg; Pkg.test()"
```

The tests compare against data from the MATLAB implementation in
[`test/reference`](test/reference):

- the state layout;
- the right-hand side at 40 random states. Half use the calibrated parameters; half use
  perturbed parameters with every constant that is zero in the calibrated set switched on;
- full simulations against a tight-tolerance MATLAB solution.

They also check mass conservation, the Jacobian sparsity pattern, parameter derivatives
and zero allocations. To regenerate the reference data, run
`test/reference/generate_reference.m` in MATLAB (takes about 3 minutes).
