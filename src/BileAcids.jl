"""
    BileAcids

Compartmental model of human bile acid metabolism by Sips et al., "In silico analysis
identifies intestinal transit as a key determinant of systemic bile acid metabolism",
ported from the MATLAB implementation in `BA_model_package`.

Seven bile acid species, each unconjugated or conjugated with taurine or glycine, are
followed through liver, gallbladder, ten small-intestinal and five colon segments,
faeces and plasma.
"""
module BileAcids

using OrdinaryDiffEq
using OrdinaryDiffEqBDF: QNDF
using SciMLBase: SciMLBase
using SparseArrays: SparseMatrixCSC
using SparseConnectivityTracer: SparseConnectivityTracer, TracerSparsityDetector
using StaticArrays
using PrecompileTools: @setup_workload, @compile_workload

include("layout.jl")
include("parameters.jl")
include("model.jl")
include("meals.jl")
include("compartments.jl")
include("rhs.jl")
include("problem.jl")
include("precompile.jl")

export Parameters, Constants, BileAcidModel, MealSchedule, GutSegment
export rhs!, simulate, with_meals, last_day, gi_reflex, jacobian_prototype, meal_callback
export compartment, state_index, state_names
export SPECIES, FORMS, COMPARTMENTS, N_STATES

end
