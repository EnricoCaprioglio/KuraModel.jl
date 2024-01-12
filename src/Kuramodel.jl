module Kuramodel

include("Kuramoto.jl")
export greet
export randθ
export Kurastep
export Kurasim
export macro_op
export Kura_obj
export Kura_step
export Kura_sim
export getsize
export Kura_step_D
export Kura_sim_D

include("Utilities.jl")
export plot_local_op
export ω_locals
export θs_locals
export θs_macro
export ω_macro
export get_splits
export fun_pattern
export macro_to_micro
export get_neighbours
export Laplacian

include("InformationTheory.jl")
export shannon
export effective_info

include("ExampleFunction.jl")

include("graphs/generators/Functions.jl")
export get_alpha
export get_ξ
export px_HRGG
export px_softHRGG
export integrand_HRGG
export integrand_softHRGG
export get_avg_k_HRGG
export get_avg_k_softHRGG
export get_R_HRGG
export get_R_softHRGG
export sample_r
export sample_ϕ
export hyper_d
export connect_nodes_HRGG
export connect_nodes_softHRGG

include("graphs/generators/GeometricGraphs.jl")
export hyperbolic_graph

include("HierarchicalChimera.jl")
export SBMvar
export _SBMvar_constructor

end