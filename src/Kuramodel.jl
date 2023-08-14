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

include("InformationTheory.jl")
export shannon
export effective_info

include("ExampleFunction.jl")

end