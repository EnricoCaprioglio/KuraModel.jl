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

include("Utilities.jl")
export plot_local_op
export ω_locals
export θs_locals
export θs_macro
export ω_macro
export get_splits
export fun_pattern

include("InformationTheory.jl")
export shannon
export effective_info

include("ExampleFunction.jl")

end