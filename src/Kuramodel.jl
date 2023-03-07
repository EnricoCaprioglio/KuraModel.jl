module Kuramodel

# export functions
export greet
export randθ
export Kurastep
export Kurasim
export macro_op
export get_splits
export fun_pattern
export θs_macro
export ω_macro

export Kura_obj
export Kura_step
export Kura_sim
export getsize

export plot_local_op
export ω_locals
export θs_locals

# where are these functions from?
include("ExampleFunction.jl")
include("Kuramoto.jl")
include("Utilities.jl")
include("CommunityStructured.jl")
end