module Kuramodel

# export functions
export greet
export randθ
export Kura_step
export Kurastep
export Kurasim
export macro_op
export get_splits
export fun_pattern
export θs_macro
export ω_macro

# where are these functions from?
include("ExampleFunction.jl")
include("Kuramoto.jl")
include("Utilities.jl")
include("CommunityStructured.jl")
end