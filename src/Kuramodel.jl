module Kuramodel

# export functions
export greet_Kuramodel
export greet
export randθ # this is to test how to add dependencies to a package
export Kura_step
export Kurasim

# where are these functions from?
include("ExampleFunction.jl")
include("Kuramoto.jl")
end
