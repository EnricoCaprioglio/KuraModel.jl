using Distributions
using Random
"""
Quick function to test package has been loaded correctly.
```jldoctest
julia> greet()
You are using the Kuramodel Package, have fun!
```
"""
function greet()
    println("You are using the Kuramodel Package, have fun!")
end

"""
Simple function to generate an array of random phases between -π and π.
"""
function randθ(N::Number)
    return rand(Uniform(-π, π),N)
end