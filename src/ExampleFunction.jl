# Example Functions (these were created to test the package was correclty loaded)
function greet_Kuramodel()
    println("Hello Kuramodel user!")
end

function greet()
    println("You are using the Kuramodel Package, have fun!")
end

using Distributions
using Random

"""
Simple function to generate array of random phases between -π and π.
"""
function randθ(N::Number)
    return rand(Uniform(-π, π),N)
end