# Example Function
function greet_Kuramodel()
    println("Hello Kuramodel user!")
end

function greet()
    println("You are using the Kuramodel Package, have fun!")
end

using Distributions

# Function to generate array of random phases between -π and π
function randθ(N::Number)
    return rand(Uniform(-π, π),N)
end