using Kuramodel
using Plots

# to call some function from the terminal, in order to sepcify the argument
# which could be some environment variable from the job description
# we can use the following:

# julia -L file.jl -e 'main(some,args)'

function main(arg)
    println("This is the input to the main function: $arg")
end

# println("Before greet()")

# greet()

# println("After greet()")

# some calculations:
a = zeros(2,2)
b = randn(5)
c = "some string"

println("Very nice calculations!\n")
println("ARGS is the global variable for the arguments: $(ARGS)\n")
println("While PROGRAM_FILE is the global variable of the name of the script: $(PROGRAM_FILE)")

plot(b)
