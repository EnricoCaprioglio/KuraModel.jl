using Kuramodel

# to call some function from the terminal, in order to sepcify the argument
# which could be some environment variable from the job description
# we can use the following:

# julia -L file.jl -e 'main(some,args)'

function main(arg)
    println("This is the input to the main function: $arg")
end

println("Before greet()")

greet()

println("After greet()")