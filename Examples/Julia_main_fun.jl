using Kuramodel
using JLD2
using Random

# to call some function from the terminal, in order to sepcify the argument
# which could be some environment variable from the job description
# we can use the following:

# function main(arg)
#     println("This is the input to the main function: $arg")
# end

# julia -L file.jl -e 'main(some,args)'
# this imitates python's main function and call from the terminal.
# However, in Julia we can just call all the functions in a file like this:
# julia file.jl arg1 arg2 ...

println("- ARGS is the global variable for the arguments: $(ARGS)\n")

println("- While PROGRAM_FILE is the global variable of the name of 
    the script: $(PROGRAM_FILE)\n"
    )

println("This is the current working directory: ", pwd())

# set path
store_data_path = "/mnt/nfs2/inf/ec627/src/Kuramodel/Examples/DataTest/"
filenameroot = "randomData"
arrayname = "$(ARGS)"
extension = ".jld2"

# create some random reproducible data
Random.seed!(123)
a = $(ARGS)

# save test object
save_object(store_data_path * filenameroot * arrayname * extension, a)