using Distributions
using Kuramodel

"""
    function fun_pattern(A,θs,t_span)

Function used to calculate the functional pattern (correlation matrix)
of a system of oscillators.

    Inputs:
        A::AbstractMatrix   adjacency matrix of the underlying network
        θs::AbstractArray   instantaneous phases at each timestep
        t_span              time window of interest (can be a tuple or vector of order 2)

    Output:
        R::AbstractMatrix   functional pattern matrix for time window t_span
"""
function fun_pattern(A::AbstractMatrix,θs::AbstractArray,t_span)
    R=zeros(size(A))
    for i in 1:length(A[1,:])
        for j in 1:length(A[1,:])
            # local corr_array=[]
            local corr_sum=0
            for t in t_span[1]:t_span[2]
                # append!(corr_array,cos(θs[t,j]-θs[t,i]))
                corr_sum+=cos(θs[t,j]-θs[t,i])
            end
            # R[i,j]=mean(corr_array)
            R[i,j]=corr_sum/(t_span[2]-t_span[1]+1)
        end
    end
    return R
end


# probably I can make this faster by subratcting the transpose of the matrix at each timestep
# not sure if I will need a fater function so for now I can keep it like this for readability

"""
    function ω_macro(ω::AbstractArray,Nc::AbstractArray)
This function finds the local mean natural frequency for each community.

    Inputs:
        ω::AbstractArray		the vector of frequencies.
        C::Integer 				the number of communities.
        M::Integer 				the number of oscillators per community.

    Output:
        means::AbstractVector   the natural frequencies of each macroscillator
"""
function ω_macro(ω::AbstractArray,Nc::AbstractArray)
    splits = get_splits(Nc)
	means=[]
	for i in 1:length(splits)
        if i%2 == 1
            push!(means,mean(ω[splits[i]:splits[i+1]]))
        end
    end
	return means
end

"""
	function θs_locals(θs::AbstractArray,C::Integer,M::Integer)
This function finds the local mean phase for each community at each time step.

    Inputs:
        θs::AbstractArray		    the matrix of instantaneous phases.
        C::Integer 				    the number of communities.
        M::Integer 				    the number of oscillators per community.
    Output:
        θs_means::AbstractArray     instantaneous community averaged phases at each timestep

The output is a matrix (timesteps × M).
"""
function θs_macro(θs::AbstractArray,Nc)
	tot_steps=length(θs[:,1])
	θs_means=zeros(tot_steps,length(Nc))
    splits = get_splits(Nc)
	for t in 1:tot_steps
		# for i in 1:M:C*M
		# 	c_temp=round(Integer,ceil(i/M))
		# 	θs_means[t,c_temp]=mean(θs[t,i:i+M-1])
		# end
        for i in 1:length(splits)
            if i%2 == 1
                c=Integer(i/2) # find community index
                θs_means[t,c_temp]=mean(θs[t,splits[i]:splits[i+1]])
            end
        end
	end
	return θs_means
end