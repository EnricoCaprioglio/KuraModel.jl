"""
    function get_splits(Nc)
This function is used in many other functions, such as ω_macro or θs_macro.
Given an Array containing the number of nodes in each community, such that
`Nc[i]` is the number of nodes in community `c_{i}` it outputs a useful iterator
to perform calculations for each community.

See functions `ω_macro` or `θs_macro` to see how it us used, or the example below
to find the average ω for each community.
## Example: 
```jldoctest
julia> Nc = [4,3,2,3]
splits = get_splits(Nc)
>> 8-element Vector{Any}:
  1
  4
  5
  7
  8
  9
 10
 12
 ω = [1,2,3,2,4,5,6,3,3,7,8,9]
 means=[]
 for i in 1:length(splits)
     if i%2 == 1
         push!(means,mean(ω[splits[i]:splits[i+1]]))
     end
 end
 means
>> 4-element Vector{Any}:
2.0
5.0
3.0
8.0
 ```
"""
function get_splits(Nc)
    splits = []
	for (i,c) in enumerate(Nc)
		append!(splits,sum(Nc[1:i-1])+1, sum(Nc[1:i-1])+Nc[i])
	end
    return splits
end


"""
```jldoctest
	function plot_local_op(θs::AbstractArray,Nc::AbstractArray)
```
Function used to plot the local order parameters for each community `c`.
To plot use `plot(plot_local_op(θs, Nc)..., kwargs)`.
"""
function plot_local_op(θs::AbstractArray, Nc::AbstractArray, steps_to_plot=nothing)
	ps = []
	splits = get_splits(Nc)

	for i in 1:length(splits)
		if i%2 == 1
			c = round(Integer,ceil(i/2))
			if steps_to_plot === nothing
				push!(ps,
					plot(
						macro_op(θs[:,splits[i]:splits[i+1]]),
						label = "",
						ylims = [0.,1.],
						bg = RGB(0.2, 0.2, 0.2),
						title = L"R_{%$(c)}(t)"
					)
				)
			else
				push!(ps,
					plot(
						steps_to_plot,
						macro_op(θs[:,splits[i]:splits[i+1]]),
						label = "",
						ylims = [0.,1.],
						bg = RGB(0.2, 0.2, 0.2),
						title = L"R_{%$(c)}(t)"
					)
				)
			end
		end
	end
	return ps
end

"""
```jldoctest
	function ω_locals(ω::AbstractArray, Nc::AbstractArray)
```
this function finds the local mean natural frequency for each community.
"""
function ω_locals(ω::AbstractArray, Nc::AbstractArray)

	splits = get_splits(Nc)
	means=[]

	for i in 1:length(splits)
		if i%2 == 1
			push!(means, mean(ω[splits[i]:splits[i+1]]))
		end
	end
	return means
end

"""
```jldoctest
	function θs_locals(θs::AbstractArray, Nc::AbstractArray)
```
this function finds the local mean phase for each community at each time step.
The output is a matrix (timesteps × M)
"""
function θs_locals(θs::AbstractArray, Nc::AbstractArray)
	
	splits = get_splits(Nc)
	tot_steps = length(θs[:,1])
	θs_means = zeros(tot_steps, length(Nc))
	
	for t in 1:tot_steps
		for i in 1:length(splits)
			if i%2 == 1
				c_temp=round(Integer, ceil(i/2))
				θs_means[t,c_temp] = mean(θs[t, splits[i]:splits[i+1]])
			end
		end
	end
	return θs_means
end