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