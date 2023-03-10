"""
```jldoctest
    shannon(p::AbstractArray)
```
Simple function to implement Shannon entropy given an array of probabilities `p`.
Requires `sum(p) == 1`.
"""
function shannon(p::AbstractArray)

    # check requirement
    if sum(p) != 0
        error("Probability distribution needs to add up to 1")
    end

    # check for dimension mismatch
	if size(A) != (N,N)
		error("Dimension mismatch with adj matrix")
	end

    # compute H(p)
	out = 0
	for pi in p
		if pi > 0
			out += pi * log2(1 / pi)
		end
	end

	return out
end

"""
```jldoctest
    effective_info(A::AbstractMatrix)
```
Simple function to calculate the effective entropy of a weighted graph.
"""
function effective_info(A::AbstractMatrix)
    
    N = size(A)[1]  # size of the system
    
    # check for dimension mismatch
	if size(A) != (N,N)
		error("Dimension mismatch with adj matrix")
	end

    det = log2(4) - mean(map(shannon, eachrow(A)))  # determinism
    deg = log2(4) - shannon(sum(A, dims = 1)./ N)  # degeneracy
    EI = det - deg

    return EI, deg, det
end