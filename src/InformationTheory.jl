"""
Simple function to implement Shannon entropy given an array of probabilities `p`.
"""
function shannon(p::AbstractArray)
	out = 0
	for pi in p
		if pi > 0
			out += pi*log2(1/pi)
		end
	end
	return out
end

"""
Simple function to calculate the effective entropy of a weighted graph.
"""
function effective_info(A::AbstractMatrix)
    
    N = size(A)[1]  # size of the system
    det = log2(4) - mean(map(shannon, eachrow(A)))  # determinism
    deg = log2(4) - shannon(sum(A, dims = 1)./ N)  # degeneracy
    EI = det - deg

    return EI, deg, det
end