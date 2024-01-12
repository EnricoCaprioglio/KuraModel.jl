using Kuramodel

"""
	function SBMvar(n::AbstractArray, B::AbstractArray, H::Number, k::Number; K = 1)

Function to crate the nested SBM variation with 3 hierarchical levels.

Inputs:

	n::AbstractArray 		vector of length 3, number of l-1 blocks in each l block
	B::AbstractArray 		number of blocks in layer l
	H::Number 				value of hierarchical parameter
	k::Number 				desired average degree, k ∈ [n₁-1, n₁n₂-1]
	K::Number 				set the coupling (weighted matrix output)

Outputs:

	A::AbstractMatrix 		 		adjacency matrix
	P::AbstractMatrix 				partition matrix
	[p₁, p₂, p₃]::AbstractArray 	probabilities of connection at layer l

"""
function SBMvar(n::AbstractArray, B::AbstractArray, H::Number, k::Number; K = 1)

	γ = (k + 1 - n[1]) / (n[1]*n[2] - n[1])

	return _SBMvar_constructor(n, B, H, γ; K = K)
	
end

"""
	function _SBMvar_constructor(n::AbstractArray, B::AbstractArray, H::Number, γ::Number; K = 1)

Function to crate the nested SBM variation with 3 hierarchical levels.

Inputs:

	n::AbstractArray 		vector of length 3, number of l-1 blocks in each l block
	B::AbstractArray 		number of blocks in layer l
	H::Number 				value of hierarchical parameter
	γ::Number 				parameter to control the average degree, γ ∈ [0,1]
	K::Number 				set the coupling (weighted matrix output)

Outputs:

	A::AbstractMatrix 		 		adjacency matrix
	P::AbstractMatrix 				partition matrix
	[p₁, p₂, p₃]::AbstractArray 	probabilities of connection at layer l

"""
function _SBMvar_constructor(n::AbstractArray, B::AbstractArray, H::Number, γ::Number; K = 1)

	if length(n) != length(B)
		error("vectors n and B mismatch")
	end
	if γ < 0 || γ > 1
		error("require 0 ≤ γ ≤ 1")
	end
	if B[1] != B[2]*n[2]
		error("B₁ != B₂ × n₂, check model conditions")
	end
	if B[2] != B[3]*n[3]
		error("B₂ != B₃ × n₃, check model conditions")
	end

	N = prod(n) # number of nodes
	L = length(n) # number of layers
	
	# store partition and adjaecncy matrix
	P = zeros(L - 1, N)
	A = zeros(N, N)

	# compute edge probabilities
	p₁ = (n[1]*H*γ) / (2*(n[1]-1)) + (1 - (n[1]*1*γ) / (2*(n[1]-1)))
	p₂ = (1/2 + H/2) * γ
	p₃ = (1/2 - H/2) * γ

	# create partition
	for i in 1:(L-1)
		P[i, :] = vcat(
			[repeat([j], prod(n[1:i])) for j in 1:B[i]]...
		)
	end

	for i in 1:N
		for j in i+1:N

			# connect if in the same l = 1 block only
			if P[:, i] == P[:, j] && rand() < p₁
				A[i, j] = K
				A[j, i] = K
				
			# connect if in the same l = 2 block only
			elseif P[2, i] == P[2, j] && rand() < p₂
				A[i, j] = K
				A[j, i] = K

			# connect if in the same l = 3 block only
			elseif P[2, i] != P[2, j] && rand() < p₃
				A[i, j] = K
				A[j, i] = K
				
			end
		end
	end

	return A, P, [p₁, p₂, p₃]
end