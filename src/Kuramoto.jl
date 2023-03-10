using Distributions
using Random

## TODO
#1 Probably no need to include the seedvalue in the step function as well, since it is
# already in the Kurasim function
#2 check if there is a better way to add options to functions, this is okay for now
#3 

"""
```jldoctest
	Kurastep(σ::AbstractArray,ω::AbstractArray,A::AbstractMatrix,dt::Number;
        θ=nothing::AbstractArray,noise_scale=nothing,seedval=nothing,τ=nothing)
```

Simple Euler's method to compute the step evolution of N Kuramoto oscillators.
    
    Inputs:
        σ::AbstractArray    array of couplings
        ω::AbstractArray    array of natural frequencies
        A::AbstractMatrix   adjacency matrix of the underlying network
        Δt::Number          time-step for Euler's method
    Optional Inputs:
        θ::AbstractArray    instantaneous phases
        noise_scale::Number noise coefficient
        τ::Number           phase delay
    
    Output:
        θ::AbstractArray    updated phases (1 Δt step)

## Example:
```jldoctest
julia> Kurastep(
        [1,1],
        [1,1],
        [0 1; 1 0],
        0.02,
        θ=[2,0.5]
    )
2-element Vector{Float64}:
2.0000501002679187
0.5399498997320811
```
"""
function Kurastep(σ::AbstractArray,ω::AbstractArray,A::AbstractMatrix,Δt::Number;
	θ=nothing,τ=nothing,noise_scale=nothing,seedval=nothing)
	# number of nodes N
	N=length(ω)
	
	# set optional arguments
	if θ === nothing
		θ=rand(Uniform(-π, π),N)
	end
	if noise_scale === nothing
		noise_scale=0
	end
	if seedval !== nothing
		# println("Seed value: $(round(Int,seedval))")
		Random.seed!(round(Int,seedval));
	end
	
	# check for dimension mismatch
	if size(A) ≢ (N,N)
		error("Dimension mismatch with adj matrix")
	end
	if length(θ) ≢ N
		error("Dimension mismatch with phases array")
	end
	if length(σ) ≢ N
		error("Dimension mismatch with couplings array")
	end

	# update dθ
	if τ === nothing
		θj_θi_mat=repeat(θ',N)-repeat(θ',N)' # matrix of θ_j-θ_i
	else
		θj_θi_mat=(repeat(θ',N)-repeat(θ',N)').-τ
		setindex!.(Ref(θj_θi_mat), 0.0, 1:N, 1:N) # set diagonal elements to zero (self-interaction terms)
	end
	k1 = map(sum,eachrow(A.*sin.(θj_θi_mat))).*σ
	θ += Δt.*(ω + k1) + noise_scale*(rand(Normal(0,1),N))*sqrt(Δt)		
	return θ
end

"""
	function Kurasim(σ::AbstractArray,ω::AbstractArray,A::AbstractMatrix,t_tot::Integer,Δt::Number;
		θ0=nothing,noise_scale=nothing,seedval=nothing,τ=nothing)

Function to execute a Kuramoto simulation.

    Inputs:
        σ::AbstractArray    array of couplings
        ω::AbstractArray    array of natural frequencies
        A::AbstractMatrix   adjacency matrix of the underlying network
        t_tot::Integer      total time for the simulation
        Δt::Number          time-step for Euler's method
    Optional Inputs:
        θ::AbstractArray    initial instantaneous phases
        noise_scale::Number noise coefficient
        τ::Number           phase delay

    Output:
        θs::AbstractArray   updated phases at each timestep (tot number of steps=t_tot/Δt)

## Example of Kuramoto simulation
```jldoctest
julia> seedvalue=4; Random.seed!(seedvalue); N=2;
# time-steps
t=10; Δt=0.2;
steps=collect(0:Δt:t-Δt) # total number of steps calculated
# initialize natural frequencies
ω=rand(-2:0.0000001:2,N);
# couplings
σ=[1,1];
# initialize network
A=[0 1;1 0];
# start simulation
θs=Kuramodel.Kurasim(σ,ω,A,t,Δt)
50×2 Matrix{Float64}:
  1.11332   0.0218031
  1.21944   0.545724
  1.37825   1.01696
  1.59114   1.43411
  ⋮        
 14.9796   15.1375
 15.2946   15.4525
 15.6097   15.7675
 15.9247   16.0825
```
"""
function Kurasim(σ::AbstractArray,ω::AbstractArray,A::AbstractMatrix,t_tot::Integer,Δt::Number;
	θ0=nothing,noise_scale=nothing,seedval=nothing,τ=nothing)

	# number of nodes N
	N=length(ω)

	# set optional arguments
	if θ0 === nothing
		θ0 = rand(Uniform(-π, π), N)
	end

	if noise_scale === nothing
		noise_scale = 0
	end

	if seedval !== nothing
		# println("Seed value: $(round(Int,seedval))")
		Random.seed!(round(Int,seedval));
	end

	# check for dimension mismatch
	if size(A) != (N,N)
		error("Dimension mismatch with adj matrix")
	end

	if length(θ0) != N
		error("Dimension mismatch with phases array")
	end

	if length(σ) != N
		error("Dimension mismatch with couplings array")
	end

	# Initialize matrix to store phases: (t x N) matrix
	# each row stores the instantaneous phases of all N oscillators
	tot_steps = round(Integer, t_tot / Δt)
	θs = zeros(tot_steps, N)
	# record initial phases at t=1
	θs[1,:] = θ0

	for t in 2:tot_steps
		# update phases
		θs[t,:] = Kurastep(σ, ω, A, Δt; θ = θs[t-1,:], τ = τ, noise_scale = noise_scale, seedval = seedval)
	end

    println("Simulation completed, total steps: $(tot_steps)")

	return θs
end

"""
```jldoctest
	macro_op(θs::AbstractMatrix)
```
Function to calculate the instantaneous macroscopic parameter for each row
of the instantaneous phases stored in θs.

	Input:
		θs::AbstractMatrix	(tot_num_timesteps × N) matrix, each row
							contains the instantaneous phases of all
							N oscillators.
	Output:
		R::AbstractArray	Instantaneous oreder parameter for each 
							timestep.

"""
function macro_op(θs::AbstractMatrix)
	# num of cols = N
	N = length(θs[1, : ])
    map(θ -> abs(sum(exp.(im * θ) / N)), eachrow(θs))
end

mutable struct Kura_obj
	σ::AbstractArray
	ω::AbstractArray
	A::AbstractMatrix
	θ::AbstractArray
	# instantiate
	Kura_obj(σ, ω, A, θ) = new(σ, ω, A, θ)
	# random initial phases
	Kura_obj(σ, ω, A) = new(σ, ω, A, rand(Uniform(-π, π), N))
	# random initial phases and complete graph
	Kura_obj(σ, ω) = new(σ, ω, ones(length(σ),length(σ)), rand(Uniform(-π, π),N))
end

getsize(x::Kura_obj) = size(x.A)

function Kura_step(kobj::Kura_obj, Δt::Number; τ = 0.0, noise_scale = 0.0)
	# get size
	N = getsize(kobj)[1];
	# update dθ
	if τ === 0.0
		θj_θi_mat = repeat(kobj.θ',N)-repeat(kobj.θ',N)' # matrix of θ_j-θ_i
	else
		θj_θi_mat = (repeat(kobj.θ',N)-repeat(kobj.θ',N)').-τ
		# setindex!.(Ref(θj_θi_mat), 0.0, 1:N, 1:N) # set diagonal elements to zero (self-interaction terms)
	end
	k1 = *(kobj.A.*sin.(θj_θi_mat), kobj.σ)
	kobj.θ += Δt.*(kobj.ω + k1) + noise_scale*(rand(Normal(0,1),N))*sqrt(Δt)
end

function Kura_sim(kobj::Kura_obj, t_tot, Δt::Number; message = true)
	N = getsize(kobj)[1];
	θs = zeros(t_tot,N);  # store thetas
	θs[1,:] = kobj.θ;  # initial thetas
	for i in 2:t_tot
		Kura_step(kobj,Δt)
		θs[i,:] = kobj.θ;  # update current thetas
	end
	if message === true
		println("Simulation completed, total steps: $(t_tot)")
	end
	return θs
end

# """
# 	Kura_step(σ::AbstractArray,ω::AbstractArray,A::AbstractMatrix,dt::Number;
#         θ=nothing::AbstractArray,noise_scale=nothing,seedval=nothing,τ=nothing)

# Simple Euler's method to compute the step evolution of N Kuramoto oscillators.
    
#     Inputs:
#         σ::AbstractArray    array of couplings
#         ω::AbstractArray    array of natural frequencies
#         A::AbstractMatrix   adjacency matrix of the underlying network
#         Δt::Number          time-step for Euler's method
#     Optional Inputs:
#         θ::AbstractArray    instantaneous phases
#         noise_scale::Number noise coefficient
#         τ::Number           phase delay
    
#     Output:
#         θ::AbstractArray    updated phases (1 Δt step)

# ## Example:
# ```jldoctest
# julia> Kura_step(
#         [1,1],
#         [1,1],
#         [0 1; 1 0],
#         0.02,
#         θ=[2,0.5]
#     )
# 2-element Vector{Float64}:
# 2.0000501002679187
# 0.5399498997320811
# ```

# """
# function Kura_step(σ::AbstractArray,ω::AbstractArray,A::AbstractMatrix,Δt::Number;
# 	θ=nothing,τ=nothing,noise_scale=nothing,seedval=nothing)
# 	# number of nodes N
# 	N=length(ω)
	
# 	# set optional arguments
# 	if θ === nothing
# 		θ = rand(Uniform(-π, π),N)
# 	end
# 	if noise_scale === nothing
# 		noise_scale = 0
# 	end
# 	if seedval ≢ nothing
# 		# println("Seed value: $(round(Int,seedval))")
# 		Random.seed!(round(Int,seedval));
# 	end
	
# 	# check for dimension mismatch
# 	if size(A) ≢ (N,N)
# 		error("Dimension mismatch with adj matrix")
# 	end
# 	if length(θ) ≢ N
# 		error("Dimension mismatch with phases array")
# 	end
# 	if length(σ) ≢ N
# 		error("Dimension mismatch with couplings array")
# 	end

# 	# update dθ
# 	if τ === nothing
# 		θj_θi_mat = repeat(θ',N)-repeat(θ',N)' # matrix of θ_j-θ_i
# 	else
# 		θj_θi_mat = (repeat(θ',N)-repeat(θ',N)').-τ
# 		setindex!.(Ref(θj_θi_mat), 0.0, 1:N, 1:N) # set diagonal elements to zero (self-interaction terms)
# 	end
# 	k1 = *(A.*sin.(θj_θi_mat), σ)
# 	θ += Δt.*(ω + k1) + noise_scale*(rand(Normal(0,1),N))*sqrt(Δt)		
# 	return θ
# end