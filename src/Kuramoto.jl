using Distributions
using Random

"""
Kura_step(σ::AbstractArray,ω::AbstractArray,A::Matrix,dt::Number;
        θ=nothing::AbstractArray,noise_scale=nothing,seedval=nothing,τ=nothing)

Simple Euler's method to compute the step evolution of N Kuramoto oscillators.
    
    Inputs:
        σ::AbstractArray    array of couplings
        ω::AbstractArray    array of natural frequencies
        A::Matrix           adjacency matrix of the underlying network
        Δt::Number          time-step for Euler's method
    Optional Inputs:
        θ::AbstractArray    instantaneous phases
        noise_scale::Number noise coefficient
        τ::Number           phase delay
    
    Output:
        θ::AbstractArray    updated phases (1 Δt step)

## Example:
```jldoctest
julia> Kura_step(
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
function Kura_step(σ::AbstractArray,ω::AbstractArray,A::Matrix,Δt::Number;
	θ=nothing::AbstractArray,τ=nothing,noise_scale=nothing,seedval=nothing)
	# number of nodes N
	N=length(ω)
	
	# set optional arguments
	if θ==nothing
		θ=rand(Uniform(-π, π),N)
	end
	if noise_scale == nothing
		noise_scale=0
	end
	if seedval≢nothing
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
	if τ == nothing
		θj_θi_mat=repeat(θ',N)-repeat(θ',N)' # matrix of θ_j-θ_i
	else
		θj_θi_mat=(repeat(θ',N)-repeat(θ',N)').-τ
		setindex!.(Ref(θj_θi_mat), 0.0, 1:N, 1:N) # set diagonal elements to zero (self-interaction terms)
	end
	k1 = *(A.*sin.(θj_θi_mat), σ)
	θ += Δt.*(ω + k1) + noise_scale*(rand(Normal(0,1),N))*sqrt(Δt)		
	return θ
end

# Probably no need to include the seedvalue in the step function as well, since it is
# already in the Kurasim function

"""
    Kurasim(σ::AbstractArray,ω::AbstractArray,A::Matrix,t_tot::Integer,Δt::Number;
        θ0=nothing::AbstractArray,τ=nothing,noise_scale=nothing,seedval=nothing)

Function to execute a Kuramoto simulation.

    Inputs:
        σ::AbstractArray    array of couplings
        ω::AbstractArray    array of natural frequencies
        A::Matrix           adjacency matrix of the underlying network
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
t=10; Δt=0.02; # time-steps
steps=collect(0:Δt:t-Δt) # total number of steps calculated
ω=rand(-2:0.0000001:2,N); # intrinsic frequencies
σ=[1,1]; # couplings
A=[0 1;1 0]; # initialize network
θs=Kuramodel.Kurasim(σ,ω,A,t,Δt,θ0=nothing,noise_scale=0.0,τ=nothing,seedval=nothing)

50×2 Matrix{Float64}:
 -1.89658  -0.200719
 -1.97067  -0.400048
 -2.0432   -0.600938
 -2.11738  -0.800179
  ⋮        
 -7.71061  -6.96407
 -7.84732  -7.10078
 -7.98403  -7.23749
 -8.12074  -7.3742
 ```
"""
function Kurasim(σ::AbstractArray,ω::AbstractArray,A::Matrix,t_tot::Integer,Δt::Number;
	θ0=nothing::AbstractArray,noise_scale=nothing,seedval=nothing,τ=nothing)
	# number of nodes N
	N=length(ω)
	# set optional arguments
	if θ0==nothing
		θ0=rand(Uniform(-π, π),N)
	end
	if noise_scale == nothing
		noise_scale=0
	end
	if seedval≢nothing
		# println("Seed value: $(round(Int,seedval))")
		Random.seed!(round(Int,seedval));
	end
	# check for dimension mismatch
	if size(A) ≢ (N,N)
		error("Dimension mismatch with adj matrix")
	end
	if length(θ0) ≢ N
		error("Dimension mismatch with phases array")
	end
	if length(σ) ≢ N
		error("Dimension mismatch with couplings array")
	end
	# Initialize matrix to store phases: (t x N) matrix
	# each row stores the instantaneous phases of all N oscillators
	tot_steps=round(Integer,t_tot/Δt)
	θs=zeros(tot_steps,N);
	# record initial phases at t=1
	θs[1,:]=θ0
	for t in 2:tot_steps
		# update phases
		θs[t,:]=Kura_step(σ,ω,A,Δt;θ=θs[t-1,:],noise_scale,seedval,τ)
	end
    println("Simulation completed, total steps: $(tot_steps)")
	return θs
end