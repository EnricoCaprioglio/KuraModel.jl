# Example Kuramoto
using Kuramodel
using Random
greet()

# set seedvalue
seedvalue=4
Random.seed!(seedvalue)

# Example of Kura_step
noise_scale=0.2; τ=0.3;
Kura_step([1,1],[1,1],[0 1; 1 0],0.02;θ=randθ(N),noise_scale=noise_scale,τ=τ)

## Example of Kuramoto simulation
N=2
# time-steps
t=10
Δt=0.2
steps=collect(0:Δt:t-Δt) # total number of steps calculated
# initialize natural frequencies
ω=rand(-2:0.0000001:2,N)
# couplings
σ=[1,1]
# initialize network
A=[0 1;1 0]
θs=Kuramodel.Kurasim(σ,ω,A,t,Δt,θ0=nothing,noise_scale=0.0,τ=nothing,seedval=nothing);