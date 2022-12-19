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

seedvalue=4; Random.seed!(seedvalue); N=2;
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
θs=Kuramodel.Kurasim(σ,ω,A,t,Δt,seedval=seedvalue)
θs

Nc = [4,3,2,3]
using Distributions
splits=get_splits(Nc)
ω = [1,2,3,2,4,5,6,3,3,7,8,9]
means=[]
for i in 1:length(splits)
   if i%2 == 1
       push!(means,mean(ω[splits[i]:splits[i+1]]))
   end
end
means