using Kuramodel
using JLD2
using Random
using Distributions
using LinearAlgebra

"""
Function to chekc if, given the required macroscopic cpuolings `δ` and some community size,
the microscopic coupling `β` satisfies the condition.
"""
function coupling_condition(β, δ, M)
	
	condition = δ ./ β .< M
	
	for i in condition
        if !i
            return false
        end
    end
    return true
	
end

"""
```jldoctest
	function m_mat(size, m, β)
```
Creates a matrix of size `size`, and with `m` elements equal to `β`. \\
All other elements are set to `0.0`.
"""
function m_mat(size, m, β)
	out = zeros(size)
	if m > size[1]*size[2]
		throw(error("too many connections!"))
	end
	while count(x -> x == β, out) < m
		local to_add = [rand((1:size[1])), rand((1:size[2]))]
		if out[to_add[1],to_add[2]] == 0.0
			out[to_add[1],to_add[2]] = β
		end
	end
	return out
end

function micro_structure(δ, β, α, C, M)
	
	sp = get_splits(repeat([M], C))
	
	# find number of connections required
	m = ones(length(δ))
	m[1] = ceil(Int64, δ[1] * M / β)
	m[2] = ceil(Int64, δ[2] * M / β)
	m[3] = ceil(Int64, δ[3] * M / β)
	
	# intra-community connectivity
	microA = zeros(M*C, M*C)
	for i in 1:lastindex(sp)
		if i % 2 == 1
			local mat = ones(M,M).*α
			microA[sp[i]:sp[i+1], sp[i]:sp[i+1]] = mat
		end
	end

	# inter-community connectivity
	mmat_ab = m_mat((M,M), m[1], β)
	mmat_bc = m_mat((M,M), m[2], β)
	mmat_cd = m_mat((M,M), m[3], β)

	microA[sp[1]:sp[2], sp[3]:sp[4]] = mmat_ab
	microA[sp[3]:sp[4], sp[1]:sp[2]] = mmat_ab'
	
	microA[sp[3]:sp[4], sp[5]:sp[6]] = mmat_bc
	microA[sp[5]:sp[6], sp[3]:sp[4]] = mmat_bc'
	
	microA[sp[5]:sp[6], sp[7]:sp[8]] = mmat_cd
	microA[sp[7]:sp[8], sp[5]:sp[6]] = mmat_cd'

	return microA
end

######### fixed parameters #########

# seedval
seedval = 123
Random.seed!(seedval)

# size of the system
# **Convert array job input to a number**
M = parse(Int, ARGS[1])  # oscillators per community 
# M = 10
C = 4  # number of communities
Nc = repeat([M], C)  # utils
sp = get_splits(Nc)  # utils

# system fixed parameters
macroω = [-2.0, -1.0, 1.0, 2.0]  # natural frequencies
x_desired = [π/8, π/6, π/4]  # desired phase differences
δ = [5.226251859505506, 6, 2 * sqrt(2)]  # macroscopic couplings δ

# simulation settings
t_f = 10^2  # simulation time
Δt = 0.01  # time-step size
noise = 0.1  # noise value
steps = collect(0+Δt : Δt : t_f)  # collect steps
tot_steps = length(steps)  # total number of steps

# frequencies, each community has a different average natural frequency
microω = []           # Normal(i, 2.0), M) or Uniform(i - 1.0, i + 1.0)
for i in macroω
    global microω = vcat(microω, rand(Uniform(i - 1.0, i + 1.0), M))
end

# initial random phases
microθ = rand(Uniform(-2, 2), C*M)

# coupling normalization # TODO: remove this requirement
microσ = repeat([1], C*M)

####################################

#### parameter space to explore ####

# α already normalized
αs = collect(1:2:100) ./ M
Np = lastindex(αs)

# ratios: β =  α * r
r = collect(1:20/50:20+(40/50)) .^-1

# β already normalized:
βs = zeros(Np, Np)
for i in 1:Np
    βs[i, :] = αs .* r[i]
end

# conditions test
check_conditions_mat = zeros(Np, Np)
for i in 1:Np
    for j in 1:Np

        # if α and β haven't been normalized
        # β = β / M
        
        # check condtion
        if coupling_condition(βs[i,j], δ, M)
            check_conditions_mat[i, j] = 1
        else
            check_conditions_mat[i, j] = 0
        end

    end
end

results = zeros(Np, Np)

for i in 1:Np
    for j in 1:Np

        if check_conditions_mat[i, j] == 1.0

            α = αs[j]
            β = βs[i, j]

            microA = micro_structure(δ, β, α, C, M)
                
            ##############
            # simulation #
            ##############

            microKura = Kura_obj(microσ, microω, microA, microθ)  # init Kura object

            microθs = zeros(tot_steps, C*M)  # preallocate phases storage

            for i in 1:tot_steps
                Kura_step(microKura, Δt, noise_scale = 0.0)
                microθs[i, :] = microKura.θ
            end

            # find average phases in the last 100 steps for each community
            local temp_vec = []
            for i in 1:lastindex(sp)
                if i % 2 == 1
                    append!(temp_vec, mean(microθs[end-100:end, sp[i]:sp[i+1]]))
                end
            end
            
            # calculate simulation phase difference
            x_simulation = abs.([temp_vec[2] - temp_vec[1], temp_vec[3] - temp_vec[2], temp_vec[4] - temp_vec[3]])

            # calculate distance to desired phase differences
            results[i,j] = norm(x_desired - x_simulation, 2)
        end
    end
end

####################################

############ save data #############

# set path MANUALLY
filepath = "/its/home/ec627/data/StructuralControl/CouplingRatio/test/"

fileroot = "CouplingRatio"
fileseed = "Seed" * string(seedval)
fileM = "M" * string(M)

# ratio to add MANUALLY
fileratio = "Ratio" * "1to20"
# alpha to add MANUALLY
filealpha = "Alpha" * "1to100"

filename = filepath * fileroot * fileseed * fileM * fileratio * filealpha * ".jld2"

save_object(filename, results)
# save_object("LOCAL-M10.jld2", results)

####################################
