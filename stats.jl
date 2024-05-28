include("fcts.jl")

##########################################################
#                                                        #
# functions for calculating statistics of simulated data #
#                                                        #
##########################################################


"""
computes trait values from ψsystem and ψparameters
"""
function traits(sys::ψsystem,π::ψparameters)
    @unpack_ψsystem sys
    @unpack_ψparameters π

    # compute additive microbial values
    m = vec(sum(reshape((Mₘ.α' * Mₘ.κ)', Jₕ, K*Nₕ), dims=1))

    # compute additive genetic values
    g = (Mₕ.α' * Mₕ.κ)'

    # compute trait values
    z = g₀ .+ g .+ m

    return z
end

"""
computes allele frequencies from a (sparse) matrix
- `ℓ` = loci to compute allele frequencies for
- `ι` = individuals to compute allele frequencies of

returns a vector of length `length(ℓ)`
"""
function freqs(;ω::SparseMatrixCSC{Bool,Int64},ℓ::Vector{Int64} = Vector{Int64}(),ι::Vector{Int64} = Vector{Int64}())

    if (length(ℓ) == 0) & (length(ι) == 0)
        
        n = sum(ω,dims=2)
        return n

    elseif length(ℓ) == 0

        n = sum(ω[:,ι],dims=2)
        return n

    elseif length(ι) == 0

        n = sum(ω[ℓ,:],dims=2)
        return n

    else

        n = size(ω[:,ι])[2]
        p = sum(ω[ℓ,ι],dims=2)
        return p

    end

    print("apparent edge case. no conditions were met.")
    return 0

end

"""
computes allele frequencies from a (sparse) matrix
- `ℓ` = loci to compute allele frequencies for
- `ι` = individuals to compute allele frequencies of

returns a vector of length `length(ℓ)`
"""
# function freqs(;𝓜::SparseMatrixCSC{Bool,Int64},ℓ::UnitRange{Int64} = 0:-1,ι::UnitRange{Int64} = 0:-1)
#     freqs(𝓜=𝓜,ℓ=collect(ℓ),ι=collect(ι))
# end

"""
computes allele frequencies for local populations
- `ω` = `𝓛×n` mutation matrix, where `n` is total number individuals
- `u` = subpopulation size (must divide `n`)
"""
function local_freq_loop(ω::SparseMatrixCSC{Bool,Int64},U::Int64)
    
    # total number of individuals
    N = size(ω)[2]

    if mod(N,U) != 0
        print("supopulation size does not divide total population size")
        return
    end

    K = trunc(Int64,N/U)

    𝓛 = size(ω)[1]

    n = zeros(𝓛,K)

    for k in 1:K

        ι = (1:U) .+ (k-1)*U

        n[:,k] = freqs(ω=ω,ι=collect(ι))

    end

    return n

end

"""
computes local allele frequencies among hosts at each location
- returns two matrices of allele frequencies
    - one for neutral loci
    - one for causal loci
"""
function local_hfreqs(sys::ψsystem,π::ψparameters)
    @unpack_ψsystem sys
    @unpack_ψparameters π

    nⁿ = local_freq_loop(Mₕ.ν,Nₕ)

    nᶜ = local_freq_loop(Mₕ.κ,Nₕ)

    return nⁿ, nᶜ
end

"""
computes global allele frequencies among all host
individuals pooled together as a single population
- returns two vectors of allele frequencies
    - one for neutral loci
    - one for causal loci
"""
function global_hfreqs(sys::ψsystem)
    @unpack_ψsystem sys

    nⁿ = freqs(ω=Mₕ.ν)
    nᶜ = freqs(ω=Mₕ.κ)

    return nⁿ, nᶜ
end

"""
computes global allele frequencies among all host
individuals pooled together as a single population
- returns a vector of length `𝓛ₘ`
"""
function global_mfreqs(sys::ψsystem)
    
end

"""
computes local allele frequencies among hosts at each location
- returns a `𝓛ₘ×K` matrix of allele frequencies
"""
function local_mfreqs(sys::ψsystem,π::ψparameters)
    @unpack_ψsystem sys
    @unpack_ψparameters π

    nⁿ = local_freq_loop(Mₘ.ν,Jₕ)

    nᶜ = local_freq_loop(Mₘ.κ,Jₕ)

    return nⁿ, nᶜ
    
end

function sfs(p::Vector{Float64}; bins::Vector{Float64})

    

end

function local_adaptation()

    # measure by transplanting entire populations? since fitnesses are relative and thus depend on population state.

    # or maybe use "absolute" fitness? or even log-absolute fitness (so that baseline fitness coefficient has no effect)?
    #   /\ so far this has my vote

end