using Parameters, LinearAlgebra, Random, Distributions, DataFrames, CSV, Documenter, DocumenterTools, SparseArrays, ProgressMeter

"""
mutation tables and probabilities
- `μⁿ` = probability of neutral mutation
- `μᶜ` = probability of  causal mutation
- `ν`  = table of neutral mutations
- `κ`  = table of causal  mutations
- `α`  = vector of additive genetic effects
- `σ`  = std dev of new additive genetic effects
"""
@with_kw mutable struct mutations

    μⁿ::Float64 = 1e-3  # neutral mutation probability
    μᶜ::Float64 = 1e-3  # causal  mutation probability

    # neutral mutation table 
    # for hosts: 𝓛ⁿₕ×KNₕ (note: \bscrL -> 𝓛)
    # for microbes: 𝓛ⁿₘ×K(NₕJₕ+Jₑ)
    ν::SparseMatrixCSC{Bool,Int64}=sparse(Matrix{Bool}(undef,0,0))

    # causal mutation table
    # for hosts: 𝓛ᶜₕ×KNₕ (note: we use 𝓛 instead of S because these sites not necessarily polymorphic)
    # for host-associated microbes: 𝓛ᶜₘ×KNₕJₕ
    # for environmental microbes: 𝓛ᶜₘ×KJₑ
    κ::SparseMatrixCSC{Bool,Int64}=sparse(Matrix{Bool}(undef,0,0))

    # additive genetic fx
    α::Vector{Float64} = Vector{Float64}() # length = 𝓛ᶜₛ for s=h,m

    # std dev of new additive genetic fx
    σ::Float64 = 0.01

    # expected pairwise genetic distance of a propagule
    # to regional microbial ancestor
    D::Float64 = 100

    # fraction propagule mutations expected to be causal
    χ::Float64 = 0.1

    # additive genetic values calculated as κₕ*αₕ .+ κₘ*αₘ

    # z = g₀ .+ reshape(κₕ*αₕ .+ κₘ*αₘ, K, Nₕ)

end

"""
additive value vector, mutation probability, and mutational variance
- `μ` = probability of neutral mutation
- `α` = vector of additive value
- `σ` = std dev of mutation
"""
@with_kw mutable struct additives

    μ::Float64 = 1e-3  # mutation probability

    # additive values
    α::Vector{Float64} = Vector{Float64}() # length = 𝓛ᶜₛ for s=h,m

    # std dev of mutation
    σ::Float64 = 1.0

end

"""
model parameters

system parameters
- `Γₘ`  = microbe generations per host generation
- `s`   = host selection strength
- `Nₕ`  = number of hosts
- `Jₑ`  = number of microbes in environment
- `Jₕ`  = number of microbes in a host
- `K`   = number of locations

simulation parameters
- `Γₕ` = number of host generations to run model
- `R` = number of replicates
- `B` = host gens burnin period
- `fldr` = folder to save data in
- `fname` = file name for saving data and parameters

simulation variable
- `k` = replicate or parameter combination index

note: mutation parameters are stored in struct `mutations`
"""
@with_kw mutable struct parameters

    # sizes
    Nₕ::Int64   = 100   # number of hosts per location
    Jₕ::Int64	= 10^3  # number of microbes in a host
    Jₑ::Int64   = 10^3  # number of microbes in local environment
    K::Int64    = 10    # number of locations

    # movement probabilities
    dₕ::Float64 = 1e-3  # host dispersal probability
    dₘ::Float64 = 1e-12 # microbe dispersal probability
    ʂ::Float64  = 1e-6  # env shedding probability
    ε::Float64 = 1e-6   # env acquisition probability
    t::Float64 = 1e-6   # social transmission probability
    p::Float64 = 1e-12  # distal microbe propagule probability
    
    # host selection parameters
    s̄::Float64 = 0.0    # expected selection strength
    θ̄::Float64 = 0.0    # expected host optima
    vₛ::Float64 = 0.01  # variance of selection strengths
    vₜ::Float64 = 1.0   # variance of host optima
    s::Vector{Float64} = abs.(rand(Normal(s̄,√vₛ),K)) # local host selection strengths
    θ::Vector{Float64} =      rand(Normal(θ̄,√vₜ),K)  # local host optima

    # system parameters
    Γₘ::Int64 = 30          # microbe generations per host generation
    Γₕ::Int64 = 50          # number of host generations to run model
    τ::Int64 = 10           # number of host generations between purge check
    R::Int64 = 100          # number of replicates
    B::Int64 = 0            # host gens burnin period
    fname::String = "000"   # file name for saving data/parameters
    fldr::String = "tst/"   # folder to save data in
    verbal::Bool = false    # whether to print host generation each purge check
    
    # simulation variable
    r::Int64 = 1            # current replicate

end

"""
model parameters
- special case where the environment is ignored

system parameters
- `Γₘ`  = microbe generations per host generation
- `s`   = host selection strength
- `Nₕ`  = number of hosts
- `Jₕ`  = number of microbes in a host
- `K`   = number of locations

simulation parameters
- `Γₕ` = number of host generations to run model
- `R` = number of replicates
- `B` = host gens burnin period
- `fldr` = folder to save data in
- `fname` = file name for saving data and parameters

simulation variable
- `k` = replicate or parameter combination index

note: mutation parameters are stored in struct `mutations`
"""
@with_kw mutable struct ψparameters

    # sizes
    Nₕ::Int64   = 100   # number of hosts per location
    Jₕ::Int64	= 10^3  # number of microbes in a host
    K::Int64    = 10    # number of locations

    # movement probabilities
    dₕ::Float64 = 1e-3  # host dispersal probability
    dₘ::Float64 = 1e-12 # microbe dispersal probability
    t::Float64 = 1e-6   # social transmission probability
    p::Float64 = 1e-12  # distal microbe propagule probability
    
    # host selection parameters
    s̄::Float64 = 0.0    # expected selection strength
    θ̄::Float64 = 0.0    # expected host optima
    vₛ::Float64 = 0.01  # variance of selection strengths
    vₜ::Float64 = 1.0   # variance of host optima
    s::Vector{Float64} = abs.(rand(Normal(s̄,√vₛ),K)) # local host selection strengths
    θ::Vector{Float64} =      rand(Normal(θ̄,√vₜ),K)  # local host optima

    # system parameters
    Γₘ::Int64 = 30          # microbe generations per host generation
    Γₕ::Int64 = 50          # number of host generations to run model
    τ::Int64 = 10           # number of host generations between purge check
    R::Int64 = 100          # number of replicates
    B::Int64 = 0            # host gens burnin period
    fname::String = "000"   # file name for saving data/parameters
    fldr::String = "tst/"   # folder to save data in
    verbal::Bool = false    # whether to print host generation each purge check

    # simulation variable
    r::Int64 = 1            # current replicate

end

"""
state of the system
- special case where the environment is ignored
"""
@with_kw mutable struct system

    # host mutations
    Mₕ = mutations()

    # host-associated microbe mutations
    Mₘ = mutations()

    # environmental microbe mutations
    Mₑ = mutations()

    # additive genetic accumulator
    # to accumulate additive effects of fixed loci
    # only used when p = 0 (so absence of propagule pressure)
    g₀ = 0

    # # host trait table (K×Nₕ)
    # z = Matrix{Float64}()

    # host generation number
    τₕ::Int64 = 0

    # microbe generation number
    τₘ::Int64 = 0

end

"""
state of the system
- special case that ignores environmental microbiomes
"""
@with_kw mutable struct ψsystem

    # host mutations
    Mₕ = mutations()

    # host-associated microbe mutations
    Mₘ = mutations()

    # additive genetic accumulator
    # only used when p = 0
    g₀ = 0

    # # host trait table (K×Nₕ)
    # z = Matrix{Float64}()

    # host generation number
    τₕ::Int64 = 0

    # microbe generation number
    τₘ::Int64 = 0

end

"""
state of the system
- special case that ignores environmental microbiomes and 
uses additive values instead of genomes
"""
@with_kw mutable struct ψΣsystem

    # host additive values
    Σₕ = additives()

    # host-associated microbial additive values
    Σₘ = additives()

    # host generation number
    τₕ::Int64 = 0

    # microbe generation number
    τₘ::Int64 = 0

end
