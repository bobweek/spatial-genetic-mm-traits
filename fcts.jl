include("structs.jl")

############################
#                          #
# functions for simulating #
#                          #
############################

"""
draw neutral microbe mutations
- `M₁` = focal mutation table (env or ham)
- `M₂` = focal mutation table (ham or env)

the operation of mutation appends a matrix containing _neutral_ mutations
onto the focal mutation table `M₁.ν`, and appends an associated empty matrix onto the
non-focal mutation table `M₂.ν`. because the operation is symmetric for environmental
microbiomes (env) and host-associated microbiomes (ham), the identities of the focal
and non-focal mutation mutation tables are left vague.
"""
function neutral_mmutations(M₁::mutations,M₂::mutations)
    @unpack_mutations M₁

    # draw mutating individuals
    n = size(ν)[2]
    𝓜 = rand(Binomial(n,μⁿ))
    mutants = sample(1:n,𝓜,replace=false)

    # append new loci to non-focal table M₂.ν
    M₂.ν = sparse_vcat(M₂.ν,spzeros(Bool,(𝓜,n)))
    dropzeros!(M₂.ν)
    
    # create matrix of novel mutations
    νᵥ = spzeros(Bool,(𝓜,n))
    νᵥ[:,mutants] = spdiagm(0=>fill(true,𝓜))

    # append new mutations to focal neutral mutation table
    ν = sparse_vcat(ν,νᵥ)
    dropzeros!(ν)

    @pack_mutations! M₁
    return M₁, M₂
end

"""
draw causal microbe mutations
- `M₁` = focal mutation table (env or ham)
- `M₂` = focal mutation table (ham or env)

the operation of mutation appends a matrix containing _causal_ mutations
onto the focal mutation table `M₁.ν`, and appends an associated empty matrix onto the
non-focal mutation table `M₂.ν`. because the operation is symmetric for environmental
microbiomes (env) and host-associated microbiomes (ham), the identities of the focal
and non-focal mutation mutation tables are left vague.
"""
function causal_mmutations(M₁::mutations,M₂::mutations)
    @unpack_mutations M₁

    # draw mutating individuals
    n = size(κ)[2]
    𝓜 = rand(Binomial(n,μᶜ))
    mutants = sample(1:n,𝓜,replace=false)

    # append new loci to non-focal table M₂.κ
    M₂.κ = sparse_vcat(M₂.κ,spzeros(Bool,(𝓜,n)))
    dropzeros!(M₂.κ)

    # create matrix of novel mutations
    κᵥ = spzeros(Bool,(𝓜,n))
    κᵥ[:,mutants] = spdiagm(0=>fill(true,𝓜))

    # append new mutations to focal causal mutation table
    κ = sparse_vcat(κ,κᵥ)
    dropzeros!(κ)

    # append new effect sizes
    αₚ = rand(Normal(0,σ),𝓜)
    α = vcat(α,αₚ)          # to focal mutations
    M₂.α = vcat(M₂.α,αₚ)    # to non-focal mutations

    @pack_mutations! M₁
    return M₁, M₂
end

"""
draw microbial mutations for either ham or env microbes
"""
function mmutate(M₁::mutations,M₂::mutations)

    M₁, M₂ = neutral_mmutations(M₁,M₂)
    M₁, M₂ = causal_mmutations(M₁,M₂)
   
    return M₁, M₂

end

"""
draw all microbial mutations
"""
function mmutates(sys::system)
    @unpack_system sys

    Mₘ, Mₑ = mmutate(Mₘ,Mₑ)
    Mₑ, Mₘ = mmutate(Mₑ,Mₘ)

    @pack_system! sys
    return sys
end

"""
draw neutral host mutations
"""
function neutral_hmutations(M::mutations)
    @unpack_mutations M

    # draw mutating individuals
    n = size(ν)[2]
    𝓜 = rand(Binomial(n,μⁿ))
    mutants = sample(1:n,𝓜,replace=false)

    # create matrix of novel mutations
    νᵥ = spzeros(Bool,(𝓜,n))
    νᵥ[:,mutants] = spdiagm(0=>fill(true,𝓜))

    # append new mutations to neutral mutation table
    ν = sparse_vcat(ν,νᵥ)
    dropzeros!(ν)

    @pack_mutations! M
    return M
end

"""
draw causal host mutations
"""
function causal_hmutations(M::mutations)
    @unpack_mutations M

    # draw mutating individuals
    n = size(κ)[2]
    𝓜 = rand(Binomial(n,μᶜ))
    mutants = sample(1:n,𝓜,replace=false)

    # create matrix of novel mutations
    κᵥ = spzeros(Bool,(𝓜,n))
    κᵥ[:,mutants] = spdiagm(0=>fill(true,𝓜))

    # append new mutations to causal mutation table
    κ = sparse_vcat(κ,κᵥ)
    dropzeros!(κ)

    # append new effect size
    α = vcat(α,rand(Normal(0,σ),𝓜))

    @pack_mutations! M
    return M
end

"""
draw host mutations
"""
function hmutate(sys::system)
    @unpack_system sys

    neutral_hmutations(Mₕ)
    causal_hmutations(Mₕ)
   
    @pack_system! sys
    return sys
end

"""
a microbial dispersal event
- `i` = dispersing microbe
"""
function mdisp(M::mutations,i::Int64,Jₑ::Int64)
    @unpack_mutations M
    
    # location index of focal individual
    k = ceil(Int64,i/Jₑ)

    # pick rnd individual in different location
    before = 1:(Jₑ*(k-1))
    after  = (Jₑ*k+1):size(ν)[2]
    j = sample(vcat(before,after))

    # copy focal individuals genome onto rnd ind
    ν[:,j] = ν[:,i]
    κ[:,j] = κ[:,i]

    dropzeros!(ν)
    dropzeros!(κ)
    
    @pack_mutations! M
    return M
end

"""
draw microbial dispersal events
- only environmental microbes
"""
function mdisps(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    if K>1 # dispersal only happens when there's multiple locations

        n = K*Jₑ
        𝓓 = rand(Binomial(n,dₘ))
        who = sample(1:n,𝓓,replace=false)
        for i in who
            mdisp(Mₑ,i,Jₑ)
        end

    end

    @pack_system! sys
    return sys
end

"""
a host dispersal event
- `i` = dispersing host
"""
function hdisp(Mₕ::mutations,Mₘ::mutations,i::Int64,Nₕ::Int64,Jₕ::Int64)
    
    # location index of focal host
    k = ceil(Int64,i/Nₕ)

    # pick rnd host in different location
    before = 1:(Nₕ*(k-1))
    after  = (Nₕ*k+1):size(Mₕ.ν)[2]
    j = sample(vcat(before,after))

    # copy focal hosts genome onto rnd host
    Mₕ.ν[:,j] = Mₕ.ν[:,i]
    Mₕ.κ[:,j] = Mₕ.κ[:,i]
    dropzeros!(Mₕ.ν)
    dropzeros!(Mₕ.κ)
    
    # get focal host microbiome indices
    f = (k-1)*Nₕ + 1 # first host in focal host pop
    𝔦 = (k-1)*Nₕ*Jₕ .+ (i-f)*Jₕ .+ (1:Jₕ)

    # got rnd host microbiome indices
    k = ceil(Int64,j/Nₕ) # population of rnd host
    f = (k-1)*Nₕ + 1     # first host in rnd host pop
    𝔧 = (k-1)*Nₕ*Jₕ .+ (j-f)*Jₕ .+ (1:Jₕ)

    # copy focal host mbiome onto rnd host mbiome
    Mₘ.ν[:,𝔧] = Mₘ.ν[:,𝔦]
    Mₘ.κ[:,𝔧] = Mₘ.κ[:,𝔦]
    dropzeros!(Mₘ.ν)
    dropzeros!(Mₘ.κ)
    
    return Mₕ, Mₘ
end

"""
draw host dispersal events
"""
function hdisps(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    if K>1 # dispersal only happens when there's multiple locations

        n = K*Nₕ
        𝓓 = rand(Binomial(n,dₕ))
        who = sample(1:n,𝓓,replace=false)
        for i in who
            Mₕ, Mₘ = hdisp(Mₕ,Mₘ,i,Nₕ,Jₕ)
        end

    end

    @pack_system! sys
    return sys
end

"""
draw host dispersal events
"""
function hdisps(sys::ψsystem,π::ψparameters)
    @unpack_ψsystem sys
    @unpack_ψparameters π

    if K>1 # dispersal only happens when there's multiple locations

        #
        # pick emitters
        #

        # get emitting host indices
        n = K*Nₕ
        𝓓 = rand(Binomial(n,dₕ))
        emₕ = sample(1:n,𝓓,replace=false,ordered=true)

        # get emitting host ham indices
        emₘ = Vector{Int64}() # target microbe indices
        for e in emₕ
            
            # microbe indices for eth host
            𝔦ₘ = (e-1)*Jₕ .+ (1:Jₕ)

            # append host target ham indices
            emₘ = vcat(emₘ,𝔦ₘ)

        end

        #
        # pick targets
        #

        tgₕ = Vector{Int64}() # target host indices
        tgₘ = Vector{Int64}() # target microbe indices

        for k in 1:K

            # host indices at location k
            𝔦ₖ = (k-1)*Nₕ .+ (1:Nₕ)

            # all other host indices
            𝔦ₒ = vcat(1:(𝔦ₖ[1]-1),(𝔦ₖ[Nₕ]+1):(K*Nₕ))

            # number of host emmitters at location k
            𝓔 = sum(emₕ .∈ [𝔦ₖ])

            # draw non-local host target indices
            𝓣 = sample(𝔦ₒ,𝓔,replace=false)

            # append non-local host targets
            tgₕ = vcat(tgₕ,𝓣)

            # get host target ham indices
            for 𝓽 in 𝓣

                # microbe indices for 𝓽th host
                𝔦ₘ = (𝓽-1)*Jₕ .+ (1:Jₕ)

                # append host target ham indices
                tgₘ = vcat(tgₘ,𝔦ₘ)

            end

        end

        #
        # copy mutations
        #

        # host mutations
        Mₕ.ν[:,tgₕ] = Mₕ.ν[:,emₕ]
        Mₕ.κ[:,tgₕ] = Mₕ.κ[:,emₕ]
        
        # microbe mutations
        Mₘ.ν[:,tgₘ] = Mₘ.ν[:,emₘ]
        Mₘ.κ[:,tgₘ] = Mₘ.κ[:,emₘ]

    end

    @pack_ψsystem! sys
    return sys
end

"""
draw host dispersal events
"""
function hdisps(sys::ψΣsystem,π::ψparameters)
    @unpack_ψΣsystem sys
    @unpack_ψparameters π

    if K>1 # dispersal only happens when there's multiple locations

        #
        # pick emitters
        #

        # get emitting host indices
        n = K*Nₕ
        𝓓 = rand(Binomial(n,dₕ))
        emₕ = sample(1:n,𝓓,replace=false,ordered=true)

        # get emitting host ham indices
        emₘ = Vector{Int64}() # target microbe indices
        for e in emₕ
            
            # microbe indices for eth host
            𝔦ₘ = (e-1)*Jₕ .+ (1:Jₕ)

            # append host target ham indices
            emₘ = vcat(emₘ,𝔦ₘ)

        end

        #
        # pick targets
        #

        tgₕ = Vector{Int64}() # target host indices
        tgₘ = Vector{Int64}() # target microbe indices

        for k in 1:K

            # host indices at location k
            𝔦ₖ = (k-1)*Nₕ .+ (1:Nₕ)

            # all other host indices
            𝔦ₒ = vcat(1:(𝔦ₖ[1]-1),(𝔦ₖ[Nₕ]+1):(K*Nₕ))

            # number of host emmitters at location k
            𝓔 = sum(emₕ .∈ [𝔦ₖ])

            # draw non-local host target indices
            𝓣 = sample(𝔦ₒ,𝓔,replace=false)

            # append non-local host targets
            tgₕ = vcat(tgₕ,𝓣)

            # get host target ham indices
            for 𝓽 in 𝓣

                # microbe indices for 𝓽th host
                𝔦ₘ = (𝓽-1)*Jₕ .+ (1:Jₕ)

                # append host target ham indices
                tgₘ = vcat(tgₘ,𝔦ₘ)

            end

        end

        #
        # copy additive values
        #

        # host additive values
        Σₕ.α[:,tgₕ] = Σₕ.α[:,emₕ]
        
        # microbe additive values
        Σₘ.α[:,tgₘ] = Σₘ.α[:,emₘ]

    end

    @pack_ψΣsystem! sys
    return sys
end

"""
a shedding event
- `i` = individual being shed
"""
function shed(sys::system,π::parameters,i::Int64)
    @unpack_system sys
    @unpack_parameters π

    # compute location
    k = ceil(Int64,i/(Nₕ*Jₕ))

    # draw microbe from local environment
    j = sample(1:Jₑ) + (k-1)*Jₑ
    # this should be benchmarked against rand(DiscreteUniform(((k-1)*Jₑ+1),(k*Jₑ)))

    # replace environmental microbe mutations with shedded microbe
    Mₑ.ν[:,j] = Mₘ.ν[:,i]
    Mₑ.κ[:,j] = Mₘ.κ[:,i]
    dropzeros!(Mₑ.ν)
    dropzeros!(Mₑ.κ)

    @pack_system! sys
    return sys
end

"""
draw shedding events
"""
function sheds(sys::system,π::parameters)
    @unpack_parameters π

    n = K*Nₕ*Jₕ
    𝓢 = rand(Binomial(n,ʂ))
    who = sample(1:n,𝓢,replace=false)
    for i in who
        shed(sys,π,i)
    end

    return sys
end

"""
a environmental acquisition event
"""
function acquire(sys::system,π::parameters,i::Int64)
    @unpack_system sys
    @unpack_parameters π

    # compute location
    k = ceil(Int64,i/(Nₕ*Jₕ))

    # draw microbe from local environment
    j = sample(((k-1)*Jₑ+1):(k*Jₑ)) # this should be benchmarked against rand(DiscreteUniform(((k-1)*Jₑ+1),(k*Jₑ)))

    # replace host microbe mutations with acquired microbe
    Mₘ.ν[:,i] = Mₑ.ν[:,j]
    Mₘ.κ[:,i] = Mₑ.κ[:,j]
    
    @pack_system! sys
    return sys
end

"""
draw environmental acquisition events
"""
function acquires(sys::system,π::parameters)
    @unpack_parameters π

    n = K*Nₕ*Jₕ
    𝓔 = rand(Binomial(n,ε))
    who = sample(1:n,𝓔,replace=false)
    for i in who
        acquire(sys,π,i)
    end

    return sys
end

"""
a social transmission event
- `i` = individual chosen to be transmitted
"""
function social(Mₘ::mutations,π::parameters,i::Int64)
    @unpack_mutations Mₘ
    @unpack_parameters π

    # compute location
    k = ceil(Int64,i/(Nₕ*Jₕ))

    # first crobe in pop k
    fₘ = (k-1)*Nₕ*Jₕ + 1

    # compute local host index
    h = ceil(Int64,(i-fₘ)/Jₕ)

    # compute index of first microbe in focal host
    ϕ = (k-1)*Nₕ*Jₕ + (h-1)*Jₕ

    # microbes in focal host
    hᵩ = ϕ:(ϕ+Jₕ)

    # draw microbe from any local host other than focal one
    j = sample(setdiff(fₘ:(Nₕ*Jₕ+fₘ-1),hᵩ))

    # transmission
    ν[:,j] = ν[:,i]
    κ[:,j] = κ[:,i]

    @pack_mutations! Mₘ
    return Mₘ
end

"""
a social transmission event
- `i` = individual chosen to be transmitted
- exactly the same for regular `system`
"""
function social(Mₘ::mutations,π::ψparameters,i::Int64)
    @unpack_mutations Mₘ
    @unpack_ψparameters π

    # compute location
    k = ceil(Int64,i/(Nₕ*Jₕ))

    # first crobe in pop k
    fₘ = (k-1)*Nₕ*Jₕ + 1

    # compute local host index
    h = ceil(Int64,(i-fₘ)/Jₕ)

    # compute index of first microbe in focal host
    ϕ = (k-1)*Nₕ*Jₕ + (h-1)*Jₕ

    # microbes in focal host
    hᵩ = ϕ:(ϕ+Jₕ)

    # draw microbe from any local host other than focal one
    j = sample(setdiff(fₘ:(Nₕ*Jₕ+fₘ-1),hᵩ))

    # transmission
    ν[:,j] = ν[:,i]
    κ[:,j] = κ[:,i]

    @pack_mutations! Mₘ
    return Mₘ
end

"""
draw social transmission events
"""
function socials(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    n = K*Nₕ*Jₕ
    𝓣 = rand(Binomial(n,t))
    who = sample(1:n,𝓣,replace=false)
    for i in who
        social(Mₘ,π,i)
    end

    @pack_system! sys
    return sys
end

"""
social transmission of ham microbes
"""
function socials(sys::ψsystem,π::ψparameters)
    @unpack_ψsystem sys
    @unpack_ψparameters π

    # pick emitters
    
    n = K*Nₕ*Jₕ
    𝓣 = rand(Binomial(n,t))
    em = sample(1:n,𝓣,replace=false,ordered=true)

    # pick targets

    tg = Vector{Int64}()
    for k in 1:K

        # first ham microbe index at location k
        fₘ = (k-1)*Nₕ*Jₕ + 1

        for i in 1:Nₕ

            # microbe indices for ith host at location k
            𝔦 = (k-1)*Nₕ*Jₕ .+ (i-1)*Jₕ .+ (1:Jₕ)

            # all other ham microbe indices at location k
            𝔦ₒ = vcat(fₘ:(𝔦[1]-1),(𝔦[Jₕ]+1):(Nₕ*Jₕ+fₘ-1))
            
            # number of emmitters in host i
            𝓔 = sum(em .∈ [𝔦])
            
            # append randomly drawn targets
            tg = vcat(tg,sample(𝔦ₒ,𝓔,replace=false))

        end

    end

    # copy emitter neutral mutations onto targets
    Mₘ.ν[:,tg] = Mₘ.ν[:,em]

    # copy emitter causal mutations onto targets
    Mₘ.κ[:,tg] = Mₘ.κ[:,em]

    @pack_ψsystem! sys
    return sys
end


"""
social transmission of ham microbes
"""
function socials(sys::ψΣsystem,π::ψparameters)
    @unpack_ψΣsystem sys
    @unpack_ψparameters π

    # pick emitters
    
    n = K*Nₕ*Jₕ
    𝓣 = rand(Binomial(n,t))
    em = sample(1:n,𝓣,replace=false,ordered=true)

    # pick targets

    tg = Vector{Int64}()
    for k in 1:K

        # first ham microbe index at location k
        fₘ = (k-1)*Nₕ*Jₕ + 1

        for i in 1:Nₕ

            # microbe indices for ith host at location k
            𝔦 = (k-1)*Nₕ*Jₕ .+ (i-1)*Jₕ .+ (1:Jₕ)

            # all other ham microbe indices at location k
            𝔦ₒ = vcat(fₘ:(𝔦[1]-1),(𝔦[Jₕ]+1):(Nₕ*Jₕ+fₘ-1))
            
            # number of emmitters in host i
            𝓔 = sum(em .∈ [𝔦])
            
            # append randomly drawn targets
            tg = vcat(tg,sample(𝔦ₒ,𝓔,replace=false))

        end

    end

    # copy additive values
    Σₘ.α[:,tg] = Σₘ.α[:,em]

    @pack_ψΣsystem! sys
    return sys
end

"""
a distal propagule event (env microbes only!)
- `i` = individual replaced by propagule
"""
function dstl_ppgl(Mₑ::mutations,Mₘ::mutations,i::Int64)
    @unpack_mutations Mₑ
    # a distant microbial propagule lands in the environment
    # and displaces microbe individual i

    # individual i's mutations are wiped clean
    ν[:,i] .= false
    κ[:,i] .= false
    dropzeros!(ν)
    dropzeros!(κ)

    # new propagule mutations are drawn
    𝓜  = rand(Poisson(D))
    𝓜ᶜ = rand(Binomial(𝓜,χ))
    𝓜ⁿ = 𝓜 - 𝓜ᶜ

    # new loci appended to host-associated mutation tables
    n = size(Mₘ.ν)[2]
    νₚ = sparse(fill(false,(𝓜ⁿ,n)))
    κₚ = sparse(fill(false,(𝓜ᶜ,n)))
    Mₘ.ν = sparse_vcat(Mₘ.ν,νₚ)
    Mₘ.κ = sparse_vcat(Mₘ.κ,κₚ)

    # new mutations appended to env mutation tables
    n = size(ν)[2]
    νₚ = sparse(fill(false,(𝓜ⁿ,n)))
    κₚ = sparse(fill(false,(𝓜ᶜ,n)))
    νₚ[:,i] .= true
    κₚ[:,i] .= true
    ν = sparse_vcat(ν,νₚ)
    κ = sparse_vcat(κ,κₚ)
    
    # new causal mutation effects are appended
    αₚ = rand(Normal(0,σ),𝓜ᶜ)
    α = vcat(α,αₚ)
    Mₘ.α = vcat(Mₘ.α,αₚ)

    @pack_mutations! Mₑ
    return Mₑ, Mₘ
end

"""
draw distal propagule events
"""
function dstl_ppgls(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    # draw propagules
    n = K*Jₑ
    𝓟 = rand(Binomial(n,p))
    who = sample(1:n,𝓟,replace=false)
    for i in who
        Mₑ, Mₘ = dstl_ppgl(Mₑ,Mₘ,i)        
    end

    @pack_system! sys
    return sys
end

"""
microbiome drift
- `M` = mutation table containing drifting microbiome
- `m` = indices of microbes in drifting microbiome
"""
function mdrift(M::mutations,m::UnitRange{Int64})
    @unpack_mutations M

    # number of microbes
    J = length(m)

    # draw parental indices
    𝔓 = sample(m,J)

    # inheritance of neutral mutations
    ν[:,m] = ν[:,𝔓]

    # inheritance of causal mutations
    κ[:,m] = κ[:,𝔓]

    @pack_mutations! M
    return M
end

"""
independent drift for each microbiome
"""
function mdrifts(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    for k in 1:K

        for h in 1:Nₕ

            # indices of local host microbes
            m = (1:Jₕ) .+ (k-1)*Nₕ*Jₕ .+ (h-1)*Jₕ

            # host-associated microbiome drift
            mdrift(Mₘ,m)

        end

        # indices of local environmental microbes
        m = (1:Jₑ) .+ (k-1)*Jₑ
        
        # environmental microbiome drift
        mdrift(Mₑ,m)

    end
    
    @pack_system! sys
    return sys
end

"""
independent drift for each microbiome
"""
function mdrifts(sys::ψsystem,π::ψparameters)
    @unpack_ψsystem sys
    @unpack_ψparameters π

    𝔓 = Vector{Int64}()
    for k in 1:K
        for h in 1:Nₕ

            # indices of local host microbes
            m = (1:Jₕ) .+ (k-1)*Nₕ*Jₕ .+ (h-1)*Jₕ

            # append random parental indices
            𝔓 = vcat(𝔓,sample(m,Jₕ))

        end
    end

    # inheritance of neutral mutations
    Mₘ.ν = Mₘ.ν[:,𝔓]

    # inheritance of causal mutations
    Mₘ.κ = Mₘ.κ[:,𝔓]
    
    @pack_ψsystem! sys
    return sys
end

"""
independent drift for each microbiome
"""
function mdrifts(sys::ψΣsystem,π::ψparameters)
    @unpack_ψΣsystem sys
    @unpack_ψparameters π

    𝔓 = Vector{Int64}()
    for k in 1:K
        for h in 1:Nₕ

            # indices of local host microbes
            m = (1:Jₕ) .+ (k-1)*Nₕ*Jₕ .+ (h-1)*Jₕ

            # append random parental indices
            𝔓 = vcat(𝔓,sample(m,Jₕ))

        end
    end

    # inheritance of additive values
    Σₘ.α = Σₘ.α[:,𝔓]

    @pack_ψΣsystem! sys
    return sys
end

"""
one microbial generation
"""
function mgen(sys::system,π::parameters)

    # microbiome drift
    mdrifts(sys,π)
    
    # host-associated microbe mutations
    mmutates(sys)

    # distal propagules land in the environment
    dstl_ppgls(sys,π)

    # microbial transmission via social contact
    socials(sys,π)

    # environmental microbial dispersal
    mdisps(sys,π)

    # host shedding of microbes into the environment
    sheds(sys,π)

    # host acquisition of microbes from the environment
    acquires(sys,π)

    # update microbe generation number
    sys.τₘ += 1

    return sys
end

"""
initialize system and parameters
"""
function sys_init(;
    K::Int64 = 3,
    Jₑ::Int64 = 10,
    Jₕ::Int64 = 10,
    Nₕ::Int64 = 10,
    dₕ::Float64 = 0.2,
    dₘ::Float64 = 0.2,
    ʂ::Float64 = 0.1,
    ε::Float64 = 0.1,
    t::Float64 = 0.1,
    p::Float64 = 0.1,
    μⁿₕ::Float64 = 0.2,
    μᶜₕ::Float64 = 0.2,
    μⁿₘ::Float64 = 0.2,
    μᶜₘ::Float64 = 0.2,
    σ::Float64 = 0.01,
    χ::Float64 = 0.5,
    Γₘ::Int64 = 2,
    Γₕ::Int64 = 10,
    τ::Int64 = 5)

    π = parameters(
        K=K, 
        Jₑ=Jₑ, 
        Jₕ=Jₕ, 
        Nₕ=Nₕ, 
        dₕ=dₕ, 
        dₘ=dₘ, 
        ʂ=ʂ, 
        ε=ε, 
        t=t, 
        p=p,
        Γₘ=Γₘ,
        Γₕ=Γₕ,
        τ=τ)

    @unpack_parameters π

    Mₕ = mutations()
    Mₘ = mutations()
    Mₑ = mutations()    

    Mₕ.μⁿ = μⁿₕ
    Mₕ.μᶜ = μᶜₕ
    Mₘ.μⁿ = μⁿₘ
    Mₘ.μᶜ = μᶜₘ
    Mₑ.μⁿ = μⁿₘ
    Mₑ.μᶜ = μᶜₘ
    
    Mₕ.σ = σ
    Mₘ.σ = σ
    Mₑ.σ = σ

    Mₕ.χ = χ
    Mₘ.χ = χ
    Mₑ.χ = χ
    
    # microbe mutations

    Mₘ.ν = sparse(fill(false,(1,K*Nₕ*Jₕ)))
    Mₘ.κ = sparse(fill(false,(1,K*Nₕ*Jₕ)))
    Mₘ.α = zeros(1)

    Mₑ.ν = sparse(fill(false,(1,K*Jₑ)))
    Mₑ.κ = sparse(fill(false,(1,K*Jₑ)))
    Mₑ.α = zeros(1)

    Mₘ, Mₑ = mmutate(Mₘ,Mₑ)
    Mₑ, Mₘ = mmutate(Mₑ,Mₘ)

    # host mutations

    Mₕ.ν = sparse(fill(false,(1,K*Nₕ)))
    Mₕ.κ = sparse(fill(false,(1,K*Nₕ)))
    Mₕ.α = zeros(1)

    neutral_hmutations(Mₕ)
    causal_hmutations(Mₕ)
    
    sys = system(Mₕ=Mₕ,Mₘ=Mₘ,Mₑ=Mₑ)

    mpurge_check(sys)
    hpurge_check(sys)

    return sys, π
end

"""
checks for purged microbial mutations
"""
function mpurge_check(sys::system)
    @unpack_system sys

    # do this for microbes every host generation
    # and at end of sim

    #
    # neutral loci
    #

    # get indices for neutral loci
    𝓛 = eachindex(axes(Mₘ.ν,1))

    # will contain empty loci
    remove = Vector{Int64}(undef,0)

    # iterate over neutral loci
    for ℓ in 𝓛
        # if mutation has been purged
        if !any(Mₘ.ν[ℓ,:]) & !any(Mₑ.ν[ℓ,:])
            # then append locus to remove
            push!(remove,ℓ)
        end
    end

    # drop empty loci/columns from neutral mutation table
    Mₘ.ν = Mₘ.ν[setdiff(1:end,remove),:]
    Mₑ.ν = Mₑ.ν[setdiff(1:end,remove),:]
    dropzeros!(Mₘ.ν)
    dropzeros!(Mₑ.ν)

    #
    # causal loci
    #

    # get indices for causal loci
    𝓛 = eachindex(axes(Mₘ.κ,1))

    # will contain empty loci
    remove = Vector{Int64}(undef,0)

    # iterate over causal loci
    for ℓ in 𝓛
        # if mutation has been purged
        if !any(Mₘ.κ[ℓ,:]) & !any(Mₑ.κ[ℓ,:])
            # then append locus to remove
            push!(remove,ℓ)
        end
    end

    # drop empty loci/columns from neutral mutation table
    Mₘ.κ = Mₘ.κ[setdiff(1:end,remove),:]
    Mₑ.κ = Mₑ.κ[setdiff(1:end,remove),:]
    dropzeros!(Mₘ.κ)
    dropzeros!(Mₑ.κ)

    # remove associated entry from α
    deleteat!(Mₘ.α,remove)
    deleteat!(Mₑ.α,remove)

    @pack_system! sys
    return sys
end

"""
checks for purged host mutations
"""
function hpurge_check(sys::system)
    @unpack_system sys
    @unpack_mutations Mₕ

    # do this for hosts every ten generations
    # and at end of sim

    #
    # neutral loci
    #

    # get indices for neutral loci
    𝓛 = eachindex(axes(ν,1))

    # will contain empty loci
    remove = Vector{Int64}(undef,0)

    # iterate over neutral loci
    for ℓ in 𝓛
        # if mutation has been purged
        if !any(ν'[:,ℓ])
            # then append locus to remove
            push!(remove,ℓ)
        end
    end

    # drop empty loci/columns from neutral mutation table
    ν = ν'[:,setdiff(1:end,remove)]'

    #
    # causal loci
    #

    # get indices for causal loci
    𝓛 = eachindex(axes(κ,1))

    # will contain empty loci
    remove = Vector{Int64}(undef,0)

    # iterate over causal loci
    for ℓ in 𝓛
        # if mutation has been purged
        if !any(κ'[:,ℓ])
            # then append locus to remove
            push!(remove,ℓ)
        end
    end

    # drop empty loci/columns from neutral mutation table
    κ = κ'[:,setdiff(1:end,remove)]'

    # remove associated entry from α
    deleteat!(α,remove)

    @pack_mutations! Mₕ
    @pack_system! sys
    return sys
end

"""
checks for fixed mutations
- needs testing!!!
"""
function fixation_check(M::mutations)
    @unpack_mutations M

    # only do this when
    #   p = 0 because we need to
    #       compare w propagule genomes
    # and/or at the end of the sim
    #   (although it may useful to keep
    #       for reconstructing ancestry)
    # based on these considerations
    #   this fct may never get implemented

    #
    # neutral loci
    #

    # get indices for neutral loci
    𝓛 = eachindex(axes(ν,1))

    # will contain empty loci
    remove = Vector{Int64}(undef,0)

    # iterate over neutral loci
    for ℓ in 𝓛
        # if mutation has been fixed
        if all(ν'[:,ℓ])
            # then append locus to remove
            push!(remove,ℓ)
        end
    end

    # drop fixed loci/columns from neutral mutation table
    ν = ν'[:,setdiff(1:end,remove)]'

    #
    # causal loci
    #

    # get indices for causal loci
    𝓛 = eachindex(axes(κ,1))

    # will contain empty loci
    remove = Vector{Int64}(undef,0)

    # iterate over causal loci
    for ℓ in 𝓛
        # if mutation has been purged
        if all(κ'[:,ℓ])
            # then append locus to remove
            push!(remove,ℓ)
        end
    end

    # drop empty loci/columns from neutral mutation table
    κ = κ'[:,setdiff(1:end,remove)]'

    g₊ = sum(α[remove])

    # remove associated entry from α
    deleteat!(α,remove)
    
    @pack_mutations! M
    return M, g₊
end

function fixation_check(sys::system)
    @unpack_system sys

    Mₘ, g₊ = fixation_check(Mₘ)
    g₀ += g₊
    Mₑ, g₊ = fixation_check(Mₑ)
    g₀ += g₊
    Mₕ, g₊ = fixation_check(Mₕ)
    g₀ += g₊

    @pack_system! sys
    return sys
end

function fixation_check(sys::ψsystem)
    @unpack_ψsystem sys

    Mₘ, g₊ = fixation_check(Mₘ)
    g₀ += g₊
    Mₕ, g₊ = fixation_check(Mₕ)
    g₀ += g₊

    @pack_ψsystem! sys
    return sys
end

"""
computes expected relative fitnesses of host individuals
"""
function fitness(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    # compute additive microbial values
    m = vec(sum(reshape((Mₘ.α' * Mₘ.κ)', Jₕ, K*Nₕ), dims=1))

    # compute additive genetic values
    g = (Mₕ.α' * Mₕ.κ)'

    # compute trait values
    z = g₀ .+ g .+ m

    # compute relative fitness matrix
    # an Nₕ×K matrix with columns that sum to 1
    w = Matrix{Float64}(undef,Nₕ,K)
    for k in 1:K
        Z = z[((k-1)*Nₕ+1):(k*Nₕ)]
        w[:,k] = normalize(exp.(-s[k] .* (θ[k] .- Z).^2 ./ 2),1)
    end

    return w
end


"""
computes expected relative fitnesses of host individuals
- exact same as normal `system`
"""
function fitness(sys::ψsystem,π::ψparameters)
    @unpack_ψsystem sys
    @unpack_ψparameters π

    # compute additive microbial values
    m = vec(sum(reshape((Mₘ.α' * Mₘ.κ)', Jₕ, K*Nₕ), dims=1))

    # compute additive genetic values
    g = (Mₕ.α' * Mₕ.κ)'

    # compute trait values
    z = g₀ .+ g .+ m

    # compute relative fitness matrix
    # an Nₕ×K matrix with columns that sum to 1
    w = Matrix{Float64}(undef,Nₕ,K)
    for k in 1:K
        Z = z[((k-1)*Nₕ+1):(k*Nₕ)]
        w[:,k] = normalize(exp.(-s[k] .* (θ[k] .- Z).^2 ./ 2),1)
    end

    return w
end

"""
local host selection, drift, and reproduction
"""
function hseldift(Mₕ::mutations,Mₘ::mutations,k::Int64,w::Matrix{Float64})
    
    # number of local host parents
    Nₕ = length(w[:,k])

    # draw host parent realized fitnesses
    𝓦 = rand(Multinomial(Nₕ,w[:,k]))

    # index of first host in focal location
    fₕ = 1 + (k-1)*Nₕ

    # get relative host parent indices
    𝒫 = findall(x->x>0,𝓦)

    # first make vector length N of repeated parent indices
    ℌ = Vector{Int64}()
    for iᵣ in 𝒫

        # compute absolute host parent index
        iₐ = iᵣ .+ (fₕ-1)

        # append absolute host parent index 𝓦[iᵣ] times
        ℌ = vcat(ℌ,fill(iₐ,𝓦[iᵣ]))

    end

    # inheritance of neutral host mutations
    Mₕ.ν[:,fₕ:(fₕ+Nₕ-1)] = Mₕ.ν[:,ℌ]

    # inheritance of causal host mutations
    Mₕ.κ[:,fₕ:(fₕ+Nₕ-1)] = Mₕ.κ[:,ℌ]

    # do things for host mbiomes
    K  = floor(Int64,size(Mₕ.ν)[2]/Nₕ)
    Jₕ = floor(Int64,size(Mₘ.ν)[2]/(K*Nₕ))

    # microbe indices
    𝔐 = Vector{Int64}()
    for iᵣ in 𝒫

        # compute absolute host parent index
        iₐ = iᵣ .+ (fₕ-1)

        # microbe indices for host iₐ
        𝔦 = (iₐ-1)*Jₕ .+ (1:Jₕ)

        # copy these indices 𝓦[iᵣ] times
        #   where 𝓦[iᵣ] is the realized fitness of host iₐ
        𝔐 = vcat(𝔐,repeat(𝔦,outer=𝓦[iᵣ]))

    end

    # first local ham microbe index
    fₘ = (fₕ-1)*Jₕ + 1

    # inheritance of neutral microbe mutations
    Mₘ.ν[:,fₘ:(fₘ+Nₕ*Jₕ-1)] = Mₘ.ν[:,𝔐]

    # inheritance of causal microbe mutations
    Mₘ.κ[:,fₘ:(fₘ+Nₕ*Jₕ-1)] = Mₘ.κ[:,𝔐]

    return Mₕ, Mₘ
end

"""
independent host drift, selection, and reproduction for each location
"""
function hseldifts(sys::system,π::parameters)
    @unpack_system sys
    @unpack_parameters π

    w = fitness(sys,π)

    for k in 1:K

        # do local host drift
        hseldift(Mₕ,Mₘ,k,w)

    end

    @pack_system! sys
    return sys
end

"""
independent host drift and selection for each location
"""
function hseldifts(sys::ψsystem,π::ψparameters)
    @unpack_ψsystem sys
    @unpack_ψparameters π

    w = fitness(sys,π)

    for k in 1:K

        # do local host drift
        hseldift(Mₕ,Mₘ,k,w)

        # can be made faster by doing index calculations 
        # across all locations simultaneously

    end


    # # draw host parent realized fitnesses
    # 𝓦 = rand(Multinomial(Nₕ,w[:,k]))

    # # index of first host in focal location
    # fₕ = 1 + (k-1)*Nₕ

    # # get relative host parent indices
    # 𝒫 = findall(x->x>0,𝓦)

    # # first make vector length N of repeated parent indices
    # ℌ = Vector{Int64}()
    # for iᵣ in 𝒫

    #     # compute absolute host parent index
    #     iₐ = iᵣ .+ (fₕ-1)

    #     # append absolute host parent index 𝓦[iᵣ] times
    #     ℌ = vcat(ℌ,fill(iₐ,𝓦[iᵣ]))

    # end

    # # inheritance of neutral host mutations
    # Mₕ.ν[:,fₕ:(fₕ+Nₕ-1)] = Mₕ.ν[:,ℌ]

    # # inheritance of causal host mutations
    # Mₕ.κ[:,fₕ:(fₕ+Nₕ-1)] = Mₕ.κ[:,ℌ]

    # # do things for host mbiomes
    # K  = floor(Int64,size(Mₕ.ν)[2]/Nₕ)
    # Jₕ = floor(Int64,size(Mₘ.ν)[2]/(K*Nₕ))

    # # microbe indices
    # 𝔐 = Vector{Int64}()
    # for iᵣ in 𝒫

    #     # compute absolute host parent index
    #     iₐ = iᵣ .+ (fₕ-1)

    #     # microbe indices for host iₐ
    #     𝔦 = (iₐ-1)*Jₕ .+ (1:Jₕ)

    #     # copy these indices 𝓦[iᵣ] times
    #     #   where 𝓦[iᵣ] is the realized fitness of host iₐ
    #     𝔐 = vcat(𝔐,repeat(𝔦,outer=𝓦[iᵣ]))

    # end

    # # first local ham microbe index
    # fₘ = (fₕ-1)*Jₕ + 1

    # # inheritance of neutral microbe mutations
    # Mₘ.ν[:,fₘ:(fₘ+Nₕ*Jₕ-1)] = Mₘ.ν[:,𝔐]

    # # inheritance of causal microbe mutations
    # Mₘ.κ[:,fₘ:(fₘ+Nₕ*Jₕ-1)] = Mₘ.κ[:,𝔐]

    @pack_ψsystem! sys
    return sys
end

"""
performs one iteration of simulation model
- equals one host generation
- includes:
    - host dispersal
    - host selection
    - host drift
    - host mutation
    - host development
        - `Γₘ` microbe generations
"""
function hgen(sys::system,π::parameters)
    
    # host dispersal
    hdisps(sys,π)
    
    # host selection and drift
    hseldifts(sys,π)

    # host mutation
    hmutate(sys)

    # host development
    # want this to happen at the very end
    # so that this fct outputs data on adult hosts
    for i in 1:π.Γₘ
        sys = mgen(sys,π)
    end

    # check for purged mutations in microbe genomes
    mpurge_check(sys)

    # if its time to check hosts for purged mutations
    if mod(sys.τₕ,π.τ) == 0
        hpurge_check(sys)
        fixation_check(sys)
    end

    # update host generation number
    sys.τₕ += 1

    # reset microbe generation number
    sys.τₘ = 0

    return sys
end

"""
run a simulation for `Γₕ` host generations
"""
function sim(sys::system,π::parameters)
    for i in 1:π.Γₕ
        hgen(sys,π)
    end
    hpurge_check(sys)
    fixation_check(sys)
    return sys    
end


###############################################################################################
#                                                                                             #
#  the following are functs tailored specifically toward lineal and social transmission only  #
#                                                                                             #
###############################################################################################


"""
a _neutral_ microbial mutation event for hams only
- `i` = mutating individual
- `M` = mutation table (ham)

the operation of mutation inserts a locus vector containing the _neutral_ mutation
into the mutation table `M.ν` of host-associated microbes (ham).
"""
function neutral_mmutation(M::mutations,i::Int64)

    # make locus vector with new mutation
    νᵥ = sparse(fill(false,size(M.ν)[2])')
    νᵥ[i] = true

    # append locus vector
    M.ν = sparse_vcat(M.ν,νᵥ)
    dropzeros!(M.ν)

    return M
end

"""
a _causal_ microbial mutation event for hams only
- `i` = mutating individual
- `M` = mutation table (ham)

the operation of mutation inserts a locus vector containing the _causal_ mutation
into the mutation table `M.ν` of host-associated microbes (ham).
"""
function causal_mmutation(M::mutations,i::Int64)

    # make locus vector with new mutation
    κᵥ = sparse(fill(false,size(M.κ)[2])')
    κᵥ[i] = true

    # append locus vector
    M.κ = sparse_vcat(M.κ,κᵥ)
    dropzeros!(M.κ)

    # append new effect size
    αₚ = rand(Normal(0,M.σ))
    M.α = vcat(M.α,αₚ)

    return M
end


"""
draw neutral microbe mutations
"""
function neutral_mmutations(M)
    
    n = size(M.ν)[2]
    𝓜 = rand(Binomial(n,M.μⁿ))
    who = sample(1:n,𝓜,replace=false)
    for i in who
        M = neutral_mmutation(M,i)        
    end

    return M

end

"""
draw causal microbe mutations
"""
function causal_mmutations(M::mutations)

    n = size(M.κ)[2] # number of microbes in M₁
    𝓜 = rand(Binomial(n,M.μᶜ)) # number of mutations
    who = sample(1:n,𝓜,replace=false) # who mutates
    for i in who
        M = causal_mmutation(M,i)
    end

    return M

end

"""
draw microbial mutations
"""
function mmutate(M::mutations)

    M = neutral_mmutations(M)
    M = causal_mmutations(M)
   
    return M

end


"""
microbial mutations
"""
function mmutates(sys::ψsystem)
    @unpack_ψsystem sys

    Mₘ = neutral_mmutations(Mₘ)
    Mₘ = causal_mmutations(Mₘ)
   
    @pack_ψsystem! sys
    return sys
end

"""
draw host mutations
- exactly the same for regular `system`
"""
function hmutate(sys::ψsystem)
    @unpack_ψsystem sys

    neutral_hmutations(Mₕ)
    causal_hmutations(Mₕ)
   
    @pack_ψsystem! sys
    return sys
end

"""
a microbial dispersal event
- `i` = dispersing microbe
- needs to be adapted for inter-ham dispersal
- is being ignored as first pass
"""
function ψmdisp(M::mutations,i::Int64,Jₑ::Int64)
    @unpack_mutations M
    
    # location index of focal individual
    k = ceil(Int64,i/Jₑ)

    # pick rnd individual in different location
    before = 1:(Jₑ*(k-1))
    after  = (Jₑ*k+1):size(ν)[2]
    j = sample(vcat(before,after))

    # copy focal individuals genome onto rnd ind
    ν[:,j] = ν[:,i]
    κ[:,j] = κ[:,i]

    dropzeros!(ν)
    dropzeros!(κ)
    
    @pack_mutations! M
    return M
end

"""
draw microbial dispersal events
- needs to be adapted for inter-ham dispersal
- is being ignored as first pass
"""
function mdisps(sys::ψsystem,π::ψparameters)
    @unpack_ψsystem sys
    @unpack_ψparameters π

    if K>1 # dispersal only happens when there's multiple locations

        n = K*Nₕ*Jₘ
        𝓓 = rand(Binomial(n,dₘ))
        who = sample(1:n,𝓓,replace=false)
        for i in who
            ψmdisp(Mₘ,i,Jₘ)
        end

    end

    @pack_ψsystem! sys
    return sys
end


"""
a distal propagule event (ham microbes only!)
- `i` = individual replaced by propagule
- needs to be adapted for invading hams
- is being ignored as first pass
"""
function dstl_ppgl(Mₘ::mutations,i::Int64)

    # a distant microbial propagule lands on a host
    # and displaces microbe individual i

    # individual i's mutations are wiped clean
    Mₘ.ν[:,i] .= false
    Mₘ.κ[:,i] .= false
    dropzeros!(Mₘ.ν)
    dropzeros!(Mₘ.κ)

    # new propagule mutations are drawn
    𝓜  = rand(Poisson(Mₘ.D))
    𝓜ᶜ = rand(Binomial(𝓜,Mₘ.χ))
    𝓜ⁿ = 𝓜 - 𝓜ᶜ

    # new mutations appended to ham mutation tables
    n = size(Mₘ.ν)[2]
    νₚ = sparse(fill(false,(𝓜ⁿ,n)))
    κₚ = sparse(fill(false,(𝓜ᶜ,n)))
    νₚ[:,i] .= true
    κₚ[:,i] .= true
    Mₘ.ν = sparse_vcat(Mₘ.ν,νₚ)
    Mₘ.κ = sparse_vcat(Mₘ.κ,κₚ)
    
    # new causal mutation effects are appended
    αₚ = rand(Normal(0,Mₘ.σ),𝓜ᶜ)
    Mₘ.α = vcat(Mₘ.α,αₚ)

    return Mₘ
end

"""
draw distal propagule events
- needs to be adapted for invading hams
- is being ignored as first pass
"""
function dstl_ppgls(sys::ψsystem,π::ψparameters)
    @unpack_ψsystem sys
    @unpack_ψparameters π

    # draw propagules
    n = K*Nₕ*Jₕ
    𝓟 = rand(Binomial(n,p))
    who = sample(1:n,𝓟,replace=false)
    for i in who
        Mₘ = dstl_ppgl(Mₘ,i)        
    end

    @pack_ψsystem! sys
    return sys
end


"""
one microbial generation
"""
function mgen(sys::ψsystem,π::ψparameters)

    # microbiome drift
    mdrifts(sys,π)
    
    # host-associated microbe mutations
    mmutates(sys)

    # distal propagules land in the environment
    # dstl_ppgls(sys,π) # ignore as first pass

    # microbial transmission via social contact
    socials(sys,π)

    # ham microbial dispersal among locations
    # mdisps(sys,π) # ignore as first pass

    # update microbe generation number
    sys.τₘ += 1

    return sys
end

"""
initialize system and parameters
"""
function ψsys_init(;
    K::Int64 = 1,
    Jₕ::Int64 = 1,
    Nₕ::Int64 = 10,
    dₕ::Float64 = 0.2,
    dₘ::Float64 = 0.2,
    t::Float64 = 1e-3,
    p::Float64 = 0.1,
    μⁿₕ::Float64 = 0.2,
    μᶜₕ::Float64 = 0.2,
    μⁿₘ::Float64 = 0.1,
    μᶜₘ::Float64 = 0.1,
    σₕ::Float64 = 0.1,
    σₘ::Float64 = 0.1,
    χ::Float64 = 0.5,
    Γₘ::Int64 = 0,
    Γₕ::Int64 = 1,
    τ::Int64 = 5,
    θ̄::Float64 = 0.0,
    s̄::Float64 = 0.0,
    vₜ::Float64 = 0.0,
    vₛ::Float64 = 0.0,
    verbal::Bool = false)

    π = ψparameters(
        K=K,
        Jₕ=Jₕ, 
        Nₕ=Nₕ, 
        dₕ=dₕ,
        dₘ=dₘ, 
        t=t, 
        p=p,
        Γₘ=Γₘ,
        Γₕ=Γₕ,
        τ=τ,
        θ̄=θ̄,
        s̄=s̄,
        vₜ=vₜ,
        vₛ=vₛ,
        verbal=verbal)

    @unpack_ψparameters π

    Mₕ = mutations()
    Mₘ = mutations()

    Mₕ.μⁿ = μⁿₕ
    Mₕ.μᶜ = μᶜₕ
    Mₘ.μⁿ = μⁿₘ
    Mₘ.μᶜ = μᶜₘ
    
    Mₕ.σ = σₕ
    Mₘ.σ = σₘ

    Mₕ.χ = χ
    Mₘ.χ = χ
    
    # microbe mutations

    Mₘ.ν = sparse(fill(false,(1,K*Nₕ*Jₕ)))
    Mₘ.κ = sparse(fill(false,(1,K*Nₕ*Jₕ)))
    Mₘ.α = zeros(1)

    Mₘ = mmutate(Mₘ)

    # host mutations

    Mₕ.ν = sparse(fill(false,(1,K*Nₕ)))
    Mₕ.κ = sparse(fill(false,(1,K*Nₕ)))
    Mₕ.α = zeros(1)

    neutral_hmutations(Mₕ)
    causal_hmutations(Mₕ)
    
    # package system

    sys = ψsystem(Mₕ=Mₕ,Mₘ=Mₘ)

    mpurge_check(sys)
    hpurge_check(sys)

    return sys, π
end

function purge_check(M::mutations)
    @unpack_mutations M

    #
    # neutral loci
    #

    # get indices for neutral loci
    𝓛 = eachindex(axes(ν,1))

    # will contain indices of polymorphic loci
    poly = Vector{Int64}(undef,0)

    # iterate over neutral loci
    for ℓ in 𝓛
        # if locus is polymorphic
        if any(ν[ℓ,:])
            # then keep it
            push!(poly,ℓ)
        end
    end

    # drop empty loci/rows from neutral mutation table
    ν = ν[poly,:]
    dropzeros!(ν)

    #
    # causal loci
    #

    # get indices for causal loci
    𝓛 = eachindex(axes(κ,1))

    # will contain indices of polymorphic loci
    poly = Vector{Int64}(undef,0)

    # iterate over causal loci
    for ℓ in 𝓛
        # if locus is polymorphic
        if any(κ[ℓ,:])
            # then keep it
            push!(poly,ℓ)
        end
    end

    # drop empty loci/rows from causal mutation table
    κ = κ[poly,:]
    dropzeros!(κ)

    # remove entries from α associated with purged causal mutations
    α = α[poly]

    @pack_mutations! M
    return M
end

"""
checks for purged microbial mutations
"""
function mpurge_check(sys::ψsystem)
    @unpack_ψsystem sys
    Mₘ = purge_check(Mₘ)
    @pack_ψsystem! sys
    return sys
end

"""
checks for purged host mutations
"""
function hpurge_check(sys::ψsystem)
    @unpack_ψsystem sys
    Mₕ = purge_check(Mₕ)
    @pack_ψsystem! sys
    return sys
end


"""
performs one iteration of simulation model
- equals one host generation
- includes:
    - host dispersal
    - host selection
    - host drift
    - host mutation
    - host development
        - `Γₘ` microbe generations
"""
function hgen(sys::ψsystem,π::ψparameters)
    
    # host dispersal
    hdisps(sys,π)
    
    # host selection and drift
    hseldifts(sys,π)

    # host mutation
    hmutate(sys)

    # host development
    # want this to happen at the very end
    # so that this fct outputs data on adult hosts
    for i in 1:π.Γₘ
        sys = mgen(sys,π)
    end

    # if its time to check for purged/fixed mutations
    if mod(sys.τₕ,π.τ) == 0
        hpurge_check(sys)
        mpurge_check(sys)
        fixation_check(sys)
        if π.verbal
            print("\nHost generation:   ",sys.τₕ,"\n")
            print("Host neutral loci: ",size(sys.Mₕ.ν)[1],"\n")
            print("Host causal  loci: ",size(sys.Mₕ.κ)[1],"\n")
            print("Mcrb neutral loci: ",size(sys.Mₘ.ν)[1],"\n")
            print("Mcrb causal  loci: ",size(sys.Mₘ.κ)[1],"\n")
        end
    end

    # update host generation number
    sys.τₕ += 1

    # reset microbe generation number
    sys.τₘ = 0

    return sys
end

"""
run a simulation for `Γₕ` host generations
"""
function sim(sys::ψsystem,π::ψparameters)
    @showprogress dt=0.5 barglyphs=BarGlyphs("[=> ]") barlen=50 color=:magenta for i in 1:π.Γₕ
        hgen(sys,π)
    end
    hpurge_check(sys)
    fixation_check(sys)
    return sys    
end
