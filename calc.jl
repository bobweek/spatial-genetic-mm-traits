include("stats.jl")

using JLD # for loading system state
using Plots
theme(:dracula)

sstm = load("sys.jld", "so")
πi = load("pi.jld", "π")

histogram(traits(sstm,πi))

@unpack_ψsystem sstm
@unpack_ψparameters πi

hⁿ, hᶜ = global_hfreqs(sstm)

histogram(hⁿ)
histogram(hᶜ)

fsⁿ = Vector{Int64}()
fsᶜ = Vector{Int64}()
for i in 1:Nₕ
    push!(fsⁿ,sum(hⁿ.==i))
    push!(fsᶜ,sum(hᶜ.==i))
end
fsⁿ
fsᶜ

sstm.Mₕ.ν
sstm.Mₕ.κ

sstm.Mₘ.ν
sstm.Mₘ.κ


# garbage collection (don't run often, maybe ever?)
# GC.gc()

# provides information on variables
varinfo(sortby=:size)