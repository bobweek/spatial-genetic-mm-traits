include("fcts.jl")

# ss, π = ψsys_init(Jₕ=10)

# simulate without environmental microbiomes
ss, π = ψsys_init(K=10,
                  Nₕ=1000,
                  Jₕ=100,
                  Γₘ=10,
                  Γₕ=20,
                  μⁿₕ=0.01,
                  μᶜₕ=0.01,
                  μⁿₘ=0.001,
                  μᶜₘ=0.001,
                  σₕ=1.0,
                  σₘ=0.1,
                  θ̄=10.0,
                  s̄=0.1,
                  vₜ=0.0,
                  vₛ=0.0)

ss, π = ψsys_init(K=2,
                  Nₕ=10000,
                  Jₕ=1,
                  Γₘ=0,
                  Γₕ=50,
                  μⁿₕ=0.001,
                  μᶜₕ=0.1,
                  μⁿₘ=0.0,
                  μᶜₘ=0.0,
                  σₕ=1.0,
                  σₘ=0.0,
                  θ̄=10.0,
                  s̄=0.01,
                  vₜ=0.0,
                  vₛ=0.0,
                  dₕ=0.9)

π.θ = [-10.0,10.0]

so = sim(ss,π)

histogram(traits(ss,π))

# save simulated output
using JLD
@save "sys.jld" so
@save "pi.jld" π

# are the mutation matrices sparse?
sum(so.Mₕ.ν)/length(so.Mₕ.ν)
sum(so.Mₕ.κ)/length(so.Mₕ.κ)
sum(so.Mₘ.ν)/length(so.Mₘ.ν)
sum(so.Mₘ.κ)/length(so.Mₘ.κ)

π.θ

# using ProfileCanvas
# using Profile
# Profile.clear()
# @profile sim(ss,π)
# ProfileCanvas.view()


@unpack_ψsystem ss
@unpack_ψparameters π


using BenchmarkTools

@benchmark socials(ss,π)

@benchmark neutral_hmutations(so.Mₕ)

@benchmark hgen(ss,π)

@benchmark mgen(ss,π)
