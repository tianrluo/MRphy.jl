
#= blochsim =#
export blochsim!, blochsim
"""
    blochsim!(M, B; T1=(Inf)u"s",T2=(Inf)u"s",γ=γ¹H,dt=(4e-6)u"s",doHist=false)
Old school 𝐵-effective magnetic field, `B`, based bloch simulation. Globally or
spin-wisely apply `B` over spins, `M`. `M` will be updated by the results.

*INPUTS*:
- `M::TypeND(Real, [2])` (nM, xyz): input spins' magnetizations.
- `B::Union{TypeND(B0D, [2,3]), Base.Generator}`:
  Global, (nT,xyz); Spin-wise, (nM,xyz,nT).
*KEYWORDS*:
- `T1 & T2 ::TypeND(T0D, [0,1])`: Global, (1,); Spin-wise, (nM,1).
- `γ::TypeND(Γ0D, [0,1])`: Global, (1,); Spin-wise, (nM, 1). gyro ratio
- `dt::T0D` (1,), simulation temporal step size, i.e., dwell time.
- `doHist::Bool`, whether to output spin history through out `B`.
*OUTPUTS*:
- `M::TypeND(Real, [2])` (nM, xyz): spins after applying `B`.
- `Mhst::TypeND(Real, [3])` (nM, xyz, nT): spins history during `B`.

See also: [`applyPulse`](@ref), [`blochsim`](@ref).

# Notes:
1. Not much sanity check inside this function, user is responsible for
   matching up the dimensions.
2. Put relax at the end of each time step may still be inaccurate, since
   physically spins relax continuously, this noise/nuance may worth study
   for applications like fingerprinting simulations, etc.
"""
function blochsim!(M ::TypeND(AbstractFloat, [2]),
                   B ::Base.Generator;
                   T1::TypeND(T0D, [0,1])=(Inf)u"s",
                   T2::TypeND(T0D, [0,1])=(Inf)u"s",
                   γ ::TypeND(Γ0D, [0,1])=γ¹H,
                   dt::T0D=(4e-6)u"s", doHist=false)

  # in unit, convert relaxations into losses/recovery per step
  E1, E2 = exp.(-dt./T1), exp.(-dt./T2)
  E1₋₁ = E1 .- 1

  nM = maximum([size(x,1) for x in (T1,T2,γ,first(B))])
  size(M,1) == 1 && (M = repeat(M, nM))

  Mhst = doHist ? zeros(size(M,1), 3, length(B)) : nothing
  Mi, M1 = M, similar(M) # Mi refers to the original array `M` in memory.

  # u, ϕ = Array{Float64}(undef, nM, 3), Array{Float64}(undef, nM, 1)

  if doHist
    for (t, b) in enumerate(B)
      u, ϕ = B2UΦ(b; γ=γ, dt=dt)
      # B2UΦ!(b, u, ϕ; γ=γ, dt=dt)
      any(ϕ.!=0) && @inbounds(UΦRot!(u, view(ϕ,:,1), M, M1))
      # relaxation
      M1[:,1:2] .*= E2
      M1[:,3]   .*= E1
      M1[:,3]   .-= E1₋₁
      @inbounds(Mhst[:,:,t] = M1)
      M, M1 = M1, M
    end
  else
    for b in B
      u, ϕ = B2UΦ(b; γ=γ, dt=dt)
      # B2UΦ!(b, u, ϕ; γ=γ, dt=dt)
      any(ϕ.!=0) && @inbounds(UΦRot!(u, view(ϕ,:,1), M, M1))
      # relaxation
      M1[:,1:2] .*= E2
      M1[:,3]   .*= E1
      M1[:,3]   .-= E1₋₁
      M, M1 = M1, M
    end
  end
  M === Mi || (Mi .= M) # if `M` doesn't point to the input array, update.

  return (M=M, Mhst=Mhst)
end

function blochsim!(M::TypeND(AbstractFloat, [2]), B::TypeND(B0D, [2,3]); kw...)
  if size(B,3) == 1 && size(B,1) != size(M,1)
    B = permutedims(B[:,:,:], [3,2,1])  # best practice?
    println("B not being spin-specific, assuming global")
  end
  return blochsim!(M, @inbounds(view(B,:,:,t) for t in axes(B,3)); kw...)
end

"""
    blochsim!(M, A, B)
Hargreave's 𝐴/𝐵, mat/vec, based bloch simulation. Globally or spin-wisely apply
matrix `A` and vector `B` over spins, `M`, described in doi:10.1002/mrm.1170

*INPUTS*:
- `M::TypeND(Real, [2])` (nM, xyz): input spins' magnetizations.
- `A::TypeND(AbstractFloat,[3])` (nM, 3,3), `A[iM,:,:]` is the `iM`-th 𝐴.
- `B::TypeND(AbstractFloat,[2])` (nM, 3), `B[iM,:]` is the `iM`-th 𝐵.
*OUTPUTS*:
- `M::TypeND(Real, [2])` (nM, xyz): spins after applying `B`.
"""
blochsim!(M, A, B) = M .= blochsim(M, A, B)

"""
    blochsim(M, B; T1, T2, γ, dt, doHist)
Same as `blochsim!(M, B; T1,T2,γ,dt,doHist)`, `M` will not be updated.

See also: [`blochsim!`](@ref)
"""
blochsim(M, B; kw...) = blochsim!(copy(M), B; kw...)

"""
    blochsim(M, A, B)
Same as `blochsim(M, A, B)`, `M` will not be updated.
"""
blochsim(M::TypeND(AbstractFloat,[2]),
         A::TypeND(AbstractFloat,[3]),
         B::TypeND(AbstractFloat,[2])) =
  sum(A.*permutedims(M[:,:,:],(1,3,2));dims=3) .+ B

# No inplace operation for `M::TypeND(Integer,[2])`.
blochsim(M::TypeND(Integer,[2]), a...; kw...) = blochsim(float(M), a...; kw...)

#= freePrec =#
export freePrec!, freePrec
"""
    freePrec!(M, t; Δf=0u"Hz", T1=(Inf)u"s", T2=(Inf)u"s")
Spins, `M`, free precess by time `t`. `M` will be updated by the results.

*INPUTS*:
- `M::TypeND(Real, [2])` (nM, xyz): input spins' magnetizations.
- `t::T0D` (1,): duration of free precession.
*KEYWORDS*:
- `T1 & T2 ::TypeND(T0D, [0,1])`: Global, (1,); Spin-wise, (nM,1).
*OUTPUTS*:
- `M::TypeND(Real, [2])` (nM, xyz): output spins' magnetizations.

See also: [`freePrec`](@ref).
"""
function freePrec!(M ::TypeND(AbstractFloat,[2]),
                   t ::T0D;
                   Δf::TypeND(F0D,[0,1])=0u"Hz",
                   T1::TypeND(T0D,[0,1])=(Inf)u"s",
                   T2::TypeND(T0D,[0,1])=(Inf)u"s")

  E1, E2 = exp.(-t./T1), exp.(-t./T2)

  M[:,1:2] .*= E2
  M[:,3]   .*= E1
  M[:,3]   .+= 1 .- E1

  eΔθ = exp.(-1im*2π*Δf*t)

  M[:,1:2] .= ((view(M,:,1)+1im*view(M,:,2)).*eΔθ |> x->[real(x) imag(x)])

  return M
end

"""
    freePrec(M, t; Δf, T1, T2)
Same as `freePrec!(M, t; Δf, T1, T2)`, `M` will not be updated.

See also: [`freePrec!`](@ref).
"""
freePrec(M, t; kw...) = freePrec!(copy(M), t; kw...)

