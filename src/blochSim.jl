# include("mObjects.jl")

#= 𝐵-effective =#
# Methods for getting 𝐵-effective
export Pulse2B
"""
    B = Pulse2B(rf, gr, loc; Δf, b1Map, γ)
Turn rf, `rf`, and gradient, `gr`, into 𝐵-effective magnetic field.

*INPUTS*:
- `rf::TypeND(RF0D, [1,2])` (nT, (nCoil))
- `gr::TypeND(GR0D, [2])`   (nT, 3)
- `loc::TypeND(L0D, [2])`   (1,3) or (nM, 3), locations.
*KEYWORDS*:
- `Δf::TypeND(F0D, [0,1,2])` (1,)  or (nM,), off-resonance.
- `b1Map::TypeND(Union{Real,Complex},[0,1,2])` (1,) or (nM,(nCoils)),
   transmit sensitivity.
- `γ::TypeND(Γ0D, [0,1])` (1,)  or (nM,), gyro-ratio
*OUTPUS*:
- `B`, generator of `TypeND(B0D, [2])` (1,1,nT), 𝐵 field.
"""
function Pulse2B(rf   ::TypeND(RF0D, [1,2]),
                 gr   ::TypeND(GR0D, [2]),
                 loc  ::TypeND(L0D,  [2]) = [0 0 0]u"cm";
                 Δf   ::TypeND(F0D,  [0,1]) = 0u"Hz",
                 b1Map::TypeND(Union{Real,Complex}, [0,1,2,3]) = 1,
                 γ    ::TypeND(Γ0D,  [0,1]) = γ¹H,
                 doGen::Bool=false)

  nM = maximum(map(x->size(x,1), (loc, Δf, b1Map, γ)))

  Bxy_gen = b1Map==1 ?
    @inbounds(view(rf,t,:)       |> x->repeat([real(x) imag(x)], nM)
              for t in axes(rf,1)) :
    @inbounds(b1Map*view(rf,t,:) |> x->       [real(x) imag(x)]
              for t in axes(rf,1))

  Bz_gen = Δf==0u"Hz" ?
    @inbounds(loc*view(gr,t,:)          for t in axes(gr,1)) :
    @inbounds(loc*view(gr,t,:).+(Δf./γ) for t in axes(gr,1))

  B_gen = @inbounds([bxy bz] for (bxy, bz) in zip(Bxy_gen, Bz_gen))
  return B_gen
end

"""
    B = Pulse2B(pulse::Pulse, loc; Δf, b1Map, γ)
Turn struct `Pulse` into effective magnetic, 𝐵, field.
"""
Pulse2B(p::Pulse, loc; kw...) = Pulse2B(p.rf, p.gr, loc; kw...)

"""
    B = Pulse2B(pulse::Pulse, spa::AbstractSpinArray, loc; Δf, b1Map)
...with `γ=spa.γ`.
"""
Pulse2B(p::Pulse, spa::AbstractSpinArray, loc; kw...) =
  Pulse2B(p, loc; γ=spa.γ, kw...)

"""
    B = Pulse2B(pulse::Pulse, cb::AbstractSpinCube; b1Map)
...with `loc, Δf, γ = cb.loc, cb.Δf, cb.γ`.
"""
Pulse2B(p::Pulse, cb::AbstractSpinCube; kw...) =
  Pulse2B(p, cb.loc, cb.Δf; γ=cb.γ, kw...)

#= rotation axis/angle, U/Φ =#
export B2UΦ
"""
    B2UΦ(B::TypeND(B0D,[2,3]); γ::TypeND(Γ0D,[0,1]), dt::T0D=4e-6u"s")
Given 𝐵-effective, `B`, compute rotation axis/angle, `U`/`Φ`.

*INPUTS*:
- `B::TypeND(B0D, [2,3])` (1,3,nT) or (nM, 3, nT), 𝐵 field.
*KEYWORDS*:
- `γ::TypeND(Γ0D, [0,1])`: Global, (1,); Spin-wise, (nM, 1). gyro ratio
- `dt::T0D` (1,), simulation temporal step size, i.e., dwell time.
*OUTPUTS*:
- `U::TypeND(Real, [2,3])` (1,3,(nT)) or (nM,3,(nT)), axis.
- `Φ::TypeND(Real, [2,3])` (1,1,(nT)) or (nM,1,(nT)), angle.
"""
@inline function B2UΦ(B::TypeND(B0D, [2,3]);
                      γ::TypeND(Γ0D, [0,1]), dt::T0D=4e-6u"s")
  Bn   = sqrt.(sum(B.*B, dims=2))     # norm of B
  Bn⁻¹ = map(x-> x==0 ? 0 : -1/x, Bn) # negate to correct cross product dir

  U = isa(B[1].*Bn⁻¹[1], Real) ? B.*Bn⁻¹     : uconvert.(NoUnits, B.*Bn⁻¹)
  Φ = isa(Bn[1]*γ[1]*dt, Real) ? Bn.*γ*2π*dt : uconvert.(NoUnits, Bn.*γ*2π*dt)
  return U, Φ
end

export UΦRot
"""
    UΦRot
"""
@inline function UΦRot(U::TypeND(AbstractFloat,[2]),
                       Φ::TypeND(AbstractFloat,[1]),
                       V::TypeND(AbstractFloat,[2,3]))
  # en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
  # 𝑅 = 𝑐𝑜𝑠𝜃⋅𝐼 - (𝑐𝑜𝑠𝜃-1)⋅(𝐮𝐮ᵀ) + 𝑠𝑖𝑛𝜃⋅[𝐮]ₓ; 𝜃/𝐮, rotation angle/axis
  if any(Φ.!= 0)
    cΦ, sΦ = cos.(Φ), sin.(Φ)
    Res = cΦ.*V .+ (1 .- cΦ).*sum(U.*V,dims=2).*U .+
          sΦ.*hcat(-U[:,3].*V[:,2,:] .+U[:,2].*V[:,3,:],
                    U[:,3].*V[:,1,:] .-U[:,1].*V[:,3,:],
                   -U[:,2].*V[:,1,:] .+U[:,1].*V[:,2,:])
  else # no rotation needed
    Res = V
  end
  return Res
end

#= 𝐴, 𝐵 =#
# Methods for getting Hargreave's 𝐴/𝐵, mat/vec, defined in doi:10.1002/mrm.1170
export B2AB

function B2AB(B ::TypeND(B0D, [2,3]);
              T1::TypeND(T0D, [0,1])=(Inf)u"s",
              T2::TypeND(T0D, [0,1])=(Inf)u"s",
              γ ::TypeND(Γ0D, [0,1])=γ¹H,
              dt::T0D=(4e-6)u"s")

  nM, nT = maximum([size(x,1) for x in (T1,T2,γ)]), size(B,3)

  # in unit, convert relaxations into losses/recovery per step
  E1 = isa(dt/T1[1], Real) ? exp.(-dt./T1) : uconvert.(NoUnits,exp.(-dt./T1))
  E2 = isa(dt/T2[1], Real) ? exp.(-dt./T2) : uconvert.(NoUnits,exp.(-dt./T2))
  U, Φ = B2UΦ(B; γ=γ, dt=dt)
  return
end

#= blochSim =#
export blochSim
"""
    blochSim(M, B; T1, T2, γ, dt, doHist)
Old school 𝐵-effective magnetic field, `B`, based bloch simulation. Globally or
spin-wisely apply `B` over spins, `M`.

*INPUTS*:
- `M::TypeND(Real, [2])` (nM, xyz): input spins' magnetizations.
- `B::Union{Base.Generator, TypeND(B0D, [2,3])}`:
  Global, (nT,xyz); Spin-wise, (nM,xyz,nT).
*KEYWORDS*:
- `T1 & T2 ::TypeND(T0D, [0,1])`: Global, (1,); Spin-wise, (nM,1).
- `γ::TypeND(Γ0D, [0,1])`: Global, (1,); Spin-wise, (nM, 1). gyro ratio
- `dt::T0D` (1,), simulation temporal step size, i.e., dwell time.
- `doHist::Bool`, whether to output spin history through out `B`.
*OUTPUTS*:
- `M::TypeND(Real, [2])` (nM, xyz): spins after applying `B`.
- `Mhst::TypeND(Real, [3])` (nM, xyz, nT): spins history during `B`.

# Notes:
1. Not much sanity check inside this function, user is responsible for
   matching up the dimensions.
2. Put relax at the end of each time step may still be inaccurate, since
   physically spins relax continuously, this noise/nuance may worth study
   for applications like fingerprinting simulations, etc.
"""
function blochSim(M ::TypeND(AbstractFloat, [2]), B::TypeND(B0D, [2,3]);
                  T1::TypeND(T0D, [0,1])=(Inf)u"s",
                  T2::TypeND(T0D, [0,1])=(Inf)u"s",
                  γ ::TypeND(Γ0D, [0,1])=γ¹H,
                  dt::T0D=(4e-6)u"s", doHist=false)

  if size(B,3) == 1 && size(B,1) != size(M,1)
    B = permutedims(B[:,:,:], [3,2,1])  # best practice?
    println("B not being spin-specific, assuming global")
  end

  # in unit, convert relaxations into losses/recovery per step
  E1 = isa(dt/T1[1], Real) ? exp.(-dt./T1) : uconvert.(NoUnits,exp.(-dt./T1))
  E2 = isa(dt/T2[1], Real) ? exp.(-dt./T2) : uconvert.(NoUnits,exp.(-dt./T2))
  E1₋₁ = E1 .- 1

  U, Φ = B2UΦ(B; γ=γ, dt=dt)

  Mhst = doHist ? zeros(size(M,1), 3, size(U,3)) : nothing

  if doHist
    for t = axes(U,3)
      @inbounds(M = UΦRot(U[:,:,t], Φ[:,1,t], M))
      # relaxation
      M[:,1:2] .*= E2
      M[:,3]   .*= E1
      M[:,3]   .-= E1₋₁
      @inbounds(Mhst[:,:,t] = M)
    end
  else
    for t = axes(U,3)
      @inbounds(M = UΦRot(U[:,:,t], Φ[:,1,t], M))
      # relaxation
      M[:,1:2] .*= E2
      M[:,3]   .*= E1
      M[:,3]   .-= E1₋₁
    end
  end

  return M, Mhst
end

function blochSim(M ::TypeND(AbstractFloat, [2]), B::Base.Generator;
                  T1::TypeND(T0D, [0,1])=(Inf)u"s",
                  T2::TypeND(T0D, [0,1])=(Inf)u"s",
                  γ ::TypeND(Γ0D, [0,1])=γ¹H,
                  dt::T0D=(4e-6)u"s", doHist=false)

  # in unit, convert relaxations into losses/recovery per step
  E1 = isa(dt/T1[1], Real) ? exp.(-dt./T1) : uconvert.(NoUnits,exp.(-dt./T1))
  E2 = isa(dt/T2[1], Real) ? exp.(-dt./T2) : uconvert.(NoUnits,exp.(-dt./T2))
  E1₋₁ = E1 .- 1

  Mhst = doHist ? zeros(size(M,1), 3, length(U)) : nothing

  if doHist
    for (t, b) in enumerate(B)
      u, ϕ = B2UΦ(b; γ=γ, dt=dt)
      @inbounds(M = UΦRot(u, view(ϕ,:,1), M))
      # relaxation
      M[:,1:2] .*= E2
      M[:,3]   .*= E1
      M[:,3]   .-= E1₋₁
      @inbounds(Mhst[:,:,t] = M)
    end
  else
    for (t, b) in enumerate(B)
      u, ϕ = B2UΦ(b; γ=γ, dt=dt)
      @inbounds(M = UΦRot(u, view(ϕ,:,1), M))
      # relaxation
      M[:,1:2] .*= E2
      M[:,3]   .*= E1
      M[:,3]   .-= E1₋₁
    end
  end

  return M, Mhst
end

blochSim(M::TypeND(Integer,[2]), a...; kw...) = blochSim(float(M), a...; kw...)

## Interfaces for mObjects
### AbstractSpinArray
"""
    blochSim(spa::AbstractSpinArray, B; dt, doHist)
Apply old school 𝐵-effective based Bloch simulation on `spa::AbstractSpinArray`,
which brings its own `M, T1, T2, γ`.
"""
blochSim(spa::AbstractSpinArray, B; kw...) =
  blochSim(spa.M, B; T1=spa.T1, T2=spa.T2, γ=spa.γ, kw...)

"""
    blochSim(spa::AbstractSpinArray, p::Pulse, loc; Δf, b1Map, doHist)
Similar to `blochSim(spa::AbstractSpinArray, B; dt, doHist)`, except for that
`B = Pulse2B(p, spa, loc; Δf, b1Map)`, and that `dt = p.dt`.
"""
blochSim(spa::AbstractSpinArray, p::Pulse, loc; doHist=false, kw...) =
  blochSim(spa, Pulse2B(p, spa, loc; kw...); dt=p.dt, doHist=doHist)

### AbstractSpinCube
"""
    blochSim(cb::AbstractSpinCube, B; dt, doHist)
Apply old school 𝐵-effective based Bloch simulation on `cb::AbstractSpinCube`,
which brings its own `M, T1, T2, γ`.
"""
blochSim(cb::AbstractSpinCube,B; kw...) =blochSim(cb.spinarray,B; kw...)

"""
    blochSim(cb::AbstractSpinCube, p::Pulse; b1Map, doHist)
Similar to `blochSim(cb::AbstractSpinCube, B; dt, doHist)`, except for that
`B = Pulse2B(p, cb; b1Map)`, and that `dt = p.dt`.
"""
blochSim(cb::AbstractSpinCube, p::Pulse; doHist=false, kw...) =
  blochSim(cb, Pulse2B(p, cb.spa; kw...); dt=p.dt, doHist=doHist)

