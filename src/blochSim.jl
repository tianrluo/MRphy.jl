# include("mObjects.jl")

#= Pulse2B =#
export Pulse2B
"""
    B = Pulse2B(rf, gr, loc; Δf, b1Map, γ)
Turn RF, `rf`, and gradient, `gr`, into 𝐵-effective magnetic field.

*INPUTS*:
- `rf::TypeND(RF0D, [1,2])` (nSteps, (nCoil))
- `gr::TypeND(GR0D, [2])`   (nSteps, 3)
- `loc::TypeND(L0D, [2,3])` (1,3) or (nSpins, 3, (nSteps)), locations.
*OPTIONALS*:
- `Δf::TypeND(F0D, [0,1,2,3])` (1,)  or (nSpins, (1, (nSteps))), off-resonance.
- `b1Map::TypeND(Union{Real, Complex}, [0,1,2,3])` (1,) or
  (nSpins, (nCoils, (nSteps))), transmit sensitivity
- `γ::TypeND(Γ0D, [0,1])` (1,)  or (nSpins,), gyro-ratio
*OUTPUS*:
- `B::TypeND(B0D, [3])` (1,3,nSteps) or (nSpins, 3, nSteps), 𝐵 field.
"""
function Pulse2B(rf   ::TypeND(RF0D, [1,2]),
                 gr   ::TypeND(GR0D, [2]),
                 loc  ::TypeND(L0D,  [2,3]) = [0 0 0]u"cm";
                 Δf   ::TypeND(F0D,  [0,1,2,3]) = 0u"Hz",
                 b1Map::TypeND(Union{Real,Complex}, [0,1,2,3]) = 0,
                 γ    ::TypeND(Γ0D,  [0,1]) = γ¹H)
  rf, gr = map(x->permutedims(x[:,:,:], [3,2,1]), (rf, gr))
  Bz  = (Δf==0u"Hz") ? sum(loc.*gr, dims=2) : sum(loc.*gr, dims=2).+(Δf./γ)
  if b1Map!=0 rf = sum(b1Map.*rf, dims=2) end
  return [repeat([real(rf) imag(rf)], size(Bz,1),1,1) Bz]
end

Pulse2B(::NoUnitChk,
        rf   ::TypeND(Union{Real, Complex}, [1,2]),
        gr   ::TypeND(Real, [2]),
        loc  ::TypeND(Real, [2,3]) = [0 0 0];
        Δf   ::TypeND(Real, [0,1,2,3]) = 0,
        b1Map::TypeND(Union{Real,Complex},[0,1,2,3]) = 0,
        γ    ::TypeND(Real, [0,1]) = ustrip(γ¹H)) =
  Pulse2B(Quantity.(rf, u"Gauss"), Quantity.(gr, u"Gauss/cm"),
          Quantity.(loc, u"cm");   Δf=Quantity.(Δf, u"Hz"),
          b1Map=b1Map, γ=Quantity.(γ, u"Hz/Gauss"))

"""
    B = Pulse2B(pulse::Pulse, loc; Δf, b1Map, γ)
Turn struct `Pulse` into effective magnetic, 𝐵, field.
"""
Pulse2B(p::Pulse, loc; kargs...) = Pulse2B(p.rf, p.gr, loc; kargs...)
Pulse2B(::NoUnitChk, p::Pulse, loc; kargs...) =
  Pulse2B(NoUnitChk(), p, loc; kargs...)

"""
    B = Pulse2B(pulse::Pulse, spa::AbstractSpinArray, loc; Δf, b1Map)
...with `γ=spa.γ`.
"""
Pulse2B(p::Pulse, spa::AbstractSpinArray, loc; kargs...) =
  Pulse2B(p, loc; γ=spa.γ, kargs...)
Pulse2B(::NoUnitChk, p::Pulse, spa::AbstractSpinArray, loc; kargs...) =
  Pulse2B(NoUnitChk(), p, loc; γ=spa.γ, kargs...)

"""
    B = Pulse2B(pulse::Pulse, cb::AbstractSpinCube; b1Map)
...with `loc, Δf, γ = cb.loc, cb.Δf, cb.γ`.
"""
Pulse2B(p::Pulse, cb::AbstractSpinCube; kargs...) =
  Pulse2B(p, cb.loc, cb.Δf; γ=cb.γ, kargs...)

#= blochSim, old school =#
export blochSim
"""
    blochSim(M, B; T1, T2, γ, dt, doHist)
Old school 𝐵-effective magnetic field, `B`, based bloch simulation. Globally or
spin-wisely apply `B` over spins, `M`.

*INPUTS*:
- `M::TypeND(Real, [2])` (nSpins, xyz): input spins' magnetizations.
- `B::TypeND(B0D, [2,3])`: Global, (nSteps,xyz); Spin-wise, (nSpins,xyz,nSteps).
*OPTIONALS*:
- `T1 & T2 ::TypeND(T0D, [0,1])`: Global, (1,); Spin-wise, (nSpins,1).
- `γ::TypeND(Γ0D, [0,1])`: Global, (1,); Spin-wise, (nSpins, 1). gyro ratio
- `dt::T0D` (1,), simulation temporal step size, i.e., dwell time.
- `doHist::Bool`, whether to output spin history through out `B`.
*OUTPUTS*:
- `M::TypeND(Real, [2])` (nSpins, xyz): spins after applying `B`.
- `Mhst::TypeND(Real, [3])` (nSpins, xyz, nSteps): spins history during `B`.

# Notes:
1. Not much sanity check inside this function, user is responsible for
   matching up the dimensions.
2. Put relax at the end of each time step may still be inaccurate, since
   physically spins relax continuously, this noise/nuance may worth study
   for applications like fingerprinting simulations, etc.
"""
function blochSim(::NoUnitChk,
                  M ::TypeND(Real, [2]),
                  B ::TypeND(Real, [2,3]);
                  T1::TypeND(Real, [0,1])=Inf,
                  T2::TypeND(Real, [0,1])=Inf,
                  γ ::TypeND(Real, [0,1])=4257.6,
                  dt::Real=4e-6,
                  doHist::Bool=false)

  if size(B,3) == 1 && size(B,1) != size(M,1)
    B = permutedims(B[:,:,:], [3,2,1])  # best practice?
    println("B not being spin-specific, assuming global")
  end
  nSpins, nSteps = size(M,1), size(B,3)

  # in unit, convert relaxations into losses/recovery per step
  E1, E2 = exp.(-dt./T1), exp.(-dt./T2)

  Bmag = sqrt.(sum(B.*B, dims=2))
  niBmag = -1 ./ Bmag  # negative to correct cross product (x-prod) direction;
  niBmag[isinf.(niBmag)] .= 0 # trick, nan-check costs 3-times of inf-check
  Bn = B.*niBmag  # Bn, normalized B

  theta = Bmag.*γ.*2π.*dt
  Ct, St = cos.(theta), sin.(theta)

  StBn, CtBn_1 = Bn.*St, Bn.*(1 .- Ct)

  Mx0, My0, Mz0 = M[:,1], M[:,2], M[:,3]

  if doHist Mox, Moy, Moz = [zeros(nSpins,nSteps) for i = 1:3] end

  doFlag = any(Bmag .!= 0, dims=1)

  # Updates are calculated with rotation-then-decay, instead of the canonical
  # differential equation expression.
  # Discretizing differential equation may cause precision issue.
  # full lower-case variables are local in loop, o.w. not local
  for istep = 1:nSteps
    if doFlag[istep]
      # step-wisely extract pre-processed variables
      bn, stbn, ctbn_1 = Bn[:,:,istep], StBn[:,:,istep], CtBn_1[:,:,istep]
      ct = Ct[:, :, istep]

      ip = sum(bn .* [Mx0 My0 Mz0], dims=2)  # vector inner product

      # explicitly express cross(bn, Mo_ii_1) as a matrix vector multiplication
      mx1 =  ct       .*Mx0 .-stbn[:,3].*My0 .+stbn[:,2].*Mz0 .+ip.*ctbn_1[:,1]
      my1 =  stbn[:,3].*Mx0 .+ct       .*My0 .-stbn[:,1].*Mz0 .+ip.*ctbn_1[:,2]
      mz1 = -stbn[:,2].*Mx0 .+stbn[:,1].*My0 .+ct       .*Mz0 .+ip.*ctbn_1[:,3]
    else
      mx1, my1, mz1 = Mx0, My0, Mz0
    end
    # relaxation effects: "1" in Mz0 since M0=1 by assumption
    # also, update Mo_ii_1.
    Mx0, My0, Mz0 = mx1.*E2, my1.*E2, mz1.*E1 .+ 1 .- E1

    if doHist Mox[:,istep], Moy[:,istep], Moz[:,istep] = Mx0, My0, Mz0 end
  end

  M = [Mx0 My0 Mz0]
  if doHist Mhst = permutedims(cat(dims=3, Mox, Moy, Moz), [1,3,2])
  else      Mhst = nothing
  end

  return M, Mhst
end

blochSim(M ::TypeND(Real, [2]),
         B ::TypeND(B0D, [2,3]);
         T1::TypeND(T0D, [0,1])=(Inf)u"s",
         T2::TypeND(T0D, [0,1])=(Inf)u"s",
         γ ::TypeND(Γ0D, [0,1])=(4257.6)u"Hz/Gauss",
         dt::T0D=(4e-6)u"s", doHist=false) =
  blochSim(NoUnitChk(), M, ustrip.(u"Gauss", B);
           T1=ustrip.(u"s",T1), T2=ustrip.(u"s",T2), γ=ustrip.(u"Hz/Gauss",γ),
           dt=ustrip.(u"s",dt), doHist=doHist)

## Interfaces for mObjects
### AbstractSpinArray
"""
    blochSim(spa::AbstractSpinArray, B; dt, doHist)
Apply old school 𝐵-effective based Bloch simulation on `spa::AbstractSpinArray`,
which brings its own `M, T1, T2, γ`.
"""
blochSim(::NoUnitChk, spa::AbstractSpinArray, B; dt=4e-6, doHist=false) =
  blochSim(NoUnitChk(), spa.M, B;
           T1=ustrip.(u"s", spa.T1), T2=ustrip.(u"s",spa.T2),
           γ=ustrip.(u"Hz/Gauss", spa.γ), dt=dt, doHist=doHist)

blochSim(spa::AbstractSpinArray, B; dt=(4e-6)u"s", doHist=false) =
  blochSim(NoUnitChk(), spa, ustrip.(u"Gauss", B);
           dt=ustrip.(u"s",dt), doHist=doHist)

"""
    blochSim(spa::AbstractSpinArray, p::Pulse, loc; Δf, b1Map, doHist)
Similar to `blochSim(spa::AbstractSpinArray, B; dt, doHist)`, except for that
`B = Pulse2B(p, spa, loc; Δf, b1Map)`, and that `dt = p.dt`.
"""
blochSim(spa::AbstractSpinArray, p::Pulse, loc; doHist=false, kargs...) =
  blochSim(spa, Pulse2B(p, spa, loc; kargs...); dt=p.dt, doHist=doHist)

### AbstractSpinCube
"""
    blochSim(cb::AbstractSpinCube, B; dt, doHist)
Apply old school 𝐵-effective based Bloch simulation on `cb::AbstractSpinCube`,
which brings its own `M, T1, T2, γ`.
"""
blochSim(::NoUnitChk, cb::AbstractSpinCube, B; dt=4e-6, doHist=false) =
  blochSim(NoUnitChk(), cb.spinarray, B; dt=dt, doHist=doHist)

blochSim(cb::AbstractSpinCube, B; dt=(4e-6)u"s", doHist=false) =
  blochSim(cb.spinarray, B; dt=dt, doHist=doHist)

"""
    blochSim(cb::AbstractSpinCube, p::Pulse; b1Map, dt, doHist)
Similar to `blochSim(cb::AbstractSpinCube, B; dt, doHist)`, except for that
`B = Pulse2B(p, cb; b1Map)`, and that `dt = p.dt`.
"""
blochSim(cb::AbstractSpinCube, p::Pulse; doHist=false, kargs...) =
  blochSim(cb, Pulse2B(cb, spa; kargs...); dt=p.dt, doHist=doHist)

