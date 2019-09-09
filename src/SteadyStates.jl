export SteadyStates
"""
Some steady state properties of common sequences.
"""
module SteadyStates
using ..MRphy
using ..Unitful, ..UnitfulMR

export Signal
"""
Analytical expressions of common steady states sequences signals.
"""
module Signal
using ..MRphy
using ..Unitful, ..UnitfulMR

"""
    bSSFP(α; TR, Δf, T1, T2)
`TE=0` Steady state bSSFP signal. 10.1002/jmri.24163, eq.(4), with `ϕ=2π*Δf*TR`.

*INPUTS*:
- `α::Real` (1,), tip angle in degree;
*KEYWORDS*:
- `TR::T0D` (1,), repetition time;
- `Δf::F0D` (1,), off-resonance in Hz;
- `T1::T0D` (1,), longitudinal relaxation coefficient;
- `T2::T0D` (1,), transverse relaxation coefficient;
*OUTPUTS*:
- `sig::Complex` (1,), steady-state signal.

See also: [`SPGR`](@ref), [`STFR`](@ref).
"""
function bSSFP(α::Real;
               TR::T0D,
               Δf::F0D=0u"Hz",
               T1::T0D=1u"s",
               T2::T0D=(6e-2)u"s")
  E1, E2, ϕ = exp(-TR/T1), exp(-TR/T2), 2π*Δf*TR # jmri.24163 typoed on E₁, E₂.
  C, D = E2*(E1-1)*(1+cosd(α)), (1-E1*cosd(α)) - (E1-cosd(α))*E2^2
  sig = sind(α) * (1-E1)*(1-E2*exp(-1im*ϕ)) / (C*cos(ϕ) + D)
  return sig
end

"""
    SPGR(α; TR, T1)
`TE=0` Steady state SPGR signal. 10.1002/mrm.1910130109, eq.(1), ideal spoiling.

*INPUTS*:
- `α::Real` (1,), tip angle in degree;
*KEYWORDS*:
- `TR::T0D` (1,), repetition time;
- `T1::T0D` (1,), longitudinal relaxation coefficient;
*OUTPUTS*:
- `sig::Real` (1,), steady-state signal.

See also: [`bSSFP`](@ref), [`STFR`](@ref).
"""
SPGR(α::Real; TR::T0D, T1::T0D=1u"s") =
  sind(α) * (1-exp(-TR/T1))/(1-cosd(α)*exp(-TR/T1));

"""
    STFR(α, β; ϕ, Δf, T1, T2, Tg, Tf)
`TE=0` Steady state STFR signal. 10.1002/mrm.25146, eq.(2), ideal spoiling.

*INPUTS*:
- `α::Real` (1,), tip-down angle in degree;
- `β::Real` (1,), tip-up angle in degree;
*KEYWORDS*:
- `ϕ::Real` (1,), phase of the tip-up pulse in radians;
- `Δf::F0D` (1,), off-resonance in Hz;
- `T1::T0D` (1,), longitudinal relaxation coefficient;
- `T2::T0D` (1,), transverse relaxation coefficient;
- `Tg::T0D` (1,), duration of gradient crusher;
- `Tf::T0D` (1,), duration of free precession in each TR;
*OUTPUTS*:
- `sig::Real` (1,), steady-state signal.

See also: [`bSSFP`](@ref), [`SPGR`](@ref).
"""
function STFR(α::Real, β::Real;
              ϕ::Real=0,
              Δf::F0D=0u"Hz",
              T1::T0D=1u"s",
              T2::T0D=(6e-2)u"s",
              Tg::T0D=(2e-3)u"s",
              Tf::T0D=(1e-2)u"s")
  cΔϕ = cos(2π*Δf*Tf - ϕ) # θf ≔ 2π*Δf*Tf
  Eg, Ef1, Ef2 = exp(-Tg/T1), exp(-Tf/T1), exp(-Tf/T2)
  sα, cα, sβ, cβ = sind(α), cosd(α), sind(β), cosd(β)

  sig = sα * (Eg*(1-Ef1)*cβ + (1-Eg)) / (1-Eg*Ef2*sα*sβ*cΔϕ-Eg*Ef1*cα*cβ)
  return sig
end

end # module Signal

export RFSpoiling
"""
Tools for simulating RF spoiling in gradient echo sequences.
"""
module RFSpoiling
using ..MRphy
using ..Unitful, ..UnitfulMR
"""
    QuadPhase(nTR::Integer, C::Real, B::Real, A::Real)
Quadratically cycling phases in (°): Φ(n) = mod.(C⋅n² + B⋅n + A, 360), n=0:nTR-1
"""
QuadPhase(nTR::Integer, C::Real, B::Real=0, A::Real=0) =
  mod.(([0:nTR-1 ...]|>n->C*n.^2 .+B*n .+A), 360)

"""
    FZstates(Φ, α; TR, T1, T2, FZ)
𝐹, 𝑍 from: 10.1002/(SICI)1099-0534(1999)11:5<291::AID-CMR2>3.0.CO;2-J, eq.(7,8).
eq.(7) refined to `÷√(2)`, instead of `÷2`, as in 10.1002/mrm20736: eq.(2).

Assume constant gradient spoiling of m⋅2π dephasing in each TR, m∈ℤ.
In practice, if dephase is a constant but not exactly m⋅2π, the resulting states
can be computed by convolving a sinc with the m⋅2π dephased results.

*INPUTS*:
- `Φ::TypeND(Real,[2])` (nC,nTR), nC: #`C` as `C` in `QuadPhase`. Typically, one
  simulates a range of `C`s, picking a `C` yielding a signal intensity equals to
  that of ideal spgr spoiling.
- `α::Real` (1,), flip-angle.
*KEYWORDS*:
- `TR::T0D` (1,), repetition time;
- `T1::T0D` (1,), longitudinal relaxation coefficient;
- `T2::T0D` (1,), transverse relaxation coefficient;
- `FZ::NamedTuple`, `(Fs,Fcs,Zs)`, simulate from prescribed states if given:\\
  `Fs ::TypeND(Complex,[2])`, transversal dephasing states, 𝐹ₙ;\\
  `Fcs::TypeND(Complex,[2])`, conjugate transversal dephasing states, 𝐹₋ₙ*;\\
  `Zs ::TypeND(Complex,[2])`, longitudinal states, 𝑍ₙ;
*OUTPUTS*:
- `FZ::NamedTuple`, `(Fs,Fcs,Zs)`, simulation results.
"""
function FZstates(Φ ::TypeND(Real,[2]), α::Real;
                  TR::T0D=(50e-3)u"s",
                  T1::T0D=(1.470)u"s",
                  T2::T0D=(71e-3)u"s",
                  FZ::Union{Nothing, NamedTuple}=nothing)
  if isnothing(FZ)
    (Fs, Fcs, Zs) = (zeros(Complex{Float64}, size(Φ)) for _ in 1:3)
    Zs[:, 1] .= 1
    FZ = (Fs=Fs, Fcs=Fcs, Zs=Zs)
  end
  eltype(Φ) <: AbstractFloat || (Φ = float(Φ))
  eltype(Fs)==Complex{Float64}  || (Fs  = Complex{Float64}(Fs))
  eltype(Fcs)==Complex{Float64} || (Fcs = Complex{Float64}(Fcs))
  eltype(Zs)==Complex{Float64}  || (Zs  = Complex{Float64}(Zs))

  size(FZ.Fs)==size(FZ.Fcs)==size(FZ.Zs) || throw(DimensionMismatch)

  E1, E2 = exp(-TR/T1), exp(-TR/T2)

  return FZ = _FZstates(Φ, α; E1=E1, E2=E2, FZ...)
end

function _FZstates(Φ ::TypeND(AbstractFloat,[2]), α::Real;
                   E1::Real,
                   E2::Real,
                   Fs ::TypeND(Complex,[2]),
                   Fcs::TypeND(Complex,[2]),
                   Zs ::TypeND(Complex,[2]))
  Φ .= mod.(Φ, 360)
  Φ .*= π/180

  Fs⁻, Fcs⁻, Zs⁻ = similar(Fs), similar(Fcs), similar(Zs) # pre-rf states
  Fcs⁻[:,end] .= 0

  c½α², s½α², cα, sα = cosd(α/2)^2, sind(α/2)^2, cosd(α), sind(α)
  ieΦsαr2, e2Φs½α² = im*exp.(im*Φ)*sα/√(2), exp.(im*2*Φ)*s½α²

  @inbounds for t in axes(Φ,2)
    # relaxation and updating pre-rf states; eq.(8)
    Fs⁻ .= E2.*[conj(Fcs[:,2]) Fs[:,1:end-1]]
    Fcs⁻[:, 1:end-1] .= E2.*Fcs[:,2:end]
    Zs⁻ .= E1.*Zs
    Zs⁻[:,1] .+= (1-E1)

    # excitation; eq.(7)
    ieϕsαr2, e2ϕs½α² = view(ieΦsαr2,:,t), view(e2Φs½α²,:,t)

    Fs  .=          c½α².*Fs⁻ .+ e2ϕs½α².*Fcs⁻ .-       ieϕsαr2.*Zs⁻
    Fcs .= conj(e2ϕs½α²).*Fs⁻ .+    c½α².*Fcs⁻ .- conj(ieϕsαr2).*Zs⁻
    Zs  .= conj(ieϕsαr2).*Fs⁻ .+ ieϕsαr2.*Fcs⁻ .+            cα.*Zs⁻
  end

  return (Fs=Fs, Fcs=Fcs, Zs=Zs)
end

end # module RFSpoiling

end # module SteadyStates

