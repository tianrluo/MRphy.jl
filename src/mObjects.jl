"""Throw `ArgumentError` when `\$x` is an immutable field."""
ExceptionImmutableField(x) = ArgumentError("`$x` is an immutable field.")

#= Pulse =#

export AbstractPulse
"""
An abstract type for pulses.

See also: [`Pulse`](@ref).
"""
abstract type AbstractPulse end

## Basics
Base.isequal(a::AbstractPulse, b::AbstractPulse) =
  all([isequal(getproperty(a,s), getproperty(b,s)) for s in fieldnames(Pulse)])

export Pulse
"""
A struct for typical MR pulses: `Pulse <: AbstractPulse`.

# Fields:
*Mutable*:
- `rf::TypeND(RF0D, [1,2])` (nT,) or (nT, nCoils).
- `gr::TypeND(GR0D, [2])` (nT, 3), where 3 accounts for x-y-z channels.
- `dt::T0D` (1,), simulation temporal step size, i.e., dwell time.
- `des::String`, an description of the pulse to be constructed.

# Usages:
    Pulse(rf, gr; dt=(4e-6)u"s", des="generic pulse")

See also: [`AbstractPulse`](@ref).
"""
mutable struct Pulse <: AbstractPulse
  rf::TypeND(RF0D, [1,2])
  gr::TypeND(GR0D, [2])
  dt::T0D
  des::String

  function Pulse(rf=missing, gr=missing; dt=4e-6u"s", des="generic pulse")
    rf_miss, gr_miss = ismissing(rf), ismissing(gr)
    rf_miss&&gr_miss && ErrorException("Missing both inputs.")
    if rf_miss rf = zeros(size(gr,1))u"Gauss" end
    if gr_miss gr = zeros(size(rf,1),3)u"Gauss/cm" end
    size(gr,2)==3 || throw(DimensionMismatch)

    if isa(rf, Number) rf = [rf] end
    return new(rf, gr, dt, des)
  end

end

## set and get
Base.setproperty!(p::Pulse, sym::Symbol, x) = begin
  if (sym==:gr) size(x)==(size(p.rf,1),3)||throw(DimensionMismatch) end
  if (sym==:rf) size(x,1)==size(p.gr,1)||throw(DimensionMismatch) end
  setfield!(p, sym, x)
end

#= Spin =#

export AbstractSpinArray
"""
This type keeps the essentials of magnetic spins. Its instance struct must
contain all fields listed listed in the exemplary struct `mSpinArray`.

# Misc
Might make `AbstractSpinArray <: AbstractArray` in a future version

See also: [`mSpinArray`](@ref), [`AbstractSpinCube`](@ref).
"""
abstract type AbstractSpinArray end

## set and get
Base.setproperty!(spa::AbstractSpinArray, s::Symbol, x) = begin
  s ∈ (:dim, :mask) && throw(ExceptionImmutableField(s))

  nM = prod(spa.dim)
  if (s==:M)&&(size(x,1)==1) x=repeat(x, nM, 1) end
  if (s ∈ (:T1,:T2,:γ,:M)) size(x,1)∈(1,nM)||throw(DimensionMismatch) end

  setfield!(spa, s, x)
end

## AbstractArray-like interface
Base.size(spa::AbstractSpinArray) = spa.dim
Base.size(spa::AbstractSpinArray, d) = (d ≤ length(spa.dim)) ? spa.dim[d] : 1

Base.isequal(a::AbstractSpinArray, b::AbstractSpinArray) =
  all([isequal(getproperty.((a,b),s)...) for s in fieldnames(mSpinArray)])

## Concrete mSpinArray
export mSpinArray
"""
    mSpinArray(dim::Dims; T1=1.47u"s", T2=0.07u"s", γ=γ¹H, M=[0. 0. 1.])
    mSpinArray(mask::BitArray; T1=1.47u"s", T2=0.07u"s", γ=γ¹H, M=[0. 0. 1.])
An exemplary struct instantiating `AbstractSpinArray`.

# Fields:
*Immutable*:
- `dim::Dims` (nd,): `nM ← prod(dim)`, dimension of the object.
- `mask::BitArray` (nx,(ny,(nz))): Mask for `M`, `dim == (nx,(ny,(nz)))`
*Mutable*:
- `T1::TypeND(T0D, [0,1])` (1,) or (nM,): Longitudinal relaxation coeff.
- `T2::TypeND(T0D, [0,1])` (1,) or (nM,): Transversal relaxation coeff.
- `γ::TypeND(Γ0D, [0,1])`  (1,) or (nM,): Gyromagnetic ratio.
- `M::TypeND(Real, [2])`   (`count(mask)`, 3):  Magnetic spins, (𝑀x,𝑀y,𝑀z).

# Notes:
off-resonance, `Δf`, and locations, `loc`, are intentionally unincluded, as they
are not intrinsic to spins, and can change over time. Unincluding them allows
extensional subtypes specialized for, e.g., arterial spin labelling.

See also: [`AbstractSpinArray`](@ref).
"""
mutable struct mSpinArray <: AbstractSpinArray
  # *Immutable*:
  dim ::Dims
  mask::BitArray
  # *Mutable*:
  T1 ::TypeND(T0D, [0,1])
  T2 ::TypeND(T0D, [0,1])
  γ  ::TypeND(Γ0D, [0,1])
  M  ::TypeND(AbstractFloat, [2])

  function mSpinArray(mask::BitArray;
                      T1=1.47u"s", T2=0.07u"s", γ=γ¹H, M=[0. 0. 1.])
    dim = size(mask)
    nM = prod(dim)
    M = eltype(M)<:AbstractFloat ? M : float(M)
    if size(M,1)==1 M=repeat(M, nM, 1) end
    (all(map(x->(size(x,1) ∈ (1,nM)), [T1,T2,γ,M]))) || throw(DimensionMismatch)

    return new(dim, mask, T1, T2, γ, M)
  end

  mSpinArray(dim::Dims=(1,); kw...) = mSpinArray(trues(dim); kw...)

end

#= Cube =#

export AbstractSpinCube
"""
`AbstractSpinCube <: AbstractSpinArray`.
This type inherits `AbstractSpinArray` as a field. Its instance struct must
contain all fields listed in the exemplary struct `mSpinCube`.

See also: [`AbstractSpinArray`](@ref), [`mSpinCube`](@ref).
"""
abstract type AbstractSpinCube <: AbstractSpinArray end

## set and get
Base.setproperty!(cb::AbstractSpinCube, s::Symbol, x) = begin
  s ∈ (:spinarray,:fov,:ofst,:loc) && throw(ExceptionImmutableField(s))
  s ∈ fieldnames(typeof(cb)) ? setfield!(cb, s,x) : setfield!(cb.spinarray, s,x)
end

Base.getproperty(cb::AbstractSpinCube, s::Symbol) =
  s ∈ fieldnames(typeof(cb)) ? getfield(cb, s) : getfield(cb.spinarray, s)

## AbstractArray-like interface
Base.size(cb::AbstractSpinCube, a...) = Base.size(cb.spinarray, a...)

Base.isequal(a::AbstractSpinCube, b::AbstractSpinCube) =
  all([isequal(getproperty.((a,b),s)...) for s in fieldnames(mSpinCube)])

## Concrete mSpinCube
export mSpinCube
"""
    spincube = mSpinCube(dim::Dims{3}, fov; ofst, Δf, T1, T2, γ)
    spincube = mSpinCube(mask::BitArray{3}, fov; ofst, Δf, T1, T2, γ)
`dim`, `mask`, `T1`, `T2`, and `γ` are passed to `mSpinArray` constructors.

An exemplary struct instantiating `AbstractSpinCube`, designed to model a set of
regularly spaced spins, e.g., a volume.

# Fields:
*Immutable*:
- `spinarray::AbstractSpinArray` (1,): inherited `AbstractSpinArray` struct
- `fov ::TypeND(L0D, [2])` (1, 3): field of view.
- `ofst::TypeND(L0D, [2])` (1, 3): fov offset from magnetic field iso-center.
- `loc ::TypeND(L0D, [2])` (nM, 3): location of spins.
*Mutable*:
- `Δf::TypeND(F0D, [0,1])` (1,) or (nM,): off-resonance map.

See also: [`AbstractSpinCube`](@ref).
"""
mutable struct mSpinCube <: AbstractSpinCube
  # *Immutable*:
  spinarray::AbstractSpinArray
  fov ::TypeND(L0D, [2])
  ofst::TypeND(L0D, [2])
  loc ::TypeND(L0D, [2])
  # *Mutable*:
  Δf  ::TypeND(F0D, [0,1])

  function mSpinCube(mask::BitArray{3}, fov::TypeND(L0D, [2]);
                     ofst::TypeND(L0D, [2])=[0 0 0]u"cm", Δf=0u"Hz",
                     T1=1.47u"s", T2=0.07u"s", γ=γ¹H)
    size(fov)==size(ofst)==(1,3) || throw(DimensionMismatch)
    spa = mSpinArray(mask; T1=T1, T2=T2, γ=γ)
    loc = CartesianLocations(spa.dim)./(reshape([spa.dim...], 1,:)./fov) .+ ofst
    return new(spa, fov, ofst, loc, Δf)
  end

  mSpinCube(dim::Dims{3}, a...; kw...) = mSpinCube(trues(dim), a...; kw...)

end

#= Bolus (*Under Construction*) =#
# TODO
# export AbstractSpinBolus
"""
*UNDER CONSTRUCTION*

`AbstractSpinBolus <: AbstractSpinArray`.
This type inherits `AbstractSpinArray` as a field. Its instance struct must
contain all fields listed in the exemplary struct `mSpinBolus`.

See also: [`AbstractSpinArray`](@ref), [`mSpinBolus`](@ref).
"""
abstract type AbstractSpinBolus <: AbstractSpinArray end

## set and get

Base.getproperty(bl::AbstractSpinBolus, s::Symbol) =
  s ∈ fieldnames(typeof(bl)) ? getfield(bl, s) : getfield(bl.spinarray, s)

## AbstractArray-like interface
Base.size(bl::AbstractSpinBolus, a...) = Base.size(bl.spinarray, a...)

Base.isequal(a::AbstractSpinBolus, b::AbstractSpinBolus) =
  all([isequal(getproperty.((a,b),s)...) for s in fieldnames(mSpinBolus)])

## Concrete mSpinBolus
# export mSpinBolus
"""
*UNDER CONSTRUCTION*

An exemplary struct instantiating `AbstractSpinBolus`, designed to model a set
of moving spins, e.g., a blood bolus in ASL context.

See also: [`AbstractSpinBolus`](@ref).
"""
mutable struct mSpinBolus <: AbstractSpinBolus
  # *Immutable*:
  # *Mutable*:
end

#= mObjects utils =#

export Pulse2B
"""
    B = Pulse2B(pulse::Pulse, loc; Δf, b1Map, γ)
Create effective magnetic field, 𝐵, from input `pulse`.

See also: [`rfgr2B`](@ref), [`B2UΦ`](@ref), [`blochSim`](@ref).
"""
Pulse2B(p::Pulse, loc; kw...) = rfgr2B(p.rf, p.gr, loc; kw...)

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

export applyPulse, applyPulse!
"""
    applyPulse(spa::AbstractSpinArray, p::Pulse, loc; Δf, b1Map, doHist)
Turn `p` into 𝐵-effective and apply it on `spa.M`, using its own `M, T1, T2, γ`.

See also: [`blochSim`](@ref), [`freePrec`](@ref).
"""
applyPulse(spa::AbstractSpinArray, p::Pulse, loc; doHist=false, kw...) =
  blochSim(spa.M, Pulse2B(p, spa, loc; kw...);
           T1=spa.T1, T2=spa.T2, γ=spa.γ, dt=p.dt, doHist=doHist)

"""
    applyPulse!(spa::AbstractSpinArray, p::Pulse, loc; Δf, b1Map, doHist)
Update `spa.M` before return.
"""
applyPulse!(spa::AbstractSpinArray, p::Pulse, loc; doHist=false, kw...) = begin
  _, Mhst = blochSim!(spa.M, Pulse2B(p, spa, loc; kw...);
                      T1=spa.T1, T2=spa.T2, γ=spa.γ, dt=p.dt, doHist=doHist)
  return (M=spa.M, Mhst=Mhst)
end

"""
    applyPulse(cb::AbstractSpinCube, p::Pulse; b1Map, doHist)
Turn `p` into 𝐵-effective and apply it on `cb.M`, using its own `M, T1, T2, γ`.
"""
applyPulse(cb::AbstractSpinCube, p::Pulse; b1Map=1, doHist=false) =
  blochSim(cb.M, Pulse2B(p, cb; b1Map=b1Map);
           T1=cb.T1, T2=cb.T2, γ=cb.γ, dt=p.dt, doHist=doHist)

"""
    applyPulse!(cb::AbstractSpinCube, p::Pulse; b1Map, doHist)
Update `cb.M` before return.
"""
applyPulse!(cb::AbstractSpinCube, p::Pulse; b1Map=1, doHist=false) = begin
  _, Mhst = blochSim!(cb.M, Pulse2B(p, cb; b1Map=b1Map);
                      T1=cb.T1, T2=cb.T2, γ=cb.γ, dt=p.dt, doHist=doHist)
  return (M=cb.M, Mhst=Mhst)
end

export freePrec!, freePrec
"""
    freePrec(spa::AbstractSpinArray, t; Δf)
`spa::AbstractSpinArray` free precess by `t`. `spa.M` will not be updated.

See also: [`applyPulse`](@ref), [`freePrec!`](@ref).
"""
freePrec(spa::AbstractSpinArray, t::T0D; Δf::TypeND(F0D,[0,1])=0u"Hz") =
  freePrec(spa.M, t; Δf=Δf, T1=spa.T1, T2=spa.T2)

"""
    freePrec!(spa::AbstractSpinArray, t; Δf)
...`spa.M` will updated by the results.
"""
freePrec!(spa::AbstractSpinArray, t::T0D; Δf::TypeND(F0D,[0,1])=0u"Hz") =
  freePrec!(spa.M, t; Δf=Δf, T1=spa.T1, T2=spa.T2)

"""
    freePrec(cb::AbstractSpinCube, t)
`cb::AbstractSpinCube` free precess by `t`. `cb.M` will not be updated.

See also: [`applyPulse`](@ref), [`freePrec`](@ref).
"""
freePrec(cb::AbstractSpinCube, t::T0D) =
  freePrec(cb.M, t; Δf=cb.Δf, T1=cb.T1, T2=cb.T2)

"""
    freePrec!(cb::AbstractSpinCube, t)
...`cb.M` will be updated by the results.
"""
freePrec!(cb::AbstractSpinCube, t::T0D) =
  freePrec!(cb.M, t; Δf=cb.Δf, T1=cb.T1, T2=cb.T2)

