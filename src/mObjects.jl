# include("utils.jl")

error_ImmutableField(x) = ErrorException("`$x` is an immutable field.")

#= Pulse =#

export AbstractPulse
"""
    AbstractPulse
An abstract type for pulses.
"""
abstract type AbstractPulse end

export Pulse
"""
    Pulse <: AbstractPulse
A struct for typical MR pulses.

# Fields:
*Mutable*:
- `rf::TypeND(RF0D, [1,2])` (nT,) or (nT, nCoils).
- `gr::TypeND(GR0D, [2])` (nT, 3), where 3 accounts for x-y-z channels.
- `dt::T0D` (1,), simulation temporal step size, i.e., dwell time.
- `des::String`, an description of the pulse to be constructed.

# Usages:
    pulse = Pulse(rf, gr; dt=dt, des=des)
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
  if sym==:gr @assert(size(x,1) == size(p.rf,1))&&(size(x,2)==3) end
  if sym==:rf @assert(size(x,1) == size(p.gr,1)) end
  setfield!(p, sym, x)
end

## Other Basics
Base.isequal(a::AbstractPulse, b::AbstractPulse) =
  all([isequal(getproperty(a,s), getproperty(b,s)) for s in fieldnames(Pulse)])

#= Spin =#

export AbstractSpinArray
"""
    AbstractSpinArray
This type keeps the essentials of magnetic spins. Its instance struct must
contain all fields listed listed in the exemplary struct `mSpinArray`.

# Misc
Might make `AbstractSpinArray <: AbstractArray` in a future version
"""
abstract type AbstractSpinArray end

## set and get
Base.setproperty!(spa::AbstractSpinArray, s::Symbol, x) = begin
  s ∈ (:dim, :mask) && throw(error_ImmutableField(s))

  nM = prod(spa.dim)
  if (s==:M)&&(size(x,1)==1) x=repeat(x, nM, 1) end
  if (s ∈ (:T1,:T2,:γ,:M)) @assert(size(x,1)∈(1,nM)) end

  setfield!(spa, s, x)
end

## AbstractArray-like interface
Base.size(spa::AbstractSpinArray) = spa.dim
Base.size(spa::AbstractSpinArray, d) = (d ≤ length(spa.dim)) ? spa.dim[d] : 1

Base.isequal(a::AbstractSpinArray, b::AbstractSpinArray) =
  all([isequal(getproperty.((a,b),s)...) for s in fieldnames(mSpinArray)])

export mSpinArray
"""
    mSpinArray
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

# Usages
    spinarray = mSpinArray(dim::Dims, T1, T2, γ; M)
    spinarray = mSpinArray(mask::BitArray, T1, T2, γ; M)
"""
mutable struct mSpinArray <: AbstractSpinArray
  # *Immutable*:
  dim ::Dims
  mask::BitArray
  # *Mutable*:
  T1 ::TypeND(T0D, [0,1])
  T2 ::TypeND(T0D, [0,1])
  γ  ::TypeND(Γ0D, [0,1])
  M  ::TypeND(Real, [2])

  function mSpinArray(mask::BitArray, T1=1.47u"s", T2=0.07u"s", γ=γ¹H;
                      M=[0 0 1])
    dim = size(mask)
    nM = prod(dim)
    if size(M,1)==1 M=repeat(M, nM, 1) end
    @assert(all(map(x->(size(x,1) ∈ (1,nM)), [T1,T2,γ,M])))

    return new(dim, mask, T1, T2, γ, M)
  end

  mSpinArray(dim::Dims=(1,), a...; kw...) =
    mSpinArray(trues(dim), a...; kw...)

end

#= Cube =#

export AbstractSpinCube
"""
    AbstractSpinCube <: AbstractSpinArray
This type inherits `AbstractSpinArray` as a field. Its instance struct must
contain all fields listed in the exemplary struct `mSpinCube`.
"""
abstract type AbstractSpinCube <: AbstractSpinArray end

## set and get
Base.setproperty!(cb::AbstractSpinCube, s::Symbol, x) = begin
  s ∈ (:spinarray,:fov,:ofst,:loc) && throw(error_ImmutableField(s))
  s ∈ fieldnames(typeof(cb)) ? setfield!(cb, s,x) : setfield!(cb.spinarray, s,x)
end

Base.getproperty(cb::AbstractSpinCube, s::Symbol) =
  s ∈ fieldnames(typeof(cb)) ? getfield(cb, s) : getfield(cb.spinarray, s)

## AbstractArray-like interface
Base.size(cb::AbstractSpinCube, a...) = Base.size(cb.spinarray, a...)

Base.isequal(a::AbstractSpinCube, b::AbstractSpinCube) =
  all([isequal(getproperty.((a,b),s)...) for s in fieldnames(mSpinCube)])

export mSpinCube
"""
    mSpinCube
An exemplary struct instantiating `AbstractSpinCube`, designed to model a set of
regularly spaced spins, e.g., a volume.

# Fields:
*Immutable*:
- `spinarray::AbstractSpinArray` (1,): inherited `AbstractSpinArray` struct
- `fov::NTuple{3,L0D}` (3,): field of view.
- `ofst::NTuple{3,L0D}` (3,): fov offset from magnetic field iso-center.
- `loc::TypeND(L0D, [2])`  (nM, 3): location of spins.
*Mutable*:
- `Δf::TypeND(F0D, [0,1])` (1,) or (nM,): off-resonance map.

# Usages
    spincube = mSpinCube(dim::Dims{3}, fov; ofst, Δf, T1, T2, γ)
    spincube = mSpinCube(mask::BitArray{3}, fov; ofst, Δf, T1, T2, γ)
`dim`, `mask`, `T1`, `T2`, and `γ` are passed to `mSpinArray` constructors.
"""
mutable struct mSpinCube <: AbstractSpinCube
  # *Immutable*:
  spinarray::AbstractSpinArray
  fov ::NTuple{3,L0D}
  ofst::NTuple{3,L0D}
  loc ::TypeND(L0D, [2])
  # *Mutable*:
  Δf  ::TypeND(F0D, [0,1])

  function mSpinCube(mask::BitArray{3}, fov;
                     ofst=Quantity.((0,0,0), u"cm"), Δf=0u"Hz",
                     T1=1.47u"s", T2=0.07u"s", γ=γ¹H)
    spa = mSpinArray(mask, T1, T2, γ)
    loc = CartesianLocations(spa.dim)./reshape([(spa.dim./fov)...], 1,:)
    return new(spa, fov, ofst, loc, Δf)
  end

  mSpinCube(dim::Dims{3}, a...; kw...) =
    mSpinCube(trues(dim), a...; kw...)

end

#= Bolus (*Under Construction*) =#

# export AbstractSpinBolus
"""
    AbstractSpinBolus <: AbstractSpinArray
This type inherits `AbstractSpinArray` as a field. Its instance struct must
contain all fields listed in the exemplary struct `mSpinBolus`.
"""
abstract type AbstractSpinBolus <: AbstractSpinArray end

## set and get

Base.getproperty(bl::AbstractSpinBolus, s::Symbol) =
  s ∈ fieldnames(typeof(bl)) ? getfield(bl, s) : getfield(bl.spinarray, s)

## AbstractArray-like interface
Base.size(bl::AbstractSpinBolus, a...) = Base.size(bl.spinarray, a...)

Base.isequal(a::AbstractSpinBolus, b::AbstractSpinBolus) =
  all([isequal(getproperty.((a,b),s)...) for s in fieldnames(mSpinBolus)])

# export mSpinBolus
"""
    mSpinBolus
An exemplary struct instantiating `AbstractSpinBolus`, designed to model a set
of moving spins, e.g., a blood bolus in ASL context.
"""
mutable struct mSpinBolus <: AbstractSpinBolus
  # *Immutable*:
  # *Mutable*:
end

