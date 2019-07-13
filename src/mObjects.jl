include("utils.jl")

error_ImmutableField(x) = error("`$x` is an immutable field." )

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
- `rf::TypeND(RF0D, [1,2])` (nSteps,) or (nSteps, nCoils).
- `gr::TypeND(GR0D, [2])` (nSteps, 3), where 3 accounts for x-y-z channels.
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
    if (rf_miss && gr_miss) error("Missing both inputs.")       end
    if rf_miss              rf = zeros(size(gr,1))u"Gauss"      end
    if gr_miss              gr = zeros(size(rf,1),3)u"Gauss/cm" end
    @assert(size(gr,2) == 3)

    if isa(rf, Number) rf = [rf] end
    return new(rf, gr, dt, des)
  end

  function Pulse(::NoUnitChk,
                 rf=missing, gr=missing; dt=4e-6, des="generic pulse")
    if !ismissing(rf) rf = Quantity.(rf, u"Gauss")    end
    if !ismissing(gr) gr = Quantity.(gr, u"Gauss/cm") end
    return Pulse(rf, gr; dt=Quantity(dt, u"s"), des=des)
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
  if (s ∈ (:dim, :mask)) error_ImmutableField(s) end

  nSpins = prod(spa.dim)
  if (s==:M)&&(size(x,1)==1) x=repeat(x, nSpins, 1) end
  if (s ∈ (:T1,:T2,:γ,:M)) @assert(size(x,1)∈(1,nSpins)) end

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
- `dim::Dims` (nd,): `nSpins ← prod(dim)`, dimension of the object.
- `mask::BitArray` (nx,(ny,(nz))): Mask for `M`, `dim == (nx,(ny,(nz)))`
*Mutable*:
- `T1::TypeND(T0D, [0,1])` (1,) or (nSpins,): Longitudinal relaxation coeff.
- `T2::TypeND(T0D, [0,1])` (1,) or (nSpins,): Transversal relaxation coeff.
- `γ::TypeND(Γ0D, [0,1])`  (1,) or (nSpins,): Gyromagnetic ratio.
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
    nSpins = prod(dim)
    if size(M,1)==1 M=repeat(M, nSpins, 1) end
    @assert(all(map(x->(size(x,1) ∈ (1,nSpins)), [T1,T2,γ,M])))

    return new(dim, mask, T1, T2, γ, M)
  end

  mSpinArray(dim::Dims=(1,), args...; kargs...) =
    mSpinArray(trues(dim), args...; kargs...)

  # Impish pirate bypassing unit check, ho ho ho.
  function mSpinArray(::NoUnitChk,
                      mask::BitArray, T1=1.47, T2=0.07, γ=ustrip(γ¹H);
                      M=[0 0 1])
    T1, T2, γ = Quantity.(T1,u"s"), Quantity.(T2,u"s"), Quantity.(γ,u"Hz/Gauss")
    return mSpinArray(mask, T1, T2, γ, M=M)
  end

  mSpinArray(::NoUnitChk, dim::Dims=(1,), args...; kargs...) =
    mSpinArray(NoUnitChk(), trues(dim), args...; kargs...)

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
  if s ∈ (:spinarray,:fov,:ofst,:loc) error_ImmutableField(s) end

  if s ∈ fieldnames(typeof(cb.spinarray)) setfield!(cb.spinarray, s, x)
  else                                    setfield!(cb, s, x)
  end
end

Base.getproperty(cb::AbstractSpinCube, s::Symbol) = begin
  # No `typeof(cb.spinarray)` here, it causes self-call then stackoverflow.
  if s ∈ fieldnames(typeof(getfield(cb, :spinarray))) getfield(cb.spinarray, s)
  else                                                getfield(cb, s)
  end
end

## AbstractArray-like interface
Base.size(cb::AbstractSpinCube, args...) = Base.size(cb.spinarray, args...)

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
- `loc::TypeND(L0D, [2])`  (nSpins, 3): location of spins.
*Mutable*:
- `Δf::TypeND(F0D, [0,1])` (1,) or (nSpins,): off-resonance map.

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

  mSpinCube(dim::Dims{3}, args...; kargs...) =
    mSpinCube(trues(dim), args...; kargs...)

  function mSpinCube(::NoUnitChk,
                     mask::BitArray{3}, fov;
                     ofst=(0,0,0), Δf=0,
                     T1=1.47, T2=0.07, γ=ustrip(γ¹H))
    fov       = Quantity.(fov,u"cm")
    ofst,  Δf = Quantity.(ofst,u"cm"),Quantity.(Δf,u"Hz")
    T1, T2, γ = Quantity.(T1,u"s"), Quantity.(T2,u"s"), Quantity.(γ,u"Hz/Gauss")
    return mSpinCube(mask, fov, ofst=ofst, Δf=Δf, T1=T1, T2=T2, γ=γ)
  end

  mSpinCube(::NoUnitChk, dim::Dims{3}, args...; kargs...) =
    mSpinCube(NoUnitChk(), trues(dim), args...; kargs...)

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

Base.getproperty(bl::AbstractSpinBolus, s::Symbol) = begin
  # No `typeof(bl.spinarray)` here, it causes self-call then stackoverflow.
  if s ∈ fieldnames(typeof(getfield(bl, :spinarray))) getfield(bl.spinarray, s)
  else                                                getfield(bl, s)
  end
end

## AbstractArray-like interface
Base.size(bl::AbstractSpinBolus, args...) = Base.size(bl.spinarray, args...)

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

