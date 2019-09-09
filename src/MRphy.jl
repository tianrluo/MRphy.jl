"""
# General Comments:
- `nM`, number of spins, as magnetic spin vectors are often denoted as 𝑀.
- `nT`, number of steps/time-points.

[`Unitful.jl`](https://github.com/PainterQubits/Unitful.jl) related:
- `𝐁 = 𝐌*𝐈^-1*𝐓^-2`, dimension of magnetic field strength.
- `𝐅 = 𝐓^-1`, dimension of temporal frequency.
- `𝐊 = 𝐋^-1`, dimension of spatial frequency.
- `𝚪 = 𝐅/𝐁`, dimension of gyro ratio.
"""
module MRphy

using LinearAlgebra

using Unitful, UnitfulMR
using Unitful: 𝐋, 𝐌, 𝐈, 𝐓, NoUnits

# Magnetic field strength, Frequency, Gyro ratio in SI unit dimensions
𝐁, 𝐅 = 𝐌*𝐈^-1*𝐓^-2, 𝐓^-1
𝐊, 𝚪 = 𝐋^-1, 𝐅/𝐁

# Custom types
struct NoUnitChk end # not using, saved for future

export TypeND
"""
    TypeND(T,Ns) = Union{AbstractArray{<:T,Ns[1]}, AbstractArray{<:T,Ns[2]},...}
Sugar for creating `Union`{`<:T` typed array of different dimensions}.

# Usage
*INPUTS*:
- `T::Type` (1,), the underlying type of the union.
- `Ns::Array{Int64,1}` (# diff dims,), an array of wanted dimensions.
"""
TypeND(T::Type, Ns::Array{Int64,1}) =
  Union{(map(x->x==0 ? T : AbstractArray{D, x} where {D<:T}, Ns))...}

#=
macro TypeND(T, Ns)
  return :(Union{(map(x->x==0 ? $T : AbstractArray{D,x} where{D<:$T}, $Ns))...})
end
=#

"""
    TypeND(T, ::Colon) = AbstractArray{<:T}
Sugar for creating `<:T` typed array of arbitrary dimensions.
"""
TypeND(T::Type, ::Colon) = AbstractArray{<:T}

## Unitful types
export B0D, Γ0D, L0D, K0D, T0D, F0D, GR0D, RF0D
"""
    B0D = Quantity{<:Real, 𝐁}
Type of magetic field strength. Based on
[`Unitful.Quantity`](https://github.com/PainterQubits/Unitful.jl).

# Examples:
```julia-repl
julia> (1u"Gauss")::B0D
1 Gauss
```
"""
B0D = Quantity{<:Real, 𝐁}

"""
    Γ0D = Quantity{<:Real, 𝚪}
Type of gyro magnetic ratio. Based on
[`Unitful.Quantity`](https://github.com/PainterQubits/Unitful.jl).

# Examples:
```julia-repl
julia> (1u"Hz/Gauss")::Γ0D
1 Hz Gauss^-1
```
"""
Γ0D = Quantity{<:Real, 𝚪}

"""
    L0D = Quantity{<:Real, 𝐋}
Type of length. Based on
[`Unitful.Quantity`](https://github.com/PainterQubits/Unitful.jl).

# Examples:
```julia-repl
julia> (1u"cm")::L0D
1 cm
```
"""
L0D = Quantity{<:Real, 𝐋}

"""
    K0D =  Quantity{<:Real, 𝐊}
Type of spatial frequency. Based on
[`Unitful.Quantity`](https://github.com/PainterQubits/Unitful.jl).

# Examples:
```julia-repl
julia> (1u"cm^-1")::K0D
1 cm^-1
```
"""
K0D =  Quantity{<:Real, 𝐊}

"""
    T0D = Quantity{<:Real, 𝐓}
Type of time. Based on
[`Unitful.Quantity`](https://github.com/PainterQubits/Unitful.jl).

# Examples:
```julia-repl
julia> (1u"s")::T0D
1 s
```
"""
T0D = Quantity{<:Real, 𝐓}

"""
    F0D =  Quantity{<:Real, 𝐅}
Type of temporal frequency. Based on
[`Unitful.Quantity`](https://github.com/PainterQubits/Unitful.jl).

# Examples:
```julia-repl
julia> (1u"s^-1")::F0D
1 s^-1
```
"""
F0D =  Quantity{<:Real, 𝐅}

"""
    GR0D = Quantity{<:Real, 𝐁/𝐋}
Type of magnetic gradient. Based on
[`Unitful.Quantity`](https://github.com/PainterQubits/Unitful.jl).

# Examples:
```julia-repl
julia> (1u"Gauss/cm")::GR0D
1 Gauss cm^-1
```
"""
GR0D = Quantity{<:Real, 𝐁/𝐋}

"""
    RF0D = Quantity{<:Union{Real, Complex}, 𝐁}
Type of magnetic RF. Based on
[`Unitful.Quantity`](https://github.com/PainterQubits/Unitful.jl).

# Examples:
```julia-repl
julia> ((1+1im)u"Gauss")::RF0D
1 + 1im Gauss
```
"""
RF0D = Quantity{<:Union{Real, Complex}, 𝐁}

# const
export γ¹H
"""
    const γ¹H = 4257.6u"Hz/Gauss"
Gyromagnetic ratio of water proton.
"""
const γ¹H = 4257.6u"Hz/Gauss"

# Misc

# Other files
# Common structs functions must be defined before this line, so they can be
# called by the sub-scripts.

include("utils.jl")
using .utils
include("SteadyStates.jl")
include("blochSim.jl")
include("mObjects.jl")

end # module

