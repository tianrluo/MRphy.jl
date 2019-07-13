module MRphy

using LinearAlgebra

using Unitful, UnitfulMR
import Unitful: 𝐋, 𝐌, 𝐈, 𝐓

# Magnetic field strength, Frequency, Gyro ratio in SI unit dimensions
𝐁, 𝐅 = 𝐌*𝐈^-1*𝐓^-2, 𝐓^-1
𝚪 = 𝐅/𝐁

# Custom types
struct NoUnitChk end

export TypeND
"""
    TypeND(T,Ns) = Union{AbstractArray{<:T,Ns[1]}, AbstractArray{<:T,Ns[2]},...}
"""
TypeND(T, Ns) =
  Union{(map(x->x==0 ? T : AbstractArray{D, x} where {D<:T}, Ns))...}

#=
macro TypeND(T, Ns)
  return :(Union{(map(x->x==0 ? $T : AbstractArray{D,x} where{D<:$T}, $Ns))...})
end
=#

## Unitful types
export L0D, B0D, T0D, F0D, Γ0D, GR0D, RF0D
L0D,  B0D, T0D = Quantity{<:Real, 𝐋},   Quantity{<:Real, 𝐁}, Quantity{<:Real, 𝐓}
F0D,  Γ0D      = Quantity{<:Real, 𝐅},   Quantity{<:Real, 𝚪}
GR0D, RF0D     = Quantity{<:Real, 𝐁/𝐋}, Quantity{<:Union{Real, Complex}, 𝐁}

# const
export γ¹H
"""
    γ¹H
Gyromagnetic ratio of water.
"""
const γ¹H = 4257.6u"Hz/Gauss"
# Misc

# Other files
# Common structs functions must be defined before this line, so they can be
# called by the sub-scripts.

include("utils.jl")
include("mObjects.jl")
include("blochSim.jl")

end # module

