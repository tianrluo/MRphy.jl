"""
# General Comments:
- `nM`, number of spins, as magnetic spin vectors are often denoted as 𝑀.
- `nT`, number of steps/time-points.
"""
module MRphy

using LinearAlgebra

export TypeND
"""
    TypeND(T,Ns) = Union{AbstractArray{<:T,Ns[1]}, AbstractArray{<:T,Ns[2]},...}
Sugar for creating `Union`{`<:T` typed array of different dimensions}.

# Usage
*INPUTS*:
- `T::Type` (1,), the underlying type of the union.
- `Ns::Dims` (ndims,), a tuple of wanted dimensions.
"""
TypeND(T::Type, Ns::Dims) =
  Union{(map(x->x==0 ? T : AbstractArray{D, x} where {D<:T}, Ns))...}

"""
    TypeND(T::Type, N::Int) = AbstractArray{<:T, N}
sugar of consistency
"""
TypeND(T::Type, N::Int) = AbstractArray{<:T, N}

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

# const
export γ¹H
"""
    const γ¹H = 4257.6
Gyromagnetic ratio of water proton, Hz/Gauss
"""
const γ¹H = 4257.6

# Misc

# Other files
# Common structs functions must be defined before this line, so they can be
# called by the sub-scripts.

include("utils.jl")
using .utils
include("SteadyStates.jl")
include("beffective.jl")
include("sims.jl")
include("mobjs.jl")

end # module

