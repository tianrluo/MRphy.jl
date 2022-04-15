module MRphy

using LinearAlgebra, StaticArrays

export γ¹H
const γ¹H = 4257.6  # Gauss/hz

const _dt = 4e-6  # Sec
const _γdt = 2π*γ¹H*_dt

include("utils.jl")
include("beffective.jl")

end # module
