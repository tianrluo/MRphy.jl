module MRphy

using LinearAlgebra, StaticArrays

export γ¹H
const γ¹H = 4257.6  # Gauss/hz

const _dt = 4e-6  # Sec
const _γdt = 2π*γ¹H*_dt
const _T1 = 1.47 # Sec
const _T2 = 0.07 # Sec

include("utils.jl")
include("beffective.jl")
include("sims.jl")

end # module
