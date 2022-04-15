using Test, Random

using LinearAlgebra, StaticArrays

using MRphy
using MRphy: _γdt

const T = Float64

# scalars
beff = SVector{3, T}(randn(3))
u, ϕ = beff/norm(beff), -_γdt*norm(beff)

@testset "beffective tests" for _ = [1]
  u_res, ϕ_res = b2uϕ(beff, _γdt)
  um_res, ϕm_res = b2uϕ(MVector{3}(beff), _γdt)

  @test u === u_res  # SVector being bitstype
  @test (u, ϕ) == (u_res, ϕ_res)
  @test (u, ϕ) == (SVector(um_res), ϕm_res)
end
