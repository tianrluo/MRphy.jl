using Test, Random

using LinearAlgebra, StaticArrays

using MRphy
using MRphy: _γdt

const T = Float64

# scalars
beff = SVector{3, T}(randn(3))
u, ϕ = beff/norm(beff), -_γdt*norm(beff)

mi = SVector{3, T}(randn(3)) |> x -> x/norm(x)
# 𝑅 = 𝑐𝑜𝑠𝜑⋅𝐼 -(𝑐𝑜𝑠𝜑-1)⋅(𝐮𝐮ᵀ) +𝑠𝑖𝑛𝜑⋅[𝐮]ₓ, rotation matrix from 𝐮/𝜑, axis/angle
mo = cos(ϕ)*mi + ((1-cos(ϕ))*(u⋅mi))*u + sin(ϕ)*(u×mi)

@testset "utils tests" for _ in [1]
  mo_res = uϕrot(mi, u, ϕ)

  @test mo == mo_res
end

@testset "beffective tests" for _ = [1]
  u_res, ϕ_res = b2uϕ(beff, _γdt)
  um_res, ϕm_res = b2uϕ(MVector{3}(beff), _γdt)

  @test u === u_res  # SVector being bitstype
  @test (u, ϕ) == (u_res, ϕ_res)
  @test (u, ϕ) == (SVector(um_res), ϕm_res)
end
