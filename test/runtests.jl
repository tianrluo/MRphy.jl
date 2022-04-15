using Test, Random

using LinearAlgebra, StaticArrays

using MRphy
using MRphy: _γdt, _dt, _T1, _T2, irelax

const T = Float64

# scalars
e1, e2 = exp(-_dt/_T1), exp(-_dt/_T2)

beff = SVector{3, T}(randn(3))
u, ϕ = beff/norm(beff), -_γdt*norm(beff)

mi = SVector{3, T}(randn(3)) |> x -> x/norm(x)
mie = typeof(mi)(e2*mi[1], e2*mi[2], (1-e1)+e1*mi[3])  # sit and ...

# 𝑅 = 𝑐𝑜𝑠𝜑⋅𝐼 -(𝑐𝑜𝑠𝜑-1)⋅(𝐮𝐮ᵀ) +𝑠𝑖𝑛𝜑⋅[𝐮]ₓ, rotation matrix from 𝐮/𝜑, axis/angle
mo = cos(ϕ)*mi + ((1-cos(ϕ))*(u⋅mi))*u + sin(ϕ)*(u×mi)
moe = cos(ϕ)*mie + ((1-cos(ϕ))*(u⋅mie))*u + sin(ϕ)*(u×mie)

# vectors
vf(x) = [copy(x), copy(x)]

vmi, vmie = vf(mi), vf(mie)

#= tests =#
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

@testset "sims tests" for _ = [1]
  mie_res = relax(mi, e1, e2)

  vmie_res = @. relax(vmi, [e1], [e2])

  @test mie == mie_res
  @test vmie == vmie_res

  mi_res = irelax(mie, e1, e2)

  vmi_res = @. irelax(vmie, [e1], [e2])

  @test mi ≈ mi_res
  @test vmi ≈ vmi_res
end
