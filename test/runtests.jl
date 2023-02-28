using Test, Random
# using InteractiveUtils  # @code_warntype

using LinearAlgebra, StaticArrays, ChainRulesCore
using ChainRulesCore: rrule
using ChainRulesTestUtils

using MRphy
using MRphy: _γdt, _dt, _T1, _T2, irelax

const T = Float64
const T_alter = Float32

# scalars
e1, e2 = exp(-_dt/_T1), exp(-_dt/_T2)

beff = SVector{3, T}(randn(3))
u, ϕ = beff/norm(beff), -_γdt*norm(beff)
b = [beff]
b2 = [beff, beff]

mi = SVector{3, T}(randn(3)) |> x -> x/norm(x)
mie = typeof(mi)(e2*mi[1], e2*mi[2], (1-e1)+e1*mi[3])  # sit and ...

# 𝑅 = 𝑐𝑜𝑠𝜑⋅𝐼 -(𝑐𝑜𝑠𝜑-1)⋅(𝐮𝐮ᵀ) +𝑠𝑖𝑛𝜑⋅[𝐮]ₓ, rotation matrix from 𝐮/𝜑, axis/angle
mo = cos(ϕ)*mi + ((1-cos(ϕ))*(u⋅mi))*u + sin(ϕ)*(u×mi)
moe = typeof(mo)(e2*mo[1], e2*mo[2], (1-e1)+e1*mo[3])  # sit and ...

# vectors
vf(x) = [copy(x), copy(x)]

vmi, vmie = vf(mi), vf(mie)
vmo, vmoe = vf(mo), vf(moe)

vb = vf(b)

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

  mh_res, mhe_res = (similar([mi]) for _ in 1:2)
  mo_res = blochsim!(mi, b; γdt=_γdt, mh=mh_res)
  moe_res = blochsim!(mi, b; e1=e1, e2=e2, γdt=_γdt, mh=mhe_res)

  @test mo_res === mh_res[end]
  @test mo == mh_res[end]

  @test moe_res === mhe_res[end]
  @test moe == mhe_res[end]

  # @code_warntype blochsim!(mi, b; mh=mh_res, e1=e1, e2=e2, γdt=_γdt,)
  # @code_warntype blochsim!(mi, b; mh=mh_res, γdt=_γdt,)

  vmh_res, vmhe_res = (vf(similar([mi])) for _ in 1:2)
  vmo_res = similar(vmi)
  blochsim!(vmi, vb; mo=vmo_res, γdt=_γdt, mh=vmh_res)
  vmoe_res = blochsim!(vmi, vb; e1=e1, e2=e2, γdt=_γdt, mh=vmhe_res)

  @test vmo_res == last.(vmh_res)
  @test vmo == last.(vmh_res)

  @test vmoe_res == last.(vmhe_res)
  @test vmoe == last.(vmhe_res)

  # @code_warntype blochsim!(vmi, vb; mh=vmh_res, e1=e1, e2=e2, γdt=_γdt,)
  # @code_warntype blochsim!(vmi, vb; mh=vmh_res, γdt=_γdt,)
end

@testset "rrules tests" for _ = [1]
  test_rrule(blochsim!, mi, b; fkwargs=(γdt=T_alter(_γdt),))
  test_rrule(blochsim!, mi, b; fkwargs=(e1=e1, e2=e2, γdt=T_alter(_γdt),))
  test_rrule(blochsim!, mi, b2; fkwargs=(γdt=T_alter(_γdt),))
  test_rrule(blochsim!, mi, b2; fkwargs=(e1=e1, e2=e2, γdt=T_alter(_γdt),))

  mo_res, _ = rrule(blochsim!, mi, b; γdt=_γdt)
  moe_res, _ = rrule(blochsim!, mi, b; e1=e1, e2=e2, γdt=_γdt)

  @test mo ≈ mo_res
  @test moe ≈ moe_res
end

