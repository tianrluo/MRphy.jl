#= Utils routinely used in MR simulations =#

export uϕrot
"""
    vo = uϕrot(v::StaticVector{3, <:Real}, u::StaticVector{3, <:Real}, ϕ::Real,)

  Apply axis-angle, `u,ϕ`-based rotation on `v`.

  *Arguments*:
  - `v`, vector to be rotated;
  - `u`, 3D rotation axis, assumed unitary;
  - `ϕ` Rad, rotation angle;
  *Outputs*:
  - `vo`, vector rotated;
"""
@inline function uϕrot(
  v::StaticVector{3, <:Real},
  u::StaticVector{3, <:Real},
  ϕ::Real,
)
  # 𝑅 = 𝑐𝑜𝑠𝜑⋅𝐼 -(𝑐𝑜𝑠𝜑-1)⋅(𝐮𝐮ᵀ) +𝑠𝑖𝑛𝜑⋅[𝐮]ₓ, rotation matrix from 𝐮/𝜑, axis/angle
  return cos(ϕ)*v + ((1-cos(ϕ))*(u⋅v))*u + sin(ϕ)*(u×v)
end
