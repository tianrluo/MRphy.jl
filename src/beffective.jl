#= 𝐵-effective =#
# Methods related to 𝐵-effective

export b2uϕ
"""
    u, ϕ = b2uϕ(b::SVector{3, T}, γdt::T=_γdt,) where {T<:Real}

  Given 𝐵-effective, `b`, compute its rotation axis/angle, `u`/`ϕ`.

  *Arguments*:
  - `b` Gauss, rotating frame 𝐵-effective;
  - `γdt` Rad/Gauss, gyro-ratio times simulation temporal step size, dt;
  *Outputs*:
  - `u`, 3D rotation axis, assumed unitary;
  - `ϕ` Rad, rotation angle;
"""
@inline function b2uϕ(
  b::SVector{3, T},
  γdt::T=_γdt,
) where {T<:Real}
  nb = norm(b)
  u, ϕ = (nb == 0 ? (b, T(0)) : (b/nb, -γdt*nb))
  return u, ϕ
end

"""
    u, ϕ = b2uϕ(b::StaticVector{3, <:Real}, γdt::Real=_γdt,)

  Given 𝐵-effective, `b`, compute its rotation axis/angle, `u`/`ϕ`.

  *Arguments*:
  - `b` Gauss, rotating frame 𝐵-effective;
  - `γdt` Rad/Gauss, gyro-ratio times simulation temporal step size, dt;
  *Outputs*:
  - `u`, 3D rotation axis, assumed unitary;
  - `ϕ` Rad, rotation angle;
"""
@inline function b2uϕ(
  b::StaticVector{3, <:Real},
  γdt::Real=_γdt,
)
  _γdt_ = eltype(b)(γdt)  # fail integer `b` here
  u, ϕ = b2uϕ(SVector(b), _γdt_)
  return typeof(b)(u), ϕ
end
