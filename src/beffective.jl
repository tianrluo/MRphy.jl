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

export rfgr2b
"""
    B_gen = rfgr2b(
      rf   ::AbstractVector{Complex{D}},
      gr   ::AbstractVector{D},
      loc  ::StaticVector{3, T};
      Δf   ::Real=0.,
      γ    ::Real=γ¹H,
    ) where {D<:Real, T<:Real}
"""
function rfgr2b(
  rf   ::AbstractVector{Complex{D}},
  gr   ::AbstractVector{D},
  loc  ::StaticVector{3, T};
  Δf   ::Real=0.,
  γ    ::Real=γ¹H,
) where {D<:Real, T<:Real}

  nM = maximum(map(x->size(x,1), (loc, Δf, b1Map, γ)))

  Bxy_gen = b1Map==1 ?
    @inbounds(view(rf,t,:)       |> x->repeat([real(x) imag(x)], nM)
              for t in axes(rf,1)) :
    @inbounds(b1Map*view(rf,t,:) |> x->       [real(x) imag(x)]
              for t in axes(rf,1))

  Bz_gen = Δf==0 ?
    @inbounds(loc*view(gr,t,:)          for t in axes(gr,1)) :
    @inbounds(loc*view(gr,t,:).+(Δf./γ) for t in axes(gr,1))

  B_gen = @inbounds([bxy bz] for (bxy, bz) in zip(Bxy_gen, Bz_gen))
  return B_gen
end


function rfgr2b(
  rf   ::AbstractVector{Complex},
  gr   ::AbstractVector{Real},
  loc  ::StaticVector{3, T};
  Δf   ::Real=0.,
  b1Map::Real=0.,
  γ    ::Real=γ¹H,
) where {T<:Real}
  throw("Not Implemented")
end
