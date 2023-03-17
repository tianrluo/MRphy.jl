#= Bloch simulation related =#

export relax
"""
    mo = relax(mi::StaticVector{3, <:Real}, e1::Real, e2::Real,)

  Apply T1, T2 relaxation on `mi`, with (0., 0., 1.) as assumed equilibrium.

  *Arguments*:
  - `mi`, spin to be relaxed;
  - `e1`, T1 relaxation coeff, exp(-dt/T1);
  - `e2`, T2 relaxation coeff, exp(-dt/T2);
  *Outputs*:
  - `mo`, spin relaxed;
"""
@inline function relax(
  mi::StaticVector{3, <:Real},
  e1::Real,
  e2::Real,
)
  return typeof(mi)(mi[1]*e2, mi[2]*e2, (1-e1)+mi[3]*e1)
end

@inline function irelax(
  mo::StaticVector{3, <:Real},
  e1::Real,
  e2::Real,
)  # TODO: maybe accelerate this by passing in pre-computed reciprocal
  return typeof(mo)(mo[1]/e2, mo[2]/e2, (e1-1+mo[3])/e1)
end


export blochsim!
"""
    mo = blochsim!(
      mi ::SVector{3, <:Real},
      b  ::AbstractVector{SVector{3, T}};
      mh ::AbstractVector{SVector{3, T}}=similar(b),
      e1 ::Union{Nothing, Real}=nothing,
      e2 ::Union{Nothing, Real}=nothing,
      γdt::Real=_γdt,
    ) where {T<:Real}

  Bloch simulation.

  *Arguments*:
  - `mi`, input spin states;
  - `b` (nT,) Gauss, rotating frame 𝐵-effective sequence;
  *Keyword Arguments*:
  - `mh` (nT,), in-place spin state history under `b`;
  - `e1`, T1 relaxation coeff, exp(-dt/T1);
  - `e2`, T2 relaxation coeff, exp(-dt/T2);
  - `γdt` Rad/Gauss, gyro-ratio times simulation temporal step size, dt;
  *Outputs*:
  - `mo`, mh[end], i.e., the output spin state;

  Spin state full history, `mh`, is not explicitly returned, but only accessible
  via the input in-palce argument.
  A future version of this function may return the full history.
"""
function blochsim!(
  mi ::SVector{3, <:Real},
  b  ::AbstractVector{SVector{3, T}};  # Gs
  mh ::AbstractVector{SVector{3, T}}=similar(b),
  e1 ::Union{Nothing, Real}=nothing,
  e2 ::Union{Nothing, Real}=nothing,
  γdt::Real=_γdt,                      # Rad/Gs ⇐ 2π⋅Hz/Gs⋅S
) where {T<:Real}
  @assert length(b) == length(mh)
  @assert (dorelax = e1 !== nothing) == (e2 !== nothing)  # All or none

  γdt = eltype(b[end])(γdt)  # For type stability
  _relax = dorelax ? @inline(m -> relax(m, e1, e2)) : @inline(m -> m)
  @inbounds for i in eachindex(mh, b)
    mi = (mh[i] = _relax(uϕrot(mi, b2uϕ(b[i], γdt)...)))
  end
  return last(mh)
end
