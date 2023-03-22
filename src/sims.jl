#= Bloch simulation related
#
# `_pb` as `_pullback`.
=#
using ChainRulesCore
using StructArrays

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


struct CTX_blochsim{T}
  u ::SVector{3, T}  # b/norm(b)
  nb::T              # norm(b)
  cϕ::T              # cos(ϕ)
  sϕ::T              # sin(ϕ)
  utmi::T            # u⋅mi
end

@inline function ctx_blochsim(
  mi ::SVector{3, <:Real},
  b  ::SVector{3, T},       # Gs
  γdt::T=_γdt,                   # Rad/Gs ⇐ 2π⋅Hz/Gs⋅S
) where {T<:Real}
  # Disassemble the signature-matched `blochsim!` to reuse computations.
  # u, ϕ = b2uϕ(b, γdt) disassembled for reuse
  nb = norm(b)
  u, ϕ = (nb == 0 ? (b, T(0)) : (b/nb, -γdt*nb)) # `-γdt`: 𝐵×𝑀 → 𝑀×𝐵
  cϕ, sϕ = cos(ϕ), sin(ϕ)
  utmi = u⋅mi

  # mo = uϕrot(mi, u, ϕ) disassembled for reuse
  mo = cϕ*mi + ((1-cϕ)*utmi)*u + sϕ*(u×mi)

  ctx = CTX_blochsim(u, nb, cϕ, sϕ, utmi)
  return mo, ctx
end

@inline function ctx_blochsim_pb(
  h1,  # Duck type, stability checked by `ChainRulesTestUtils.test_rrule`
  mi ::SVector{3, <:Real},
  mo ::SVector{3, <:Real},
  γdt::T,
  ctx::CTX_blochsim{T},
) where {T}
  u, nb, cϕ, sϕ, utmi = ctx.u, ctx.nb, ctx.cϕ, ctx.sϕ, ctx.utmi
  uth1, h1xu = u⋅h1, h1×u

  h0 = cϕ*h1 + ((1-cϕ)*uth1)*u + sϕ*h1xu # uϕrot(h1, u, -ϕ)

  𝔠ϕ, 𝔰ϕ = nb == 0 ? (T(0), -γdt) : ((1-cϕ)/nb, sϕ/nb)

  nγdtdϕ = (-γdt)*(h1xu⋅mo)
  nγdtdu = 𝔰ϕ*(mi×h1) + (𝔠ϕ*utmi)*h1 + (𝔠ϕ*uth1)*mi
  db = (nγdtdϕ-u⋅nγdtdu)*u + nγdtdu  # -γdt⋅(dϕ⋅u-1/ϕ⋅(I-uuᵀ)du)

  return h0, db
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

"""
    mo = blochsim!(
      mi ::AbstractArray{SVector{3, Tmi}},
      b  ::AbstractArray{Db};  # Gs
      mh ::AbstractArray{Dmh}=similar(b),
      e1 ::Union{Nothing, Real}=nothing,
      e2 ::Union{Nothing, Real}=nothing,
      γdt::Real=_γdt,
    ) where {
      Tmi<:Real,
      Db<:AbstractVector{SVector{3, T}},
      Dmh<:AbstractVector{SVector{3, T}},
    } where {T<:Real}

  Bloch simulation.

  *Arguments*:
  - `mi` (Nd,), input spin states;
  - `b` {(Nd,) ⊻ (1,)} of (nT,) Gauss, rotating frame 𝐵-effective sequences;
  *Keyword Arguments*:
  - `mo` (Nd,), in-place output spin states under `b`;
  - `mh` (Nd,) of (nT,), in-place spin states histories under `b`;
  - `e1`, T1 relaxation coeff, exp(-dt/T1);
  - `e2`, T2 relaxation coeff, exp(-dt/T2);
  - `γdt` Rad/Gauss, gyro-ratio times simulation temporal step size, dt;
  *Outputs*:
  - `mo`, (Nd,), output spin states;

  Spin state full history, `mh`, is not explicitly returned, but only accessible
  via the input in-palce argument.
  A future version of this function may return the full history.
"""
function blochsim!(
  mi ::AbstractArray{SVector{3, Tmi}},
  b  ::AbstractArray{Db};  # Gs
  mo ::AbstractArray{SVector{3, Tmi}}=similar(mi),
  mh ::AbstractArray{Dmh}=fill(similar(b[end]), size(mi)),
  e1 ::Union{Nothing, Real}=nothing,
  e2 ::Union{Nothing, Real}=nothing,
  γdt::Real=_γdt,                      # Rad/Gs ⇐ 2π⋅Hz/Gs⋅S
) where {
  Tmi<:Real,
  Db<:AbstractVector{SVector{3, T}},
  Dmh<:AbstractVector{SVector{3, T}},
} where {T<:Real}
  @assert size(mi) == size(b) == size(mh)

  ((mi, b, mh)->blochsim!(mi, b; e1=e1, e2=e2, γdt=γdt, mh=mh)).(mi, b, mh)

  mo .= last.(mh)
  return mo
end


function ChainRulesCore.rrule(  # could we learn to live right.
  ::typeof(blochsim!),
  mi ::SVector{3, <:Real},
  b  ::AbstractVector{SVector{3, T}};  # Gs
  mh ::AbstractVector{SVector{3, T}}=similar(b),
  e1 ::Union{Nothing, Real}=nothing,
  e2 ::Union{Nothing, Real}=nothing,
  γdt::Real=_γdt,                      # Rad/Gs ⇐ 2π⋅Hz/Gs⋅S
  ctx::AbstractVector{CTX_blochsim{T}}=similar(b, CTX_blochsim{eltype(b[end])}),
  db ::AbstractVector{SVector{3, T}}=similar(b),  # Gs
) where {T<:Real}
  # Dissemble the signature-matched `blochsim!` to reuse computations.
  @assert length(b) == length(mh) == length(ctx) == length(db)
  @assert (dorelax = e1 !== nothing) == (e2 !== nothing)  # All or none

  γdt = eltype(b[end])(γdt)  # For type stability
  (_relax, _irelax, _E) = (dorelax
    ? (@inline(m -> relax(m, e1, e2)),
       @inline(m -> irelax(m, e1, e2)),
       @inline(v -> v.*(e2, e2, e1)))
    : (@inline(m -> m),
       @inline(m -> m),
       @inline(v -> v))
  )

  inds = eachindex(mh, b, db, ctx)
  _mi = mi
  @inbounds for i in inds
    mh[i], ctx[i] = ctx_blochsim(_mi, b[i], γdt)
    _mi = mh[i] = _relax(mh[i])
  end

  """
  Provided `d(mo)`, compute: `d(mi)`, `d(b)`.
  The adpulses paper denotes `{d(mo), d(mi)}` as `{h1, h0}`, respectively.
  """
  function blochsim!_pb(h1) # If we turn back time,
    @inbounds for i in @view(inds[end:-1:begin+1])
      h0 = _E(h1)  # pb of relax
      h1, db[i] = ctx_blochsim_pb(h0, mh[i-1], _irelax(mh[i]), γdt, ctx[i])
    end
    h0 = _E(h1)  # pb of relax
    h0, db[begin] = ctx_blochsim_pb(h0, mi, _irelax(mh[begin]), γdt, ctx[begin])

    return NoTangent(), h0, db  # How to @thunk and inplace???
  end

  return last(mh), blochsim!_pb
end

function ChainRulesCore.rrule(  # could we learn to live right.
  ::typeof(blochsim!),
  mi ::AbstractArray{SVector{3, Tmi}},
  b  ::AbstractArray{Db};  # Gs
  mo ::AbstractArray{SVector{3, Tmi}}=similar(mi),
  mh ::AbstractArray{Dmh}=fill(similar(b[end]), size(mi)),
  e1 ::Union{Nothing, Real}=nothing,
  e2 ::Union{Nothing, Real}=nothing,
  γdt::Real=_γdt,                      # Rad/Gs ⇐ 2π⋅Hz/Gs⋅S
  ctx::AbstractVector{Dctx}
  =fill(similar(b[end], CTX_blochsim{eltype(b[end][end])}), size(mi)),
  h0 ::AbstractArray{SVector{3, Tmi}}=similar(mi),
  db ::AbstractVector{Db}=fill(similar(b[end]), size(mi)),  # Gs
) where {
  Tmi<:Real,
  Db<:AbstractVector{SVector{3, T}},
  Dmh<:AbstractVector{SVector{3, T}},
  Dctx<:AbstractVector{CTX_blochsim{T}},
  Ddb<:AbstractVector{SVector{3, T}},
} where {T<:Real}
  @assert (size(mi) == size(b) == size(mo) == size(mh) ==
    size(ctx) == size(h0) == size(db))

  @inline(_rrule(mi, b, mh, ctx, db) =  # output only scalar pullbacks
    rrule(blochsim!, mi, b; mh=mh, e1=e1, e2=e2, γdt=γdt, ctx=ctx, db=db)[2])

  _blochsim!_pb = _rrule.(mi, b, mh, ctx, db)  # scalar pullbacks collection
  mo .= last.(mh)

  function pb_unpack(_pb, h1)
    _, _h0_, _db_ = _pb(h1)
    return (_h0=_h0_, _db=_db_)
  end
  """
  Provided `d(mo)`, compute: `d(mi)`, `d(b)`.
  The adpulses paper denotes `{d(mo), d(mi)}` as `{h1, h0}`, respectively.
  """
  function blochsim!_pb(h1) # If we turn back time,

    sa = StructArray(_h0=h0, _db=db)  # sa for struct array
    sa .= pb_unpack.(_blochsim!_pb, h1)
    return NoTangent(), sa._h0, sa._db
  end

  return mo, blochsim!_pb
end

