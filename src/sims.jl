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
