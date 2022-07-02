
#= 𝐵-effective =#
# Methods for getting 𝐵-effective
export rfgr2B
"""
    B = rfgr2B(rf, gr, loc=[0 0 0]; Δf=0, b1Map=1, γ=γ¹H)
Turn rf, `rf`, and gradient, `gr`, into 𝐵-effective magnetic field.

*INPUTS*:
- `rf ::TypeND(Complex, (1,2))` (nT, (nCoil))
- `gr ::TypeND(Real, 2)`        (nT, 3)
- `loc::TypeND(Real, 2)`        (1,3) or (nM, 3), locations.
*KEYWORDS*:
- `Δf::TypeND(Real, (0,1,2))` (1,)  or (nM,), off-resonance.
- `b1Map::TypeND(Union{Real,Complex},(0,1,2))` (1,) or (nM,(nCoils)),
   transmit sensitivity.
- `γ::TypeND(Real, (0,1))` (1,)  or (nM,), gyro-ratio
*OUTPUS*:
- `B`, generator of `TypeND(Real, 2)` (1,1,nT), 𝐵 field.

See also: [`Pulse2B`](@ref), [`blochsim`](@ref).

# TODO:
Support `loc`, `Δf`, and `b1Map` being `Base.Generators`.
"""
function rfgr2B(rf   ::TypeND(Complex, (1,2)),
                gr   ::TypeND(Real, 2),
                loc  ::TypeND(Real, 2) = [0 0 0];
                Δf   ::TypeND(Real, (0,1)) = 0,
                b1Map::TypeND(Union{Real,Complex}, (0,1,2,3)) = 1,
                γ    ::TypeND(Real, (0,1)) = γ¹H)

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

#= 𝐴, 𝐵 =#
export B2AB
"""
    B2AB(B; T1=Inf, T2=Inf, γ=γ¹H, dt=4e-6)
Turn B-effective into Hargreave's 𝐴/𝐵, mat/vec, see: doi:10.1002/mrm.1170.

*INPUTS*:
- `B::Union{TypeND(Real, (2,3)), Base.Generator}`:
  Global, (nT,xyz); Spin-wise, (nM,xyz,nT).
*KEYWORDS*:
- `T1 & T2 ::TypeND(Real, (0,1))`: Global, (1,); Spin-wise, (nM,1).
- `γ::TypeND(Real, (0,1))`: Global, (1,); Spin-wise, (nM, 1). gyro ratio
- `dt::Real` (1,), simulation temporal step size, i.e., dwell time.
*OUTPUTS*:
- `A::TypeND(AbstractFloat, 3)` (nM, 3,3), `A[iM,:,:]` is the `iM`-th 𝐴.
- `B::TypeND(AbstractFloat, 2)` (nM, 3), `B[iM,:]` is the `iM`-th 𝐵.

See also: [`rfgr2B`](@ref), [`Pulse2B`](@ref).
"""
function B2AB(B ::Base.Generator;
              T1::TypeND(Real, (0,1))=Inf,
              T2::TypeND(Real, (0,1))=Inf,
              γ ::TypeND(Real, (0,1))=γ¹H,
              dt::Real=4e-6)

  nM = maximum([size(x,1) for x in (T1,T2,γ,first(B))])
  AB = reshape([ones(nM) zeros(nM,3) ones(nM) zeros(nM,3) ones(nM) zeros(nM,3)],
               (nM,3,4)) # as if cat(A,B;dims=3), avoid constructing U, Φ twice.
  AB1 = similar(AB)

  # in unit, convert relaxations into losses/recovery per step
  E1, E2 = exp.(-dt./T1), exp.(-dt./T2)
  E1₋₁ = E1 .- 1

  # u, ϕ = Array{Float64}(undef, nM, 3), Array{Float64}(undef, nM, 1)

  for b in B
    u, ϕ = B2UΦ(b; γ=γ, dt=dt)
    # B2UΦ!(b, u, ϕ; γ=γ, dt=dt)
    any(ϕ.!=0) && @inbounds(UΦRot!(u, view(ϕ,:,1), AB, AB1))
    AB1[:,1:2,:] .*= E2
    AB1[:,3,:]   .*= E1
    AB1[:,3,4]   .-= E1₋₁
    AB, AB1 = AB1, AB
  end

  return (A=AB[:,:,1:3], B=AB[:,:,4]) # Can one avoid the array copies here?
end

function B2AB(B ::TypeND(Real, (2,3)); kw...)
  if size(B,3) == 1 && size(B,1) != size(M,1)
    B = permutedims(B[:,:,:], [3,2,1])  # best practice?
    println("B not being spin-specific, assuming global")
  end
  return B2AB(@inbounds(view(B,:,:,t) for t in axes(B,3)); kw...)
end
