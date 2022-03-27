"""
Some utilities functions routinely used in MR simulations.
"""
module utils

using ..MRphy

#= General =#
export CartesianLocations
"""
    CartesianLocations(dim::Dims, doShift::Bool=true)

Retuns a `(prod(dim),length(dim))` sized array of grid locations. `doShift`
shifts the locations to be consistent with `fftshift`.

# Examples
```julia-repl
julia> loc = [CartesianLocations((2,2)), CartesianLocations((2,2),false)]
2-element Array{Array{Int,2},1}:
 [-1 -1; 0 -1; -1 0; 0 0]
 [1 1; 2 1; 1 2; 2 2]
```
"""
CartesianLocations(dim::Dims, doShift::Bool=true) = (doShift
  ? Array(hcat(collect.(Tuple.(CartesianIndices(dim).-ctrSub(dim)))...)')
  : Array(hcat(collect.(Tuple.(CartesianIndices(dim)))...)'))

export ctrSub
"""
    ctrSub(dim::Dims) = CartesianIndex(dim .÷ 2 .+ 1)

As a separate function, ensure consistent behaviour of getting ::CartesianIndex
to the center of a Nd-Array of size `dim`.
This `center` should match `fftshift`'s `center`.

See also: [`ctrInd`](@ref).

# Notes:
The function may be removed once julia FFT packages provides this functionality.
"""
ctrSub(dim::Dims) = CartesianIndex(dim .÷ 2 .+ 1)

export ctrInd
"""
    ctrInd(dim::Dims) = sum((dim.÷2) .* [1; cumprod([dim[1:end-1]...])])+1

As a separate fn, ensure consistent behaviour of getting the linear index to the
center of a Nd-array of size `dim`.
This `center` should match `fftshift`'s `center`.

See also: [`ctrSub`](@ref).
"""
ctrInd(dim::Dims) = sum((dim .÷ 2) .* [1; cumprod([dim[1:end-1]...])])+1

#= MR =#
export k2g
"""
    k2g(k, isTx::Bool=false; dt::Real=4e-6, γ::Real=γ¹H)
Gradient, `g`, of the `TxRx` k-space, (trasmit/receive, excitation/imaging).

# Usage
*INPUTS*:
- `k` (nSteps, Nd...), cm⁻¹ Tx or Rx k-space, w/ unit.
- `isTx::Bool`, if `true`, compute transmit k-space, `k`, ends at the origin.
*KEYWORDS*:
- `dt::Real` (1,), Sec, gradient temporal step size, i.e., dwell time.
- `γ::Real` (1,), Hz/Gauss, gyro-ratio.
*OUTPUTS*:
- `g` (nSteps, Nd...), Gauss/cm, gradient.

# Note
The function asserts if `k` ends at the origin for `isTx==true`.

See also: [`g2k`](@ref), [`g2s`](@ref).
"""
k2g(k, isTx::Bool=false; dt::Real=4e-6, γ::Real=γ¹H) =
  ((isTx && any(selectdim(k, 1, size(k, 1)) .!= 0)) 
   ? error("Tx `k` must end at 0")
   : [selectdim(k, 1, 1:1); diff(k, dims=1)]/(γ*dt))

export g2k
"""
    g2k(g; isTx::Bool=false, dt::Real=4e-6, γ::Real=γ¹H)
Compute k-space from gradient.

# Usage
*INPUTS*:
- `g` (nSteps, Nd...), Gauss/cm, gradient
- `isTx::Bool`, if `true`, compute transmit k-space, `k`, ends at the origin.
*KEYWORDS*:
- `dt::Real` (1,), Sec, gradient temporal step size, i.e., dwell time.
- `γ::Real` (1,), Hz/Gauss, gyro-ratio.
*OUTPUTS*:
- `k` (nSteps, Nd...), cm⁻¹, k-space.

See also: [`k2g`](@ref), [`g2s`](@ref).
"""
g2k(g, isTx::Bool=false; dt::Real=4e-6, γ::Real=γ¹H) =
  γ*dt*cumsum(g,dims=1) |> k->isTx ? k.-selectdim(k,1,size(k,1):size(k,1)) : k

export g2s
"""
    g2s(g; dt::Real=4e-6)

Slew rate `sl`, of the gradient, `g`.

# Usage
*INPUTS*:
- `g` (nSteps, Nd...), Gauss/cm
*KEYWORDS*:
- `dt::Real` (1,), Sec, gradient temporal step size, i.e., dwell time.
*OUTPUTS*:
- `sl` (nSteps, Nd...), Gauss/cm/Sec, slew rate

See also: [`g2k`](@ref), [`k2g`](@ref).

# Note
No `s2g` is provided for the moment.
"""
g2s(g; dt::Real=4e-6) = [selectdim(g, 1, 1:1); diff(g, dims=1)]/dt

end # module utils
