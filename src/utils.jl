
export CartesianLocations
"""
    CartesianLocations(dim::Dims, doShift::Bool=true)

Retuns a `(prod(dim),length(dim))` array of grid locations. `doShift` shifts the
locations to be consistent with `fftshift`.

#Examples
```julia-repl
julia> loc = [CartesianLocations((2,2)), CartesianLocations((2,2),false)]
2-element Array{Array{Int64,2},1}:
 [-1 -1; 0 -1; -1 0; 0 0]
 [1 1; 2 1; 1 2; 2 2]
```
"""
CartesianLocations(dim::Dims, doShift::Bool=true) = doShift ?
  Array(hcat(collect.(Tuple.(CartesianIndices(dim).-ctrSub(dim)))...)') :
  Array(hcat(collect.(Tuple.(CartesianIndices(dim)))...)')

export ctrSub
"""
    ctrSub(dim::Dims) = CartesianIndex(dim .÷ 2 .+ 1)

As a separate fn, ensure consistent behaviour of getting ::CartesianIndex to the
center of a Nd-Array of size `dim`.
This `center` should match `fftshift`'s `center`.

See also: `ctrInd`
"""
ctrSub(dim::Dims) = CartesianIndex(dim .÷ 2 .+ 1)

export ctrInd
"""
    ctrInd(dim::Dims) = sum((dim.÷2) .* [1; cumprod([dim[1:end-1]...])])+1

As a separate fn, ensure consistent behariour of getting the linear index to the
center of a Nd-array of size `dim`.
This `center` should match `fftshift`'s `center`.

See also: `ctrSub`
"""
ctrInd(dim::Dims) = sum((dim.÷2) .* [1; cumprod([dim[1:end-1]...])])+1

