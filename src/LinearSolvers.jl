"""
A Julia package to iteratively solve linear problems
with custom matrix-vector operations.   

This package contains some codes from IterativeSolvers.jl: Copyright (c) 2013--2016 The Julia Language.   
For the copyright notice, please check the LICENSE file.
"""
module LinearSolvers

export bicgstab!, qmrcgstab!

using LinearAlgebra
using Random
using Printf

import Base: iterate

abstract type VectorLike end
abstract type MatrixLike end
abstract type ScalarLike end

mul!(x::VectorLike, A::MatrixLike, y::VectorLike) = begin x .= A * y end
mul!(x::VectorLike, y::ScalarLike, z::VectorLike) = begin x .= y .* z end
sub!(x::VectorLike, y::VectorLike, z::VectorLike) = begin x .= y .- z end
add!(x::VectorLike, y::VectorLike, z::VectorLike) = begin x .= y .+ z end
dot(x::Nothing, y::Nothing) = nothing # LinearAlgebra.dot(x,y)
norm(x::Nothing) = nothing # LinearAlgebra.norm(x)
add_scmul!(x::VectorLike, scalar::ScalarLike, y::VectorLike) = begin x .+= scalar .* y end
sub_scmul!(x::VectorLike, scalar::ScalarLike, y::VectorLike) = begin x .-= scalar .* y end
eltype(x::Nothing) = nothing # Base.eltype(x)
size(x::Nothing) = nothing # Base.size(x)
fill!(x::Nothing, scalar::ScalarLike) = nothing # Base.fill!(x, scalar)
rand!(x::Nothing) = nothing # Random.rand!(x)
copyto!(x::Nothing, y::Nothing) = nothing # Base.copyto!(x, y)


include("bicgstab.jl")
include("qmrcgstab.jl")

end
