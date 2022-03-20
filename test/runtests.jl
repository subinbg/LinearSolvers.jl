using LinearSolvers
using Test
using LinearAlgebra
using Random

mutable struct Custom{T}
    data::T

    function Custom{T}(::UndefInitializer, D::Dims) where T
        new{T}(T(undef, D))
    end

    function Custom(data::T) where T
        new{T}(data)
    end
end


LinearSolvers.mul!(x::Custom, A::Custom, y::Custom) = begin x.data .= A.data * y.data end
LinearSolvers.mul!(x::Custom, y::Number, z::Custom) = begin x.data .= y .* z.data end
LinearSolvers.sub!(x::Custom, y::Custom, z::Custom) = begin x.data .= y.data .- z.data end
LinearSolvers.add!(x::Custom, y::Custom, z::Custom) = begin x.data .= y.data .+ z.data end
LinearSolvers.dot(x::Custom, y::Custom) = dot(x.data, y.data)
LinearSolvers.norm(x::Custom) = norm(x.data)
LinearSolvers.add_scmul!(x::Custom, scalar::Number, y::Custom) = begin x.data .+= scalar .* y.data end
LinearSolvers.sub_scmul!(x::Custom, scalar::Number, y::Custom) = begin x.data .-= scalar .* y.data end
LinearSolvers.eltype(x::Custom) = eltype(x.data)
LinearSolvers.size(x::Custom) = size(x.data)
LinearSolvers.fill!(x::Custom, scalar::Number) = fill!(x.data, scalar)
LinearSolvers.rand!(x::Custom) = rand!(x.data)
LinearSolvers.copyto!(x::Custom, y::Custom) = copyto!(x.data, y.data)

@testset "LinearSolvers.jl" begin
    N = 20
    for solver! in (bicgstab!, qmrcgstab!)
        for T in (Float32, Float64, ComplexF32, ComplexF64)
            A = Custom(rand(T, N, N) + 15I)
            x = Custom(ones(T, N))
            b = Custom(A.data * x.data)
            LinearSolvers.rand!(x)

            history = solver!(x, A, b, verbose=false, maxiter=100)
            @test history[end] > 0 # Convergence
            @test norm(A.data * x.data .- b.data) / norm(b.data) <= sqrt(eps(real(LinearSolvers.eltype(b))))
        end
    end
end
