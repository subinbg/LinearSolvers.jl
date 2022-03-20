"""
This package contains some codes from IterativeSolvers.jl: Copyright (c) 2013--2016 The Julia Language.   
For the copyright notice, please check the LICENSE file.
"""
mutable struct BiCGStabIterable{
    opT, solT, vecT, matT, collectT, 
    realT <: Real, 
    scalarT <: Number}

    A::opT
    l::Int

    x::solT
    r_shadow::solT
    rs::collectT
    us::collectT

    maxiter::Int
    niter::Int
    tol::realT
    residual::realT

    gamma::vecT
    omega::scalarT
    sigma::scalarT
    M::matT

    function BiCGStabIterable(
        A::opT, l::Int, x::solT,
        r_shadow::solT, rs::collectT, us::collectT,
        maxiter::Int, niter::Int, tol::realT, residual::realT,
        gamma::vecT, omega::scalarT, sigma::scalarT,
        M::matT) where {opT,solT,vecT,matT,collectT,realT<:Real,scalarT<:Number}

        new{opT,solT,vecT,matT,collectT,realT,scalarT}(
            A,l,x,r_shadow,rs,us,maxiter,niter,tol,residual,gamma,omega,sigma,M)
    end
end


function init_bicgstab(
    x, A, b, l::Int=2;
    abstol::Real=zero(real(eltype(b))),
    reltol::Real=sqrt(eps(real(eltype(b)))),
    maxiter::Int=2000)

    T = typeof(x)
    N = size(x)
    niter = 0

    # Large vectors.
    r_shadow = T(undef, N)
    rand!(r_shadow)
    rs = Array{T}(undef, l+1)
    us = Array{T}(undef, l+1)
    for idx in 1:(l+1)
        rs[idx] = T(undef, N)
        us[idx] = T(undef, N)
        fill!(us[idx], zero(eltype(x)))
    end

    residual = rs[1]

    # Compute the initial residual = b - A * x
    mul!(residual, A, x)
    sub!(residual, b, residual)
    niter += 1

    gamma = zeros(eltype(x), l)
    omega = sigma = one(eltype(x))

    nrm = norm(residual)

    # For the least-squares problem
    M = zeros(eltype(x), l+1, l+1)

    # Stopping condition based on absolute and relative tolerance.
    tolerance = max(reltol * nrm, abstol)

    BiCGStabIterable(A, l, x, r_shadow, rs, us,
        maxiter, niter, tolerance, nrm,
        gamma, omega, sigma, M
    )
end

@inline converged(it::BiCGStabIterable) = it.residual ≤ it.tol
@inline start(::BiCGStabIterable) = 0
@inline done(it::BiCGStabIterable, iteration::Int) = it.niter ≥ it.maxiter || converged(it)



function iterate(it::BiCGStabIterable, iteration::Int=start(it))
    if done(it, iteration) return nothing end

    T = eltype(it.x)
    L = 2:(it.l + 1)

    it.sigma = -it.omega * it.sigma

    ## BiCG part
    for j = 1:(it.l)
        rho = dot(it.r_shadow, it.rs[j])
        beta = rho / it.sigma

        # us[:, 1 : j] .= rs[:, 1 : j] - β * us[:, 1 : j]
        for idx = 1:j
            sub!(it.us[idx], it.rs[idx], it.us[idx])
        end

        # us[:, j + 1] = Pl \ (A * us[:, j])
        mul!(it.us[j+1], it.A, it.us[j])

        it.sigma = dot(it.r_shadow, it.us[j+1])
        alpha = rho / it.sigma

        for idx = 1:j
            sub_scmul!(it.rs[idx], alpha, it.us[idx+1])
        end

        # rs[:, j + 1] = Pl \ (A * rs[:, j])
        mul!(it.rs[j+1], it.A , it.rs[j])

        # x = x + α * us[:, 1]
        add_scmul!(it.x, alpha, it.us[1])
    end

    # Bookkeeping
    it.niter += 2 * it.l

    ## MR part
    # M = rs' * rs
    for cidx in CartesianIndices(it.M)
        i1 = cidx[1]
        i2 = cidx[2]
        it.M[i1,i2] = dot(it.rs[i1], it.rs[i2])
    end

    # γ = M[L, L] \ M[L, 1]
    F = lu!(view(it.M, L, L))
    ldiv!(it.gamma, F, view(it.M, L, 1))

    for idx in L
        sub_scmul!(it.us[1], it.gamma[idx-1], it.us[idx])
        add_scmul!(it.x, it.gamma[idx-1], it.rs[idx-1])
        sub_scmul!(it.rs[1], it.gamma[idx-1], it.rs[idx])
    end

    it.omega = it.gamma[it.l]
    it.residual = norm(it.rs[1])

    it.residual, iteration + 1
end

"""
l: number of GMRES steps
"""
function bicgstab!(
    x, A, b, l::Int=2;
    abstol::Real=zero(real(eltype(b))),
    reltol::Real=sqrt(eps(real(eltype(b)))),
    maxiter::Int=2000,
    verbose::Bool=false)

    history = zeros(0)

    iterable = init_bicgstab(
        x, A, b, l;
        abstol=abstol, reltol=reltol,
        maxiter=maxiter)

    for (iteration, item) = enumerate(iterable)
        push!(history, iterable.residual)
        verbose && @printf("%3d\t%1.2e\n", iteration, iterable.residual)
    end

    if !converged(iterable)
        push!(history, -1)
        verbose && println("BiCGStab: not converged...")
    end
    # println()

    return history
end



