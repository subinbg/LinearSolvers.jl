mutable struct QMRCGSTABIterable{
    opT, solT, realT<:Real, scalarT<:Number}

    A::opT
    x::solT

    rt::solT
    r::solT
    p::solT
    v::solT
    s::solT
    d::solT
    dt::solT
    xt::solT
    t::solT

    rho::scalarT
    rho_prev::scalarT
    alpha::scalarT
    omega::scalarT
    tau::scalarT
    theta::scalarT
    eta::scalarT
    beta::scalarT
    thetat::scalarT
    etat::scalarT

    maxiter::Int
    niter::Int
    tol::realT
    residual::realT
end

function init_qmrcgstab(
    x, A, b;
    abstol::Real=zero(real(eltype(b))),
    reltol::Real=sqrt(eps(real(eltype(b)))),
    maxiter::Int=2000)

    T = typeof(x)
    N = size(x)
    niter = 0

    rt = T(undef, N)
    r = T(undef, N)
    p = T(undef, N)
    v = T(undef, N)
    s = T(undef, N)
    d = T(undef, N)
    dt = T(undef, N)
    xt = T(undef, N)
    t = T(undef, N)

    fill!(p, zero(eltype(x)))
    fill!(v, zero(eltype(x)))
    fill!(d, zero(eltype(x)))

    mul!(r, A, x)
    sub!(r, b, r)
    niter += 1

    copyto!(rt, r) # rand!(rt)
    nrm = norm(r)
    tolerance = max(reltol * nrm, abstol)

    rho = one(eltype(x))
    rho_prev = one(eltype(x))
    alpha = one(eltype(x))
    omega = one(eltype(x))
    tau = nrm * one(eltype(x))
    theta = zero(eltype(x))
    eta = zero(eltype(x))
    beta = zero(eltype(x))
    thetat = zero(eltype(x))
    etat = zero(eltype(x))

    QMRCGSTABIterable(
        A, x, 
        rt, r, p, v, s, d, dt, xt, t, 
        rho, rho_prev, alpha, omega, tau, theta, eta,
        beta, thetat, etat,
        maxiter, niter, tolerance, nrm)
end


@inline converged(it::QMRCGSTABIterable) = it.residual ≤ it.tol
@inline start(::QMRCGSTABIterable) = 0
@inline done(it::QMRCGSTABIterable, iteration::Int) = it.niter ≥ it.maxiter || converged(it)



function iterate(it::QMRCGSTABIterable, iteration::Int=start(it))
    if done(it, iteration) return nothing end

    
    it.rho_prev = it.rho
    it.rho = dot(it.rt, it.r)

    it.beta = it.rho * it.alpha / (it.rho_prev * it.omega)
    sub_scmul!(it.p, it.omega, it.v)
    mul!(it.p, it.beta, it.p)
    add!(it.p, it.r, it.p)

    mul!(it.v, it.A, it.p)
    it.alpha = it.rho / dot(it.rt, it.v)
    mul!(it.s, it.alpha, it.v)
    sub!(it.s, it.r, it.s)

    # First quasi-minimization and update
    it.thetat = norm(it.s) / it.tau
    c = one(it.thetat) / sqrt(one(it.thetat) + it.thetat*it.thetat) ####
    taut = it.tau * it.thetat * c ####

    it.etat = c*c * it.alpha
    mul!(it.dt, (it.theta * it.theta * it.eta / it.alpha), it.d)
    add!(it.dt, it.p, it.dt)
    mul!(it.xt, it.etat, it.dt)
    add!(it.xt, it.x, it.xt)

    # Update r
    mul!(it.t, it.A, it.s)
    it.omega = (it.thetat * it.tau)^2 / dot(it.s, it.t) # dot(s, t) / dot(t, t)
    mul!(it.r, it.omega, it.t)
    sub!(it.r, it.s, it.r)

    # Second quasi-minimization and update
    it.theta = norm(it.r) / taut
    c = one(it.theta) / sqrt(one(it.theta) + it.theta*it.theta)
    it.tau = taut * it.theta * c

    it.eta = c*c * it.omega
    mul!(it.d, (it.thetat * it.thetat * it.etat / it.omega), it.dt)
    add!(it.d, it.s, it.d)
    mul!(it.x, it.eta, it.d)
    add!(it.x, it.xt, it.x)

    it.niter += 2
    it.residual = sqrt(iteration+1) * abs(it.tau)
    
    it.residual, iteration + 1
end



"""
- Saad, Y. (2003). Interactive method for sparse linear system.
- Freund, W. R., & Nachtigal, N. M. (1990). QMR : for a Quasi-Minimal
  Residual Linear Method Systems. (December).
- New implementation of QMR-type algorithms (2005)
- A Quasi-Minimal Residual Variant of the Bi-CGSTAB Algorithm for Nonsymmetric Systems
"""
function qmrcgstab!(
    x, A, b;
    abstol::Real=zero(real(eltype(b))),
    reltol::Real=sqrt(eps(real(eltype(b)))),
    maxiter::Int=2000,
    verbose::Bool = false)

    history = zeros(0)

    iterable = init_qmrcgstab(
        x, A, b; 
        abstol=abstol, 
        reltol=reltol, 
        maxiter=maxiter)


    for (iteration, item) = enumerate(iterable)
        push!(history, iterable.residual)
        verbose && @printf("%3d\t%1.2e\n", iteration, iterable.residual)
    end

    if !converged(iterable)
        push!(history, -1)
        verbose && println("QMRCGStab: not converged...")
    end
    # println()

    return history
end