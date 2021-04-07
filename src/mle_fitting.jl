using DifferentialEquations, Optim, DiffEqFlux, DiffEqSensitivity, Plots
using Zygote, ForwardDiff, ReverseDiff
# using BlackBoxOptim#,GalacticOptim
using BenchmarkTools
using RecursiveArrayTools
# -------------------------------------------------------------------------

mutable struct F_t_Data
    F # forces vector
    t # lifetimes vector
    kd # kernel density of F
end

mutable struct M
    model   # model function
    u₀      # initial conditions
    u₀_type # string defining initial condtion type
    tspan   # time span
    dens    # site densities
    cons    # constant parameters
    f       # force
    solver  # ode solver
end


# function solve_prob(M,p;tᵢ=[])
#     f = (du,u,p,t) -> M.model(du,u,p,t,M.dens,M.cons)
#
#     # k = p[1:end-1]
#     # a = p[end]
#     u₀ = initial_u(p,M)
#
#     prob = ODEProblem(f,u₀,M.tspan,p)
#     sol = solve(prob,Rodas5(),tstops=tᵢ,dtmax=0.5,
#                             # callback=PositiveDomain(),
#                             isoutofdomain=(u,p,t)->any(x->x<0,u),
#                             abstol=1e-10,
#                             reltol=1e-10,
#                             verbose=false)
# end
function solve_prob(M,p;tᵢ=[])

    f = (du,u,p,t) -> M.model(du,u,p,t,M.dens,M.cons;f=M.f)
    u₀ = initial_u(p,M)

    prob = ODEProblem(f,u₀,M.tspan,p)
    sol = solve(prob,M.solver(),tstops=tᵢ,dtmax=0.5,
                            # callback=PositiveDomain(),
                            isoutofdomain=(u,p,t)->any(x->x<0,u),
                            abstol=1e-10,
                            reltol=1e-10,
                            verbose=false)
end

function mle_loss(p,data,M)

    L = []
    for (Fᵢ,tᵢ) in zip(data.F,data.t)

        M.f = Fᵢ
        sol = solve_prob(M,p;tᵢ)

        #check whether integration is successful
        if sol.retcode != :Success || sol.t[end] < tᵢ
            push!(L,0.0) # Likelihood of model is 0
        else

            ŷᵢ = sum(sol(tᵢ))

            mean_t = mean_lifetime_linear(sol)
            # Likelihood
            Lᵢ = ŷᵢ / mean_t
            if Lᵢ < 0
                @show Lᵢ, mean_t, ŷᵢ, tᵢ
            end

            push!(L,Lᵢ)
        end
    end
    # logL = logLikelihood(L)
    # logL = wLogLikelihood(L,1 ./ pdf(data.kd,data.F))
    # logL = wLogLikelihood(L,data.t)
    # logL = wLogLikelihood(L,log.(data.t .+ 1.0) ./ pdf(data.kd,data.F))
    logL = wLogLikelihood(L,data.t .^(0.7) ./ pdf(data.kd,data.F))
    # logL = wLogLikelihood(L,data.t .^(1))
    return logL
end

function logLikelihood(L)
    logL = -sum(log.(L))
end

function wLogLikelihood(L,w)
    wlogL = -sum(w .* log.(L))
end

I(P) = log(1/P)

using QuadGK
function mean_lifetime(sol) # numerical integration to be changed
    # <t> = ∫p(t) between [0,∞), in this case <t> = ∫sol
    # stop inegration when solution reaches 0
    endpt_id = findfirst(iszero,sol[:])
    if endpt_id == nothing
        endpt = sol.t[end]
    else
        endpt = sol.t[endpt_id]
    end

    t_meanᵤ ,err = quadgk(sol,0.0,endpt,atol=1e-10,rtol=1e-10)
    # if sum(t_meanᵤ) == 0.0
    #     @show sol
    # end
    t_mean = sum(t_meanᵤ)
end

using NumericalIntegration
function mean_lifetime_linear(sol)
    t = sol.t
    u = sum(sol,dims=1)'
    sum(integrate(t,u))

end

function optimization(M,data,options)

    # loss_in = (p) -> mle_loss(p,data,M)
    loss_in = (p) -> mle_loss(10 .^ p,data,M)

    bounds = options
    loss_in(pinit)

    # BBO optimization (use log bounds for easier search across orders of magnitude)
    # #= -----------------------------------------------------------------------
    log_bounds = map(x->log10.(x .+ 1e-5), bounds)
    res_bbo = bboptimize(loss_in; SearchRange = log_bounds,
    # res_bbo = bboptimize(loss_in; SearchRange = bounds,
                                NumDimensions = length(pinit),
                                Method = :adaptive_de_rand_1_bin,
                                # Method = :adaptive_de_rand_1_bin_radiuslimited,
                                # Method = :separable_nes,
                                NThreads=Threads.nthreads()-1,
                                MaxSteps = 1500)#,
    #
    p₀_bbo = best_candidate(res_bbo)
    p₀_bbo = 10 .^ p₀_bbo
    #-----------------------------------------------------------------------
    # =#

    # @show p₀_bbo
    # p₀_bbo = log10.(pinit)
    # p₀_bbo = pinit

    #------------------------------------------------------------------------
    #------------------------------------------------------------------------


    # then hone in on minimum with gradient based method
    loss_in = (p) -> mle_loss(p,data,M)

    # bounds = map(x->log10.(x), options)
    lb,ub = boundaries(bounds)

    # @show so.minimize(x ->loss_in(x), pinit,method="BFGS")
    # @show so.minimize(loss_in, pinit,method="BFGS")
    # @show so.newton(loss_in, pinit)#,method="Newton")
    # #=
    od = OnceDifferentiable(loss_in, p₀_bbo; autodiff = :forward);
    res = optimize(od,lb,ub,p₀_bbo,Fminbox(BFGS()),#,
    # res = optimize(loss_in,lb,ub,p₀_bbo,Fminbox(NelderMead()),#
                                        Optim.Options(show_trace=true,
                                        iterations = 50) )#,
                                        # f_tol = 1e-5,
                                        outer_iterations = 1))
    # =#
    # res = optimize(loss_in,p₀_bbo,lb,ub,Fminbox(BFGS()),
    #                                     iterations = 1,
    #                                     Optim.Options(iterations = 10))

    #= -----------------------------------------------------------------------
    optprob = GalacticOptim.OptimizationFunction((x,p) -> loss_in(x),p₀_bbo,GalacticOptim.AutoForwardDiff())
    # optprob = GalacticOptim.OptimizationFunction((x,p) -> loss_in(x))
    # # prob = GalacticOptim.OptimizationProblem(optprob,p₀_bbo)
    # # res = GalacticOptim.solve(prob,NelderMead(),cb = callback)
    #
    prob = GalacticOptim.OptimizationProblem(optprob,p₀_bbo;lb=lb,ub=ub)
    # )#
    res = GalacticOptim.solve(prob,BFGS(initial_stepnorm=0.01),
                                Fminbox(BFGS()),
                                cb = callback,
                                maxiters = 10)
    # res = GalacticOptim.solve(prob,BFGS(),cb = callback,maxiters = 20)
    =#

    rates = res.minimizer
    l = res.minimum

    hes = ForwardDiff.hessian(loss_in,rates)

    @show M.model
    @show l
    @show rates
    return rates, l, hes
end

callback = function (p, l)
  display(l)
  return false
end

function boundaries(opt)
    lb = [b[1] for b in opt]
    ub = [b[2] for b in opt]
    return lb,ub
end

using LinearAlgebra
hessian2σ(H) = sqrt.(abs.( diag(inv(H)) ))
