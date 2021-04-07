
# module model_fit

using DifferentialEquations, Flux, Optim, DiffEqFlux, DiffEqSensitivity, Plots
using Zygote, ForwardDiff, ReverseDiff
using BlackBoxOptim, GalacticOptim
using BenchmarkTools
using RecursiveArrayTools
# -------------------------------------------------------------------------


mutable struct M
    model   # model function
    u₀      # initial conditions
    tspan   # time span
    dens    # site densities
    cons    # constant parameters
end

mutable struct DissData
    y
    t
    f
end


model = bimolecular_diss!
u₀ = [1.0]
tspan =  (0.0,10.0)
dens = []
cons = [1]
γδM = M(model,u₀,tspan,dens,cons)
# γδM.model = bimolecular!






function optimization(M,data,options)
# fit function - given model, data, densities, options get best parameters
    unique_t = [unique(times) for times in data.t]
    loss_in = (p) -> loss(p,data,M,unique_t)
    callback_in = (p,l) -> callback(p,l,data,M,unique_t)

    lb,ub = options

    @show loss_in(pinit)

    sr = [(0.0,1.0),(0.0,100.0)]
    res_bbo = bboptimize(loss_in; SearchRange = sr,
    # res_bbo = bboptimize(loss_in; SearchRange = options,
                                NumDimensions = length(pinit),
                                Method = :adaptive_de_rand_1_bin_radiuslimited,
                                # Method = :separable_nes,
                                NThreads=Threads.nthreads()-1,
                                MaxSteps = 20000)#,

    p₀_bbo = best_candidate(res_bbo)
    @show p₀_bbo

    res = optimize(loss_in,lb,ub,p₀_bbo, Fminbox(BFGS()), autodiff = :forward,
                                Optim.Options(show_trace=true,
                                              f_tol = 1e-5,
                                              outer_iterations = 10))

    rates = res.minimizer
    l = res.minimum

    # bootstrap_rates = bootstrapping(rates,data,model,dens,cons,t)

    @show M.model, l, rates

    hes_zyg = []
    bootstrap_rates = []

    return rates, l
end



function run_model(M,p,dens,cons,t) # version of run function with multiple models
    # f = (du,u,p,t) -> model(du,u,10 .^ p,t,densities,constants)
    f = (du,u,p,t) -> M.model(du,u,p,t,dens,cons)
    tspan = (0.0,maximum(t))
    # # tmp_prob = ODEProblem(f,M.u₀,M.tspan,p) #allows to specify tspan
    # tmp_prob = ODEProblem(f,M.u₀,tspan,p)  #automatic tspan

    # # uncomment to fit the initial conditions as well...
    # # #u₀ = [(1-p[1])/2, (1-p[1])/2, p[1]]
    u₀ = [(1-p[1])*0.610, (1-p[1])*0.390, p[1]] #ratio depends on γ and δ affinities
    p = p[2:end]
    tmp_prob = ODEProblem(f,u₀,tspan,p)
    # # end fit initial conditions
    tmp_sol = solve(tmp_prob,AutoVern7(KenCarp4()),saveat=t, abstol=1e-8,reltol=1e-8,verbose=false)
end

function loss(p,data,M,unique_t)
    Δy_all = [] # keep track of errors for all curves (change this to array of arrays!)
    # Δy_all = Zygote.Buffer([[0.,0.],[0.,0.]],length(data))
    for i in 1:length(data.y) # for each bin
      # run model with densities and parameter/dens/cons combinations
      tmp_sol = run_model(M, p, M.dens[i], M.cons[i], unique_t[i])

      if tmp_sol.retcode != :Success #if there is an error return infitite loss
        Δy_sol = fill(Inf, length(data.t[i]))
      else
        Δy_sol = find_diffs(tmp_sol,data.t[i],unique_t[i],data.y[i])
      end
      Δy_all = vcat(Δy_all,Δy_sol)
    end

    t_all = reduce(vcat,data.t)

    Losses = Δy_all .^2 .* t_all .^2
    L = sum(Losses) # calculate error across all differences from all curves

    return L
end

function find_diffs(tmp_sol,t,unique_t,y) # solutions use unique time points,
# need to evaluate the solution at each experimental time point as well as
# calculating the difference between fit and experiment

    Σ_sol = sum(Array(tmp_sol),dims=1)
    Δy = []
    for (i, tᵢ) in enumerate(t)
        ind = isequal.(unique_t,tᵢ)
        yᵢ_exp = Σ_sol[ind][1]
        yᵢ_obs = y[i]
        Δyᵢ = yᵢ_obs .- yᵢ_exp
        Δy=vcat(Δy,Δyᵢ)
    end
    return copy(Δy) #Δy[1:end] # copy(Δy)
end

# function callback(p,l,tmp_sol,data,model,dens,cons,t)
function callback(p,l,data,M,unique_t)
  # @show l, p
  @show l

  false
end
