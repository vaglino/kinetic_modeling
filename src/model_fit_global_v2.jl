using Revise, BenchmarkTools
using DifferentialEquations, Optim, DiffEqFlux, DiffEqSensitivity, Plots
using Zygote, ForwardDiff, ReverseDiff
using LinearAlgebra, RecursiveArrayTools, Measurements, Statistics
using DataFrames
# using BlackBoxOptim
# -------------------------------------------------------------------------

mutable struct Data
    n # number of bonds, or p
    t # lifetimes vector
end

mutable struct M_global
    model   # model function
    u₀      # initial conditions
    # u₀_type # string defining initial condtion type
    tspan   # time span
    dens    # site densities
    cons    # constant parameters
    f       # force
    solver  # ode solver
end


function run_model(M,p,dens,cons,t) # version of run function with multiple models
    # f = (du,u,p,t) -> model(du,u,10 .^ p,t,densities,constants)
    f = (du,u,p,t) -> M.model(du,u,p,t,dens,cons)
    tmp_prob = ODEProblem(f,M.u₀,M.tspan,p)
    tmp_sol = solve(tmp_prob,M.solver(),saveat=t,
                                isoutofdomain=(u,p,t)->any(x->x<0,u),
                                abstol=1e-8,
                                reltol=1e-8,
                                verbose=false)
    # AutoVern7(KenCarp4()) works well
    # tmp_sol = solve(tmp_prob,Vern7(),saveat=t, abstol=1e-8,reltol=1e-8,verbose=false)
    tmp_sol
end


function loss(p,data,M,unique_t)
    Δy_all = [] # keep track of errors for all curves (change this to array of arrays!)
    # Δy_all = Zygote.Buffer([[0.,0.],[0.,0.]],length(data))
    for i in 1:length(data.n)
      # run model with densities and parameter/dens/cons combinations
      tmp_sol = run_model(M,p,M.dens[i],M.cons[i],unique_t[i])

      if tmp_sol.retcode != :Success #if there is an error return infitite loss
        Δy_sol = fill(Inf, length(data.n[i]))
      else
        Δy_sol = find_diffs(tmp_sol,data.t[i],unique_t[i],data.n[i])
      end
      Δy_all = vcat(Δy_all,Δy_sol)
      # Δy_all[i] = Δy_sol
    end
    # Δy = copy(Δy_all
    # Losses = [sum(Δ .^ 2) for Δ in Δy]
    # Losses = [sum( (Δy[i] .^ 2) .* t[i] .^2 ) for i in 1:length(Δy)]
    t_all = reduce(vcat,data.t)

    Losses = Δy_all .^2 .* t_all .^2
    # Losses = Δy_all .^2 .* t_all
    L = sum(Losses) # calculate error across all differences from all curves

    return L
end

function find_diffs(tmp_sol,t,unique_t,y) # solutions use unique time points,
# need to evaluate the solution at each experimental time point as well as
# calculating the difference between fit and experiment

    Σ_sol = sum(Array(tmp_sol),dims=1)
    # Σ_sol = sum(tmp_sol,dims=1)
    # Δy = Zygote.Buffer(t,length(t))
    Δy = []
    # Δy = zeros(length(t))
    for (i, tᵢ) in enumerate(t)
        ind = isequal.(unique_t,tᵢ)
        yᵢ_exp = Σ_sol[ind][1]
        yᵢ_obs = y[i]
        Δyᵢ = yᵢ_obs .- yᵢ_exp
        # Δy[i] = Δyᵢ
        # push!(Δy,Δyᵢ)
        Δy=vcat(Δy,Δyᵢ)
    end
    return copy(Δy) #Δy[1:end] # copy(Δy)
end

# function callback(p,l,tmp_sol,data,model,dens,cons,t)
function callback(p,l,data,M,unique_t)
  @show l
  # Σ_sol = run_model(model,p,data,dens,cons,unique_t)
  # fig = plot(unique_t,Σ_sol')
  # scatter!(fig,t,data,legend=:bottomright)
  # display(fig)
  false
end

function optimization(M,data,options)
# fit function - given model, data, densities, options get best parameters
  unique_t = [unique(times) for times in data.t]
  loss_in = (p) -> loss(p,data,M,unique_t)
  callback_in = (p,l) -> callback(p,l,data,M,unique_t)

  # @unpack lb,ub = options
  lb,ub = options

  # rates = [0.41751687626687767, -0.23866855433293688, 0.06689410353314687, -0.5510289313898404, -0.09606981254493874, -0.1155499030032302, -3.9999999999999964, -1.7908513448187287]
  # hes_zyg = Zygote.hessian(loss_in,rates)
  @show loss_in(pinit)
  # start = scan_param_space(pinit,data,model,dens,cons,t,unique_t,lb,ub) # !! change behavior so that it doesnt need pinit
  # start = DiffEqFlux.sciml_train(loss_in,start,ADAM(0.001),cb=callback_in,maxiters=50)
  # start = start.minimizer
  # @show start
  # res = optimize(loss_in,start,BFGS(),Optim.Options(f_tol = 1e-8))
  # res = optimize(loss_in,start,NelderMead(),Optim.Options(f_tol = 1e-8))
  res_bbo = bboptimize(loss_in; SearchRange = (lb, ub),
                                NumDimensions = length(pinit),
                                Method = :adaptive_de_rand_1_bin_radiuslimited,
                                # Method = :separable_nes,
                                NThreads=Threads.nthreads()-1,
                                MaxSteps = 10000)#,

                                # # MaxFuncEvals = 10000)#,
                                # TargetFitness=0.05)#,
                                # NThreads=Threads.nthreads()-1)
  p₀_bbo = best_candidate(res_bbo)
  @show p₀_bbo
  # p₀_bbo = pinit
  # bbo_minimizer = [0.829897, -3.39453, -3.83109, 0.215222, 0.848485, -3.08371, 0.606198, -0.906295]
  # res = optimize(loss_in,lb,ub,p₀_bbo, Fminbox(NelderMead()),
  res = optimize(loss_in,lb,ub,p₀_bbo, Fminbox(BFGS()), autodiff = :forward,
                                Optim.Options(show_trace=true,
                                              f_tol = 1e-5,
                                              outer_iterations = 10))


  # res = optimize(loss_in,p₀_bbo, LBFGS(),
  #               Optim.Options(show_trace=true,
  #                             f_tol = 1e-4))
  # res = optimize(loss_in,lb,ub,bbo_minimizer, Fminbox(NelderMead()) , Optim.Options(show_trace=true))


  # SCIML_TRAIN otpimizers
  # res = DiffEqFlux.sciml_train(loss_in,p₀_bbo,ADAM(),cb=callback_in,maxiters=100)
  # res = DiffEqFlux.sciml_train(loss_in,start,BBO(),cb=callback_in,maxiters=100)

  rates = res.minimizer
  l = res.minimum

  # bootstrap_rates = bootstrapping(rates,data,model,dens,cons,t)
  # bootstrap_rates = []

  @show M.model
  @show l
  @show rates
  # @show best_fitness(res_bbo), best_candidate(res_bbo)
  hes_for = ForwardDiff.hessian(loss_in,rates)
  # # hes_for = ForwardDiff.hessian(loss_in,rates)
  # hes_zyg = Zygote.hessian(loss_in,rates)
  # @show hes_for ≈ hes_zyg
  # hes_zyg = []
  return rates , hes_for, l
end

using LinearAlgebra
hessian2σ(H) = sqrt.(abs.( diag(inv(H)) ))


function hessian2σ(H,L,data,rates)
  m = size(data,1)
  n = size(rates,1)
  ν = m - n # degrees of freedom
  MSE = L/ν # error nomralized by degrees of freedom
  σ² = abs.( diag(inv(H)) )
  σ²_norm = σ² * MSE
  σ = sqrt.(σ²_norm)
end

# function bootstrapping(rates,data,model,dens,cons,t)
#   #create sample
#   data = n
#   B_size = 10
#   # rand_mat = rand(1:length(data),(length(data),B_size))
#
#   rand_mat = [rand(1:length(data),length(data)) for i in 1:B_size]
#
#   bootstrap_rates = zeros(size(rates))
#   # bootstrap_rates = []
#   for i in 1:B_size
#     resampled_data = data[rand_mat[i]]
#     resampled_t = t[rand_mat[i]]
#     sort_ind = sortperm(resampled_t)
#     resampled_data = resampled_data[sort_ind]
#     resampled_t = resampled_t[sort_ind]
#     unique_t = unique(resampled_t)
#
#     # loss_in = (p) -> loss_2(p,resampled_data,model,dens,cons,resampled_t,unique_t)
#     # callback_in = (p,l) -> callback(p,l,resampled_data,model,dens,cons,resampled_t,unique_t)
#     loss_in = (p) -> loss(p,data,model,dens,cons,t,unique(t))
#     callback_in = (p,l) -> callback(p,l,data,model,dens,cons,t,unique(t))
#
#     # @show loss_in(rates)
#     res = optimize(loss_in,-9, 2, rates,Fminbox(NelderMead()))
#     # @time res = DiffEqFlux.sciml_train(loss_in,rates,BFGS(),
#     #                                 lower_bounds = -9,
#     #                                 upper_bounds = 2,
#     #                                 cb = callback_in,
#     #                                 maxiters=50,
#     #                                 f_tol = 1e-4)
#
#     # @show res.minimizer
#     bootstrap_rates= hcat(bootstrap_rates,res.minimizer)
#     # push!(bootstrap_rates,res.minimizer)
#     # Σ_sol = run_model_2(model,p,resampled_data,dens,cons,unique_t)
#   end
#
#   return bootstrap_rates[:,2:end]
# end
#

# function scan_param_space(p,data,model,dens,cons,t,unique_t,lb,ub)
#       #create rand mesh
#       n = 5000
#       Δbounds = ub - lb
#       # ps = rand(n,size(p,1)) .* 10 .- 8
#       ps = rand(n,size(p,1)) .* Δbounds  .+ lb
#       # ps = rand(n,size(p,1)) .* 100
#       errors = zeros(n)
#       for i in 1:n
#           # error, a = loss(ps[i,:],data,model,dens,cons,t)
#           error = loss(ps[i,:],data,model,dens,cons,t,unique_t)
#           errors[i] = error
#           # @show error
#           if mod(i,100) == 0
#             @show i, error
#           end
#       end
#       min_ind = argmin(errors)
#       min_error = minimum(errors)
#       @show findmin(errors)
#       p_best = ps[min_ind,:]
#       @show p_best
#       if size(p,1) == 2
#           fig = scatter(ps[:,1],ps[:,2],log10.(errors),
#           zcolor = log10.(errors),camera=[0,90])
#           scatter!([p_best[1]],[p_best[2]],[log10(min_error)],markercolor=:green,markersize=5)
#           display(fig)
#       end
#       return p_best
# end
