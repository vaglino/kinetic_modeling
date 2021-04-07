
# module model_fit

using DifferentialEquations, Flux, Optim, DiffEqFlux, DiffEqSensitivity, Plots
using BlackBoxOptim
using BenchmarkTools
using RecursiveArrayTools # for VectorOfArray
# -------------------------------------------------------------------------


function generate_data(sol)
  t = collect(range(0,stop=10,length=101))

  solution = sum(sol(t),dims=1)
  # randomized = VectorOfArray([(sol(t[i]) + .02randn(size(sol,1))) for i in 1:length(t)])
  # dataset = convert(Array,randomized)
  randomized = solution .+ .01randn(size(solution))
  dataset = randomized
  display(scatter!(t,dataset'))
  return dataset
end

function generate_data(sol,t)

  solution = sum(sol(t),dims=1)
  # randomized = VectorOfArray([(sol(t[i]) + .02randn(size(sol,1))) for i in 1:length(t)])
  # dataset = convert(Array,randomized)
  randomized = solution .+ .01randn(size(solution))
  dataset = randomized
  display(scatter!(t,dataset'))
  return dataset
end


function run_model_2(model,p,data,densities,constants,t) # version of run function with multiple models
    f = (du,u,p,t) -> model(du,u,10 .^ p,t,densities,constants)
    tmp_prob = ODEProblem(f,u₀,tspan,p)
    # tmp_sol = solve(tmp_prob,AutoVern7(KenCarp4()),saveat=t, abstol=1e-8,reltol=1e-8)
    tmp_sol = solve(tmp_prob,Vern7(),saveat=t, abstol=1e-5,reltol=1e-5)
    # stack Σ of solution across n species
    tmp_sol
end


function loss_2(p,data,model,dens,cons,t,unique_t)
    # for i in 1:siz
    tmp_sol = run_model_2(model,p,data,dens,cons,unique_t)

    # @show Σ_sol
    # @show solution
    # if any((s.retcode != :Success for s in solution))
    # @show typeof(solution.retcode)
    if tmp_sol.retcode != :Success
      # @show solution
      L = Inf
    else
      Σ_sol = sum(Array(tmp_sol),dims=1)
      Δy = 0.
      for (i, tᵢ) in enumerate(t)
          ind = unique_t .== tᵢ
          yᵢ_exp = Σ_sol[ind][1]
          yᵢ_obs = data[i]
          Δyᵢ = yᵢ_obs .- yᵢ_exp
          # setindex!(Δy, Δyᵢ, i)
          # Δy[i] = Δyᵢ
          Δy = vcat(Δy,Δyᵢ)
          # @show Δy
          # @show tᵢ, yᵢ_obs, Δyᵢ
      end
      L = sum(abs2, Δy[2:end]) #, Σ_sol
    end
    return L
end


# function callback(p,l,tmp_sol,data,model,dens,cons,t)
function callback(p,l,data,model,dens,cons,t,unique_t)
  @show l, p

  # Σ_sol = run_model_2(model,p,data,dens,cons,unique_t)

  # fig = plot(unique_t,Σ_sol')
  # scatter!(fig,t,data,legend=:bottomright)
  # display(fig)

  false
end

function scan_param_space(p,data,model,dens,cons,t,unique_t,lb,ub)
      #create rand mesh
      n = 10000
      ps = rand(n,size(p,1)) .* 10 .- 8
      errors = zeros(n)
      for i in 1:n
          # error, a = loss(ps[i,:],data,model,dens,cons,t)
          error = loss_2(ps[i,:],data,model,dens,cons,t,unique_t)
          errors[i] = error
          # @show error
      end
      min_ind = argmin(errors)
      min_error = minimum(errors)
      @show findmin(errors)
      p_best = ps[min_ind,:]
      @show p_best
      if size(p,1) == 2
          fig = scatter(ps[:,1],ps[:,2],log10.(errors),
          zcolor = log10.(errors),camera=[0,90])
          scatter!([p_best[1]],[p_best[2]],[log10(min_error)],markercolor=:green,markersize=5)
          display(fig)
      end
      return p_best
end

function optimization(model,data,dens,cons,options,t)
# fit function - given model, data, densities, options get best parameters
  unique_t = unique(t)
  loss_in = (p) -> loss_2(p,data,model,dens,cons,t,unique_t)
  # callback_in = (p,l,tmp_sol) -> callback(p,l,tmp_sol,data,model,dens,cons,t)
  callback_in = (p,l) -> callback(p,l,data,model,dens,cons,t,unique_t)


  # @unpack lb,ub = options
  lb,ub = options

  start = scan_param_space(pinit,data,model,dens,cons,t,unique_t,lb,ub) # !! change behavior so that it doesnt need pinit

  res = DiffEqFlux.sciml_train(loss_in,start,BFGS(),cb=callback_in)

  rates = res.minimizer
  l = res.minimum

  # bootstrap_rates = bootstrapping(rates,data,model,dens,cons,t)
  bootstrap_rates = []
  
  f = (du,u,p,t) -> model(du,u,p,t,dens,cons) # optimize in log space
  prob = ODEProblem(f,u₀,tspan,10 .^ rates)
  sol = solve(prob,Vern7())
  h = plot!(sol)
  display(h)
  @show model, rates, start

  hes_for = ForwardDiff.hessian(loss_in,rates)
  hes_zyg = Zygote.hessian(loss_in,rates)
  @show hes_for ≈ hes_zyg

  return rates, hes_zyg, l, bootstrap_rates
end



function bootstrapping(rates,data,model,dens,cons,t)
  #create sample
  data = n
  B_size = 10
  # rand_mat = rand(1:length(data),(length(data),B_size))

  rand_mat = [rand(1:length(data),length(data)) for i in 1:B_size]

  bootstrap_rates = zeros(size(rates))
  # bootstrap_rates = []
  for i in 1:B_size
    resampled_data = data[rand_mat[i]]
    resampled_t = t[rand_mat[i]]
    sort_ind = sortperm(resampled_t)
    resampled_data = resampled_data[sort_ind]
    resampled_t = resampled_t[sort_ind]
    unique_t = unique(resampled_t)

    # loss_in = (p) -> loss_2(p,resampled_data,model,dens,cons,resampled_t,unique_t)
    # callback_in = (p,l) -> callback(p,l,resampled_data,model,dens,cons,resampled_t,unique_t)
    loss_in = (p) -> loss_2(p,data,model,dens,cons,t,unique(t))
    callback_in = (p,l) -> callback(p,l,data,model,dens,cons,t,unique(t))

    # @show loss_in(rates)
    res = optimize(loss_in,-9, 2, rates,Fminbox(NelderMead()))
    # @time res = DiffEqFlux.sciml_train(loss_in,rates,BFGS(),
    #                                 lower_bounds = -9,
    #                                 upper_bounds = 2,
    #                                 cb = callback_in,
    #                                 maxiters=50,
    #                                 f_tol = 1e-4)

    # @show res.minimizer
    bootstrap_rates= hcat(bootstrap_rates,res.minimizer)
    # push!(bootstrap_rates,res.minimizer)
    # Σ_sol = run_model_2(model,p,resampled_data,dens,cons,unique_t)
  end

  return bootstrap_rates[:,2:end]
end
