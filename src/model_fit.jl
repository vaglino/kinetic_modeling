# module model_fit

using DifferentialEquations, Flux, Optim, DiffEqFlux, DiffEqSensitivity, Plots
using BlackBoxOptim
using BenchmarkTools
using RecursiveArrayTools 
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


# function run_model(model,p,data,densities,constants,t) # version of run function with multiple models
#
#     # Σ_sol_stack = Array{Any}(undef, 1, size(data,2))
#     Σ_sol_stack = zeros(1, size(data,2))
#     # display(1:size(data,1))
#     for i in 1:size(data,1)
#         # run model with given densities
#         densities_i = densities[i,:]
#         f = (du,u,p,t) -> model(du,u,10 .^ p,t,densities_i,constants)
#         tmp_prob = ODEProblem(f,u₀,tspan,p)
#         # tmp_sol = solve(tmp_prob,Rosenbrock23(),saveat=t)
#         tmp_sol = solve(tmp_prob,Rodas4(),saveat=t)
#         # stack temp solution
#         # stack Σ of solution across n species
#         Σ_sol = sum(Array(tmp_sol),dims=1)
#         Σ_sol_stack = vcat(Σ_sol_stack,Σ_sol)
#     end
#     return Σ_sol_stack[2:end,:]
# end

function run_model(model,p,data,densities,constants,t) # version of run function with multiple models

    # Σ_sol_stack = Array{Any}(undef, 1, size(data,2))
    # Σ_sol_stack = zeros(1, size(data,2))


    f = (du,u,p,t) -> model(du,u,10 .^ p,t,densities,constants)
    tmp_prob = ODEProblem(f,u₀,tspan,p)
    # tmp_sol = solve(tmp_prob,Rosenbrock23(),saveat=t)
    # tmp_sol = solve(tmp_prob,Rodas4(),saveat=t)
    tmp_sol = solve(tmp_prob,Tsit5(),saveat=t)
    # stack temp solution
    # stack Σ of solution across n species
    Σ_sol = sum(Array(tmp_sol),dims=1)

end

function loss(p,data,model,dens,cons,t)
    Σ_sol = run_model(model,p,data,dens,cons,t)
    # display(Σ_sol)
    # sum(abs2, (Σ_sol - data) .* t), Σ_sol
    sum(abs2, (Σ_sol - data))#, Σ_sol
end

function callback(p,l,tmp_sol,data,model,dens,cons,t)
    @show l

    Σ_sol = run_model(model,p,data,dens,cons,t)

    fig = plot(t,Σ_sol')
    scatter!(fig,t,data',legend=:bottomright)
    display(fig)

    false
end

function scan_param_space(p,data,model,dens,cons,t,lb,ub)
        #create rand mesh
        n = 10000
        ps = rand(n,size(p,1)) .* 9 .- 7
        errors = zeros(n)
        for i in 1:n
            # error, a = loss(ps[i,:],data,model,dens,cons,t)
            error = loss(ps[i,:],data,model,dens,cons,t)
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

    loss_in = (p) -> loss(p,data,model,dens,cons,t)
    callback_in = (p,l,tmp_sol) -> callback(p,l,tmp_sol,data,model,dens,cons,t)

    # @unpack lb,ub = options
    lb,ub = options
    # f = (du,u,p,t) -> model(du,u,p,t,dens,cons) # optimize in log space
    # prob = ODEProblem(f,u₀,tspan,p)

    # start = scan_param_space(pinit,data,model,dens,cons,t,lb,ub) # !! change behavior so that it doesnt need pinit
    # use ADAM first to find ballparck solution
    # res = DiffEqFlux.sciml_train(loss_in,pinit,ADAM(0.001),cb=callback_in,maxiters=300)
    # res = DiffEqFlux.sciml_train(loss_in,pinit,BFGS(initial_stepnorm = 0.0001),cb=callback_in)
    start = [-4.,-1.]
    res = DiffEqFlux.sciml_train(loss_in,start,BFGS(),cb=callback_in)
    # res = DiffEqFlux.sciml_train(loss_in,start,NelderMead(),cb=callback_in)
    @show res
    rates = res.minimizer
    # rates = start

    f = (du,u,p,t) -> model(du,u,p,t,dens,cons) # optimize in log space
    prob = ODEProblem(f,u₀,tspan,10 .^ rates)
    sol = solve(prob,Rodas4())
    h = plot!(sol)
    display(h)
    display([model,rates,start,dens,cons])
    return rates
end

# end
