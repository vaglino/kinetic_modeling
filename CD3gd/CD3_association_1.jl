using Revise
include("model_fit_global.jl")
include("models.jl")
include("helper_functions.jl")


#=  Densities from Zhou
Micropipette	TCR	         CD3de  	CD3ge   	CD3mix
RBC coating	   165.6213223	 406.7805	742.7222634	550.7926
site density   198.7162679	 356.0556	576.7556	415.3837
=#
Ac = 1.0
densities_γT = [576.7556, 198.7162679, Ac]
densities_δT = [356.0556, 198.7162679, Ac]
densities_mixT = [415.3837/2, 198.7162679, 415.3837/2, Ac]
densities_γδ = [576.7556, 356.0556, Ac]
#= bimolecular rates from Zhou
        on rate   off rate    =#
rates_γT = [3.6e-5, 14.1]
rates_δT = [2.3e-5, 11.1]
rates_γδ = [1.6e-4, 14.9]


function solve_model(rates_sds)
    prob = ODEProblem(f,u₀,tspan,rates_sds)
    # prob = remake(prob, p=10 .^ rates)
    sol = solve(prob,Vern7(),saveat=0.01, abs_tol=1e-10,rel_tol=1e-10)
    return sol
end


tspan = (0.0,10.0)
t = collect(0.0:0.1:10)

gr()
# bimolecular TCR-pMHC bonds
u₀ = [0.0]

f = (du,u,p,t) -> bimolecular!(du,u,p,t,densities_γT,[])
γT = solve_model(rates_γT)
plot(γT)

f = (du,u,p,t) -> bimolecular!(du,u,p,t,densities_δT,[])
δT = solve_model(rates_δT)
plot!(δT)

f = (du,u,p,t) -> bimolecular!(du,u,p,t,densities_γδ,[])
γδ = solve_model(rates_γδ)
plot!(γδ)
