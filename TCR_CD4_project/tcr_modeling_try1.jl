include("C:/Users/stefano/Documents/Zhulab/julia_code/model_fit.jl")
include("C:/Users/stefano/Documents/Zhulab/julia_code/models.jl")
include("C:/Users/stefano/Documents/Zhulab/julia_code/helper_functions.jl")
using DifferentialEquations, Flux, Optim, DiffEqFlux, DiffEqSensitivity, Plots
using BlackBoxOptim
using BenchmarkTools
using RecursiveArrayTools
# --------------------------create dummy data-----------------------------------
f = (du,u,p,t) -> bimolecular!(du,u,p,t,dens,cons) # enclose constants

# u0 = [0.0;0.0]
u₀ = [0.0]          # Initial condition
tspan = (0.0,10.0)  # Simulation interval
p = [1e-4,0.5]      # equation parameter. p = [k₁, k₋₁]
dens = [100., 10., 1.]
cons = []

# Setup the ODE problem, then solve
prob = ODEProblem(f,u₀,tspan,p)
sol = solve(prob,Tsit5())
# sol = solve(prob,Tsit5())
plot(sol)
dataset1 = generate_data(sol)


model_opt2 = (1.0e-6,10.0)
# pinit = [2.0e-4,0.6]
pinit = [-4.0,0.0]
t = collect(range(0,stop=10,length=101))
rates = optimization(bimolecular!,dataset1,dens',[],model_opt2,t)

# -------------------------------------------------------------------------------


# u0 = [0.0;0.0]
u₀ = [0.0]          # Initial condition
tspan = (0.0,10.0)  # Simulation interval
p = [1e-4,0.5]      # equation parameter. p = [k₁, k₋₁]
dens_new = [100., 50., 1.]
cons = []

f = (du,u,p,t) -> bimolecular!(du,u,p,t,dens_new,cons) # enclose constants
# Setup the ODE problem, then solve
prob = ODEProblem(f,u₀,tspan,p)
sol = solve(prob,Tsit5())
# sol = solve(prob,Tsit5())
plot!(sol)
dataset3 = generate_data(sol)
dataset_comb = vcat(dataset1,dataset3)

dens_comb = reduce(vcat, [dens', dens_new'])
rates = optimization(bimolecular!,dataset_comb,dens_comb,cons,model_opt2,t)


#----------------------------------------------------------------------------
# try TCR-pMHC only data from Muaz

using CSV
tcr_data = CSV.read("C:/Users/stefano/Documents/Zhulab/julia_code/tcr_only_data.csv")

t, Pa, n_se_tcr = clean_data(tcr_data)
n = Pa2n(Pa)

scatter(t,Pa)
scatter!(t,n)
scatter!(n_se_tcr[:,1],n_se_tcr[:,2],yerror=n_se_tcr[:,3],markercolor=:black)

densities_tcr = [25.0, 38.0, 1.0]
tspan = (0,maximum(t)+1)
u₀ = [0.0]
f = (du,u,p,t) -> bimolecular!(du,u,p,t,densities_tcr,[]) # enclose constants

pinit = [-4.5,-0.7]

rates_tcr = optimization(bimolecular!,n',densities_tcr',[],model_opt2,t)

prob = ODEProblem(f,u₀,tspan,10 .^ rates_tcr)
# prob = remake(prob, p=10 .^ rates)
sol = solve(prob,Rosenbrock23())
plot!(sol,linecolor=:black)
summ = reduce(vcat,sum(sol,dims=1))
h = plot!(sol.t,summ,linecolor=:black)

# prob = ODEProblem(f,u₀,tspan,[1e-4,0.5])

#------------------------------------------------------------------------------
# try TCR-pMHC-CD4 fit with 2-path, 2-step

tcr_cd4_1 = CSV.read("C:/Users/stefano/Documents/Zhulab/julia_code/tcr_cd4_1.csv")
tcr_cd4_2 = CSV.read("C:/Users/stefano/Documents/Zhulab/julia_code/tcr_cd4_2.csv")
tcr_cd4_3 = CSV.read("C:/Users/stefano/Documents/Zhulab/julia_code/tcr_cd4_3.csv")
t_tc1, Pa_tc1, n_se_1 = clean_data(tcr_cd4_1)
t_tc2, Pa_tc2,  n_se_2  = clean_data(tcr_cd4_2)
t_tc3, Pa_tc3,  n_se_3 = clean_data(tcr_cd4_3)

n_tc1 = Pa2n(Pa_tc1)
n_tc2 = Pa2n(Pa_tc2)
n_tc3 = Pa2n(Pa_tc3)

# scatter(t,Pa)
scatter(t_tc1,n_tc1)
scatter!(t_tc2,n_tc2)
scatter!(t_tc3,n_tc3)


densities_tc1 = [15.0, 38.0, 10.0, 1.0] # TCR, pMHC, CD4, Ac
tspan = (0,maximum(t)+1)
u₀ = [0.0, 0.0]

cons_tcr = 10 .^ rates_tcr
# f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_tc1,cons_tcr)
# rates_tc1 = optimization(trimolecular!,n_tc1',densities_tc1',cons_tcr,model_opt2,t_tc1)


u₀ = [0.0, 0.0, 0.0]
# CD4 rates from Muaz, need to get full CD4 curve
kon2 = 7.01979340915551e-06 / 6
koff2 = 0.4000
cons_tc = vcat(cons_tcr, kon2, koff2)

# cons_tcr_from_mat = [1.4e-4, 0.5]
# cons_tc = vcat(cons_tcr_from_mat, kon2, koff2)

f = (du,u,p,t) -> two_path!(du,u,p,t,densities_tc1,cons_tc) # enclose constants
# p = [9.32e-7,1.81e-1,1.7e-2,5.23e-5]
p = [0.0, 0.0, 0.0, 0.0]
prob = ODEProblem(f,u₀,tspan,p)
# prob = remake(prob, p=10 .^ rates)
sol = solve(prob,Rosenbrock23())
summ = reduce(vcat,sum(sol,dims=1))
plot!(sol.t,summ)

pinit = [1e-2,0.5,1e-2,0.5]
rates_tc1_2paths = optimization(two_path!,n_tc1',densities_tc1',cons_tc,model_opt2,t_tc1)

densities_tc2 = [8.0, 38.0, 17.0, 1.0]
rates_tc2_2paths = optimization(two_path!,n_tc2',densities_tc2',cons_tc,model_opt2,t_tc2)

densities_tc3 = [3.0, 38.0, 22.0, 1.0]
rates_tc3_2paths = optimization(two_path!,n_tc3',densities_tc3',cons_tc,model_opt2,t_tc3)


scatter(n_se_1[:,1],n_se_1[:,2],yerror=n_se_1[:,3],markercolor=:red)
scatter!(n_se_2[:,1],n_se_2[:,2],yerror=n_se_2[:,3],markercolor=:blue)
scatter!(n_se_3[:,1],n_se_3[:,2],yerror=n_se_3[:,3],markercolor=:cyan)

f = (du,u,p,t) -> two_path!(du,u,p,t,densities_tc1,cons_tc)
prob = ODEProblem(f,u₀,tspan,10 .^ rates_tc1_2paths)
sol = solve(prob,Rosenbrock23())
summ1 = reduce(vcat,sum(sol,dims=1))
time = sol.t
plot!(time,summ1,linecolor=:red)

f = (du,u,p,t) -> two_path!(du,u,p,t,densities_tc2,cons_tc)
prob = ODEProblem(f,u₀,tspan,10 .^ rates_tc2_2paths)
sol = solve(prob,Rosenbrock23())
summ = reduce(vcat,sum(sol,dims=1))
time = sol.t
plot!(sol.t,summ,linecolor=:blue)

f = (du,u,p,t) -> two_path!(du,u,p,t,densities_tc3,cons_tc)
prob = ODEProblem(f,u₀,tspan,10 .^ rates_tc3_2paths)
sol = solve(prob,Rosenbrock23())
summ = reduce(vcat,sum(sol,dims=1))
plot!(sol.t,summ,linecolor=:cyan)

function plot_sum(rates)
    prob = ODEProblem(f,u₀,tspan,rates)
    # prob = remake(prob, p=10 .^ rates)
    sol = solve(prob,Rosenbrock23())
    summ = reduce(vcat,sum(sol,dims=1))
    plot!(sol.t,summ)
end


# dprob = DiscreteProblem(u₀, tspan, 10 .^ rates)
#
# # now we create a JumpProblem, and specify Gillespie's Direct Method as the solver:\n",
# jprob = JumpProblem(dprob, Direct(), f)
#
# # now let's solve and plot the jump process:
# sol = solve(jprob, SSAStepper())
# plot(sol, fmt=:svg, ylims = (0, 2))
#
# prob = DiscreteProblem(bimol, u₀, tspan, 10 .^ rates)
# jump_prob = JumpProblem(bimol, prob, Direct())
# sol = solve(jump_prob, SSAStepper())
#
# plot(sol)
#
# function bimol(du,u,p,t)
#     # unpack rates and constants
#     nᵣ = u[1]
#     k₁,k₋₁  = p
#     mᵣ,mₗ,A = [100,10,1]
#     # model
#     du[1] = dnᵣ = A*k₁*mᵣ*mₗ - k₋₁*nᵣ
#
# end
