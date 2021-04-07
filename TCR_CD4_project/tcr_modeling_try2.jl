using Revise
src_dir = "C:/Users/stravaglino3/Downloads/julia_code/src/"
# includet(src_dir*"model_fit_global.jl")
includet(src_dir*"model_fit_global_v2.jl")
includet(src_dir*"models.jl")
includet(src_dir*"helper_functions.jl")
includet(src_dir*"global_plotting.jl")
using Plots
plot(0:10)
# --------------------------create dummy data-----------------------------------
# f = (du,u,p,t) -> bimolecular!(du,u,p,t,dens,cons) # enclose constants
#
# # u0 = [0.0;0.0]
# u₀ = [0.0]          # Initial condition
# tspan = (0.0,10.0)  # Simulation interval
# p = [1e-4,0.5]      # equation parameter. p = [k₁, k₋₁]
# dens = [100., 10., 1.]
# cons = []
#
# # Setup the ODE problem, then solve
# prob = ODEProblem(f,u₀,tspan,p)
# sol = solve(prob,Tsit5())
# # sol = solve(prob,Tsit5())
# plot(sol)
# dataset1 = generate_data(sol)
#
# dens = [50., 35., 1.]
# prob = ODEProblem(f,u₀,tspan,p)
# sol = solve(prob,Tsit5())
# # sol = solve(prob,Tsit5())
# plot!(sol)
# dataset2 = generate_data(sol)
#
# model_opt2 = (1.0e-6,10.0)
# # pinit = [2.0e-4,0.6]
# pinit = [-4.0,0.0]
# t = collect(range(0,stop=10,length=101))
# dens_comb = [[100., 10., 1.],[50., 35., 1.]]
# rates,hes,l = optimization(bimolecular!,[vec(dataset1),vec(dataset2)],
#                             dens_comb,[[],[]],model_opt2,[t,t])


#----------------------------------------------------------------------------
# try TCR-pMHC only data from Muaz

using CSV

data_dir = "C:/Users/stravaglino3/Downloads/julia_code/TCR_CD4_project/data/"
tcr_data = CSV.read(data_dir*"tcr_only_data.csv", DataFrame)

@time t, Pa, n_se_tcr = clean_data(tcr_data)
n = Pa2n(Pa)

# scatter(t,Pa)
# using Plots
scatter(t,n)
scatter!(n_se_tcr[:,1],n_se_tcr[:,2],yerror=n_se_tcr[:,3],markercolor=:black)

densities_tcr = [[25.0, 38.0, 1.0]]
tspan = (0,maximum(t)+1)
u₀ = [0.0]
bounds = (1.0e-6,100.0)
pinit = [-4.0,0.0]
f = (du,u,p,t) -> bimolecular!(du,u,p,t,densities_tcr,[]) # enclose constants

pinit = [-4.5,-0.7]

tcr = Data([n],[t])
M_tcr = M(bimolecular!,u₀,tspan,densities_tcr,[[]],0.0,Rodas5)

# rates_tcr,hes_tcr,l_tcr = optimization(bimolecular!,[n],[densities_tcr],[[]],model_opt2,[t])
using BlackBoxOptim
rates_tcr,hes_tcr,l_tcr = optimization(M_tcr,tcr,bounds)
# σ_tcr = hessian2σ(hes_tcr)
σ_tcr = hessian2σ(hes_tcr,l_tcr,tcr.n[1],rates_tcr)
# u₀ = [0.0,0.]
# rates_tcr,hes_tcr,l_tcr = optimization(trimolecular!,n,densities_tc1',rates_tc,model_opt2,t)
f = (du,u,p,t) -> bimolecular!(du,u,p,t,densities_tcr[1],[])
plot_sum(rates_tcr)

# unique_t = [unique(times) for times in tcr.t]
# loss_in = (p) -> loss(p,tcr,M_tcr,unique_t)
# hes_for =ForwardDiff.hessian(loss_in,rates_tcr)
# σ_tcr = hessian2σ(hes_for)
# plot_curve_sds(hes_tcr,l_tcr,n,rates_tcr)
# calculate_sds(hes_tcr,l_tcr,n,rates_tcr)
# prob = ODEProblem(f,u₀,tspan,10 .^ rates_tcr)
# # prob = remake(prob, p=10 .^ rates)
# sol = solve(prob,Rodas4())
# plot!(sol)
# # fig = plot!(sol,linecolor=:black)
# summ = reduce(vcat,sum(sol,dims=1))
# h = plot!(sol.t,summ,linecolor=:black)

# prob = ODEProblem(f,u₀,tspan,[1e-4,0.5])

#------------------------------------------------------------------------------
# try TCR-pMHC-CD4 fit with 2-path, 2-step

tcr_cd4_1 = CSV.read(data_dir*"tcr_cd4_1.csv", DataFrame)
tcr_cd4_2 = CSV.read(data_dir*"tcr_cd4_2.csv", DataFrame)
tcr_cd4_3 = CSV.read(data_dir*"tcr_cd4_3.csv", DataFrame)
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
densities_tc2 = [8.0,  38.0, 17.0, 1.0]
densities_tc3 = [3.0,  38.0, 22.0, 1.0]


u₀ = [0.0, 0.0]
# cons_tcr = 10 .^ rates_tcr
cons_tcr = rates_tcr

# CD4 rates from Muaz, need to get full CD4 curve
kon2 = 7.01979340915551e-06 / 6
koff2 = 0.4000
ka = 8.09e-7 / 6
kon2 = ka*0.4000
koff2 = 0.4000
cons_tc = vcat(cons_tcr, kon2, koff2)
tspan = (0,maximum(t)+1)

pinit = [0.0,0.0]

tc1 = Data([n_tc1],[t_tc1])
tc2 = Data([n_tc2],[t_tc2])
tc3 = Data([n_tc3],[t_tc3])
M_tc1 = M(trimolecular!,u₀,tspan,[densities_tc1],[cons_tc],0.0,Rodas5)
M_tc2 = M(trimolecular!,u₀,tspan,[densities_tc2],[cons_tc],0.0,Rodas5)
M_tc3 = M(trimolecular!,u₀,tspan,[densities_tc3],[cons_tc],0.0,Rodas5)


rates_tc1,hes_tc1,l_tc1 = optimization(M_tc1,tc1,bounds)
σ_tc1 = hessian2σ(hes_tc1,l_tc1,tc1.n[1],rates_tc1)
rates_tc2,hes_tc2,l_tc2 = optimization(M_tc2,tc2,bounds)
σ_tc2 = hessian2σ(hes_tc2,l_tc2,tc2.n[1],rates_tc2)
rates_tc3,hes_tc3,l_tc3 = optimization(M_tc3,tc3,(1e-6,1e3))
σ_tc3 = hessian2σ(hes_tc3,l_tc3,tc3.n[1],rates_tc3)

# rates_tc1,hes_tc1,l_tc1,boot_tc1 = optimization(trimolecular!,n_tc1,densities_tc1',
#                                     10 .^ rates_tcr,model_opt2,t_tc1)

# boot = bootstrapping(rates_tc1,n_tc1,trimolecular!,densities_tc1,10 .^rates_tcr,t_tc1)

f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_tc1, rates_tcr)
plot_sum(rates_tc1)
plot_sum([0., 0.])
f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_tc2, rates_tcr)
plot_sum(rates_tc2)
plot_sum([0., 0.])
f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_tc3, rates_tcr)
plot_sum(rates_tc3)
plot_sum([0., 0.])


display_rates([rates_tc1,rates_tc2,rates_tc3])
# plot_curve_sds(hes_tc1,l_tc1,n_tc1,rates_tc1)
# rates_tc1_sds = calculate_sds(hes_tc1,l_tc1,n_tc1,rates_tc1)

# rates_tc1,hes_tc2,l_tc2  = optimization(trimolecular!,n_tc2,densities_tc2',cons_tc,model_opt2,t_tc2)
# rates_tc1,hes_tc2,l_tc2  = optimization(trimolecular!,n_tc3,densities_tc3',cons_tc,model_opt2,t_tc3)


scatter(n_se_1[:,1],n_se_1[:,2],yerror=n_se_1[:,3],markercolor=:red)
scatter!(n_se_2[:,1],n_se_2[:,2],yerror=n_se_2[:,3],markercolor=:blue)
scatter!(n_se_3[:,1],n_se_3[:,2],yerror=n_se_3[:,3],markercolor=:cyan)

f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_tc1,10 .^rates_tcr)
prob = ODEProblem(f,u₀,tspan,10 .^ rates_tc1_sds)
sol = solve(prob,Vern7(),saveat=1)
summ1 = reduce(vcat,sum(sol,dims=1))
time = sol.t
plot!(summ1,linecolor=:red)

f = (du,u,p,t) -> two_path!(du,u,p,t,densities_tc2,cons_tc)
prob = ODEProblem(f,u₀,tspan,10 .^ rates_tc2_2paths)
sol = solve(prob,Rodas4())
summ = reduce(vcat,sum(sol,dims=1))
time = sol.t
plot!(sol.t,summ,linecolor=:blue)

f = (du,u,p,t) -> two_path!(du,u,p,t,densities_tc3,cons_tc)
prob = ODEProblem(f,u₀,tspan,10 .^ rates_tc3_2paths)
sol = solve(prob,Rodas4())
summ = reduce(vcat,sum(sol,dims=1))
plot!(sol.t,summ,linecolor=:cyan)

function plot_sum(rates)
    prob = ODEProblem(f,u₀,tspan,rates)
    # prob = remake(prob, p=10 .^ rates)
    sol = solve(prob,Rosenbrock23())
    summ = reduce(vcat,sum(sol,dims=1))
    plot!(sol.t,summ)
end

function display_rates(rates)
    for rate in rates
        display(10 .^ rate')
    end
end

display_rates([rates_tc1_2paths,rates_tc2_2paths,rates_tc3_2paths])


# function find_stds(model,rates,data)
# end

# rates = optimization(bimolecular!,dataset1,dens',[],model_opt2,t)
# loss(p,data,model,dens,cons,t)
using ForwardDiff
using Zygote
using ReverseDiff
# loss_in = (p) -> loss(p,data,model,dens,cons,t)
loss_in = (p) -> loss_2(p,dataset1,bimolecular!,dens',[],t)
loss(p,dataset1,bimolecular!,dens',[],t)
grad_zyg = Zygote.gradient(loss_in,[-4,-0.5])
grad_for = ForwardDiff.gradient(loss_in,[-4,-0.5])
grad_rev = ReverseDiff.gradient(loss_in,[-4,-0.5])
hes_for = ForwardDiff.hessian(loss_in,[-4,-0.297])
hes_rev = ReverseDiff.hessian(loss_in,[-4,-0.297])
grad_zyg = Zygote.gradient(loss_in,rates)
grad_for = ForwardDiff.gradient(loss_in,rates)
grad_rev = ReverseDiff.gradient(loss_in,rates)
hes_for = ForwardDiff.hessian(loss_in,rates)
hes_rev = ReverseDiff.hessian(loss_in,rates)
L = loss_in(rates)
covrp = inv(hes_rev)


m = size(dataset1,2)
n = size(rates,1)
ν = m - n
MSE = L/ν
stds = sqrt.(diag(covrp).*MSE)
using LinearAlgebra
# stds = sqrt.(diag(covrp).*L)

using Measurements
rates_std = rates .± stds

prob = ODEProblem(f,u₀,tspan,10 .^ rates_std)
# prob = remake(prob, p=10 .^ rates)
sol = solve(prob,Tsit5())
plot(sol,errorbars = :ribbon, alpha=0.2)
# fig = plot!(sol,linecolor=:black)


loss_in = (p) -> loss(p,n',bimolecular!,densities_tcr',[],t)
loss(rates_tcr,n',bimolecular!,densities_tcr',[],t)

grad_zyg = Zygote.gradient(loss_in,rates_tcr)
grad_for = ForwardDiff.gradient(loss_in,rates_tcr)
grad_rev = ReverseDiff.gradient(loss_in,rates_tcr)
hes_for = ForwardDiff.hessian(loss_in,rates_tcr)
hes_rev = ReverseDiff.hessian(loss_in,rates_tcr)
hes_zyg = Zygote.hessian(loss_in,rates_tcr)

L = loss_in(rates_tcr)
covrp = inv(hes_for)
using LinearAlgebra
function plot_curve_sds(hes,l,data,rates)
    covrp = inv(hes)
    m = size(data,1)
    n = size(rates,1)
    ν = m - n
    MSE = l/ν
    stds = sqrt.(diag(covrp).*MSE)
    rates_std = rates .± stds

    prob = ODEProblem(f,u₀,tspan,10 .^ rates_std)
    # prob = remake(prob, p=10 .^ rates)
    sol = solve(prob,Vern7(),abs_tol=1e-10,rel_tol=1e-10)
    summ = reduce(vcat,sum(sol,dims=1))
    # time = sol.t
    plot!(sol,errorbars = :ribbon, alpha=0.2)
    # plot!(time,summ,errorbars = :ribbon, alpha=0.2)
    # plot!(time,summ)
end

function calculate_sds(hes,l,data,rates)
    covrp = inv(hes)
    m = size(data,1)
    n = size(rates,1)
    ν = m - n
    MSE = l/ν
    stds = sqrt.(diag(covrp).*MSE)
    rates_std = rates .± stds
end
