## modeling of TCR-MHC-CD4 kinetics from adhesion frequency experiments
##
using Revise
src_dir = "C:/Users/stravaglino3/Downloads/julia_code/src/"
# includet(src_dir*"model_fit_global.jl")
includet(src_dir*"model_fit_global_v2.jl")
includet(src_dir*"models.jl")
includet(src_dir*"helper_functions.jl")
includet(src_dir*"global_plotting.jl")
using Plots
plot(0:10)

# experimental conditions
##
# Site densities (from flow cytometry)
densities_T =   [25.0, 38.0, 1.0]       # TCR, pMHC, Ac
densities_TC1 = [15.0, 38.0, 10.0, 1.0] # TCR, pMHC, CD4, Ac
densities_TC2 = [8.0,  38.0, 17.0, 1.0]
densities_TC3 = [3.0,  38.0, 22.0, 1.0]
densities_C =   [1700., 850., 1.0]       # CD4, pMHC, A

# Importing experimental data
##
using CSV
data_dir = "C:/Users/stravaglino3/Downloads/julia_code/TCR_CD4_project/data/"
tcr_data = CSV.read(data_dir*"tcr_only_data.csv", DataFrame)
tcr_cd4_1 = CSV.read(data_dir*"tcr_cd4_1.csv", DataFrame)
tcr_cd4_2 = CSV.read(data_dir*"tcr_cd4_2.csv", DataFrame)
tcr_cd4_3 = CSV.read(data_dir*"tcr_cd4_3.csv", DataFrame)
cd4_data = CSV.read(data_dir*"cd4_data_from_muaz.csv", DataFrame)
t_t,   Pa_t,   n_se_t   = clean_data(tcr_data)
t_tc1, Pa_tc1, n_se_tc1 = clean_data(tcr_cd4_1)
t_tc2, Pa_tc2, n_se_tc2 = clean_data(tcr_cd4_2)
t_tc3, Pa_tc3, n_se_tc3 = clean_data(tcr_cd4_3)
t_c,   Pa_c,   n_se_c   = clean_data(cd4_data)

n_t   = Pa2n(Pa_t)
n_tc1 = Pa2n(Pa_tc1)
n_tc2 = Pa2n(Pa_tc2)
n_tc3 = Pa2n(Pa_tc3)
n_c   = Pa2n(Pa_c)

scatter(n_se_t[:,1],n_se_t[:,2],yerror=n_se_t[:,3],markercolor=:black,markersize=5)
scatter!(n_se_tc1[:,1],n_se_tc1[:,2],yerror=n_se_tc1[:,3],markercolor=:blue,markersize=5)
scatter!(n_se_tc2[:,1],n_se_tc2[:,2],yerror=n_se_tc2[:,3],markercolor=:red,markersize=5)
scatter!(n_se_tc3[:,1],n_se_tc3[:,2],yerror=n_se_tc3[:,3],markercolor=:cyan,markersize=5,legend=false)
# scatter!(n_se_c[:,1],n_se_c[:,2],yerror=n_se_c[:,3],markercolor=:orange,markersize=5)

# TCR-MHC fit (bimolecular)
##
using BlackBoxOptim, Measurements

tspan = (0,maximum(t_t)+1)
u₀ = [0.0]
bounds = (1.0e-6,100.0)
pinit = [-4.0,0.0]

tcr = Data([n_t],[t_t])
M_tcr = M_global(bimolecular!,u₀,tspan,[densities_T],[[]],0.0,Rodas5)

rates_tcr,hes_tcr,l_tcr = optimization(M_tcr,tcr,bounds)
σ_tcr = hessian2σ(hes_tcr,l_tcr,tcr.n[1],rates_tcr)
@show k_t = rates_tcr .± σ_tcr

f = (du,u,p,t) -> bimolecular!(du,u,p,t,densities_T, [])
Color = :black; plot_sum(rates_tcr,:black)

# CD4-MHC fit (bimolecular)
##
tspan = (0,maximum(t_c)+1)
u₀ = [0.0]
bounds = (1.0e-9,100.0)
pinit = [-4.0,0.0]

cd4 = Data([n_c],[t_c])
M_cd4 = M_global(bimolecular!,u₀,tspan,[densities_C],[[]],0.0,Rodas5)

rates_cd4,hes_cd4,l_cd4 = optimization(M_cd4,cd4,bounds)
σ_cd4 = hessian2σ(hes_cd4,l_cd4,cd4.n[1],rates_cd4)
@show k_c = rates_cd4 .± σ_cd4

# TCR-MHC-CD4 fit (1 path)
##
u₀ = [0.0, 0.0]
cons_tcr = rates_tcr
cons_tc = vcat(cons_tcr, rates_cd4)
tspan = (0,maximum(t_tc1)+1)
pinit = [0.0,0.0]

tc1 = Data([n_tc1],[t_tc1])
tc2 = Data([n_tc2],[t_tc2])
tc3 = Data([n_tc3],[t_tc3])
M_tc1 = M_global(trimolecular!,u₀,tspan,[densities_TC1],[cons_tcr],0.0,Rodas5)
M_tc2 = M_global(trimolecular!,u₀,tspan,[densities_TC2],[cons_tcr],0.0,Rodas5)
M_tc3 = M_global(trimolecular!,u₀,tspan,[densities_TC3],[cons_tcr],0.0,Rodas5)


rates_tc1,hes_tc1,l_tc1 = optimization(M_tc1,tc1,bounds)
σ_tc1 = hessian2σ(hes_tc1,l_tc1,tc1.n[1],rates_tc1)
rates_tc2,hes_tc2,l_tc2 = optimization(M_tc2,tc2,bounds)
σ_tc2 = hessian2σ(hes_tc2,l_tc2,tc2.n[1],rates_tc2)
rates_tc3,hes_tc3,l_tc3 = optimization(M_tc3,tc3,(1e-6,1e3))
σ_tc3 = hessian2σ(hes_tc3,l_tc3,tc3.n[1],rates_tc3)

@show k_tc1 = rates_tc1 .± σ_tc1
@show k_tc2 = rates_tc2 .± σ_tc2
@show k_tc3 = rates_tc3 .± σ_tc3

# plot 1 path trimolecular fits
##
f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_TC1, rates_tcr)
Color = :blue; plot_sum(rates_tc1,:blue)
# plot_sum([0., 0.])
f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_TC2, rates_tcr)
Color = :red; plot_sum(rates_tc2,:red)
# plot_sum([0., 0.])
f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_TC3, rates_tcr)
Color = :cyan; plot_sum(rates_tc3,:cyan)

# plot_sum([0., 0.])
#
# TCR-MHC-CD4 fit (1 path, with diffusion?)
##
# u₀ = [0., 0., 0.]
# f = (du,u,p,t) -> trimolecular_1path_diff!(du,u,p,t,densities_TC1, rates_tcr)
# plot_sum([1e-3,1e-3,1,1e-3])
# pinit = [0.0,0.0,0.0,0.0]
# M_tc1_diff = M_global(trimolecular_1path_diff!,u₀,tspan,[densities_TC1],[cons_tcr],0.0,Rodas5)
# rates_tc1_diff,hes_tc1_diff,l_tc1_diff = optimization(M_tc1_diff,tc1,(1e-8,1e3))
# f = (du,u,p,t) -> trimolecular_1path_diff!(du,u,p,t,densities_TC1, rates_tcr)
# plot_sum(rates_tc1_diff)
