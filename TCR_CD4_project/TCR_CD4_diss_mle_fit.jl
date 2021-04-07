## modeling of TCR-MHC-CD4 kinetics from adhesion frequency experiments
##
using Revise
src_dir = "C:/Users/stravaglino3/Downloads/julia_code/src/"
# includet(src_dir*"model_fit_global.jl")
includet(src_dir*"model_fit_global_v2.jl")
includet(src_dir*"models.jl")
includet(src_dir*"models_dissociation.jl")
includet(src_dir*"helper_functions.jl")
includet(src_dir*"binning_helper_fs.jl")
includet(src_dir*"global_plotting.jl")
includet(src_dir*"mle_fitting.jl")
using Plots
plot(0:10)
using StatsBase, KernelDensity, Distributions
using BlackBoxOptim
using DataFrames
# load data
data_dir = "C:/Users/stravaglino3/Downloads/julia_code/TCR_CD4_project/dissociation data/"
tm_j = CSV.read(data_dir*"jurkat_TCR_MHC.csv", DataFrame)
cm_j = CSV.read(data_dir*"jurkat_CD4_MHC.csv", DataFrame)
tmc_j = CSV.read(data_dir*"jurkat_TCR_MHC_CD4.csv", DataFrame)

tm_b = CSV.read(data_dir*"Bimolecular polylink interactions.txtbinned bond lifetime.csv", DataFrame)
cm_b = CSV.read(data_dir*"CD4HLA only.txtbinned bond lifetime.csv", DataFrame)
tmc_b = CSV.read(data_dir*"All from 2018-05-1_10_15.txtbinned bond lifetime.csv", DataFrame)

# clean data from 0s
tm_j_binned = clean_dissociation_data(tm_j)
cm_j_binned = clean_dissociation_data(cm_j)
tmc_j_binned = clean_dissociation_data(tmc_j)

tm_b = clean_dissociation_data(tm_b)
cm_b = clean_dissociation_data(cm_b)
tmc_b = clean_dissociation_data(tmc_b)
# cobine bins to complete dataset
F_tm_j, t_tm_j = combine_bins(tm_j_binned)
F_cm_j, t_cm_j = combine_bins(cm_j_binned)
F_tmc_j, t_tmc_j = combine_bins(tmc_j_binned)

F_tm_b, t_tm_b = combine_bins(tm_b)
F_cm_b, t_cm_b = combine_bins(cm_b)
F_tmc_b, t_tmc_b = combine_bins(tmc_b)

gr()
s1 = scatter(F_tm_j, t_tm_j,label="TCR-MHC");
s2 = scatter(F_cm_j, t_cm_j,label="CD4-MHC");
s3 = scatter(F_tmc_j, t_tmc_j,label="TCR-MHC-CD4");
plot(s1,s2,s3); title!("Jurkats")
s1 = scatter(F_tm_b, t_tm_b,label="TCR-MHC");
s2 = scatter(F_cm_b, t_cm_b,label="CD4-MHC");
s3 = scatter(F_tmc_b, t_tmc_b,label="TCR-MHC-CD4");
plot(s1,s2,s3); title!("Beads")

plotly()
plot()
Color=:black; plot_F_vs_t(tm_j)
Color=:blue; plot_F_vs_t(cm_j)
Color=:red; plot_F_vs_t(tmc_j)
title!("Jurkats")
plot()
Color=:black; plot_F_vs_t(tm_b)
Color=:blue; plot_F_vs_t(cm_b)
Color=:red; plot_F_vs_t(tmc_b)
title!("Beads")

##
# load new data from Kaitao and compare to Muaz data
cm_b_new = CSV.read(data_dir*"CD4_LT_bead_KL.csv", DataFrame)
F_cm_b_new = cm_b_new[:,1]; t_cm_b_new = cm_b_new[:,2]
scatter(F_cm_b, log.(t_cm_b),label="Muaz")
scatter!(F_cm_b_new, log.(t_cm_b_new), label="Kaitao")
xlabel!("F"); ylabel!("ln(t)"); title!("CD4-MHC (beads)")

## fit bimolecular TCR-MHC (on Jurkats)

tm_j  = F_t_Data(F_tm_b, t_tm_b, InterpKDE(kde(F_tm_b)))
cm_j   = F_t_Data(F_cm_b, t_cm_b, InterpKDE(kde(F_cm_b)))
tmc_j = F_t_Data(F_tmc_b, t_tmc_b, InterpKDE(kde(F_tmc_b)))

model = bimolecular_diss_bell_activated_state_w_cons!
u₀ = [] # don't provide initial states
u₀_type = "equilibrium" # calculate initial states based on equilibrium
tspan = (0.0,maximum(tm_j.t))
tspan = (0.0,10)
k₋ᵒf_tm = 0.74
M_tm_j = M(model,u₀,u₀_type,tspan,[],[k₋ᵒf_tm],0.0,Rodas5)

pinit = [   0.01,
        0.1,0.6,
        3.,0.2,
        150.0,4.0]
opt =  [        (0.,10.),
        (0.,1.),(0.,10.),
        (0.,10.),(0.,10.),
        (0.,500.),(0.,10.)]


# pinit = rand(8)*10
@time mle_loss(pinit,tm_j,M_tm_j)
using BlackBoxOptim

@time rates_tm_j, l_tm_j, hes_tm_j = optimization(M_tm_j,tm_j,opt)
rates_tm_j =[0.13772078920646594, 0.9999999882398696, 0.2216191952783784, 6.454133898001537e-5, 1.6550629048008287, 499.08524377458565, 9.935988576599074]
ks_tm_j = rates_tm_j
σ_tm_j = hessian2σ(hes_tm_j)
@show k_tm_j = rates_tm_j .± σ_tm_j

includet(src_dir*"mle_plotting.jl")
plotly()
plot_ft_fit(M_tm_j,ks_tm_j,tm_j)

# fitting is nonsensical

# try fitting again, this time allowing k₋ᵒf_tm to change as well
##
model = bimolecular_diss_bell!
u₀ = [] # don't provide initial states
u₀_type = "equilibrium" # calculate initial states based on equilibrium
tspan = (0.0,10)
M_tm_j  = M(model,u₀,u₀_type,tspan,[],[],0.0,Rodas5)

pinit = [14.1,0.01,
        0.1,0.6,
        3.,0.2,
        150.0,4.0]

opt =  [(0.,10.),(0.,10.),
        (0.,1.),(0.,10.),
        (0.,10.),(0.,10.),
        (0.,500.),(0.,10.)]

@time mle_loss(pinit,tm_j,M_tm_j )

@time rates_tm_j, l_tm_j, hes_tm_j = optimization(M_tm_j,tm_j,opt)
ks_tm_j = rates_tm_j
σ_tm_j = hessian2σ(hes_tm_j)
@show k_tm_j = rates_tm_j .± σ_tm_j

includet(src_dir*"mle_plotting.jl")
plotly()
plot_ft_fit(M_tm_j,ks_tm_j,tm_j)
plot(); Color = :blue
plot_mean_ft_fit(M_tm_j,ks_tm_j,tm_j,tm_j_binned;species="T-M")

## try fitting kaitao cd4 data
