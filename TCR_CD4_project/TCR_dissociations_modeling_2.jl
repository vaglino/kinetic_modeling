include("C:/Users/stravaglino3/Downloads/julia_code/julia_code/model_fit_global.jl")
include("C:/Users/stravaglino3/Downloads/julia_code/julia_code/models.jl")
include("C:/Users/stravaglino3/Downloads/julia_code/julia_code/helper_functions.jl")
using DataFrames
tcr_mhc_j = CSV.read("C:/Users/stravaglino3/Downloads/julia_code/dissociation data/jurkat_TCR_MHC.csv", DataFrame)
cd4_mhc_j = CSV.read("C:/Users/stravaglino3/Downloads/julia_code/dissociation data/jurkat_CD4_MHC.csv", DataFrame)
tcr_mhc_cd4_j = CSV.read("C:/Users/stravaglino3/Downloads/julia_code/dissociation data/jurkat_TCR_MHC_CD4.csv", DataFrame)

tcr_mhc = clean_dissociation_data(tcr_mhc_j)
cd4_mhc = clean_dissociation_data(cd4_mhc_j)
tcr_mhc_cd4 = clean_dissociation_data(tcr_mhc_cd4_j)

# pyplot()
plot()
plot_dissociations(tcr_mhc)
plot_dissociations(cd4_mhc)
plot_dissociations(tcr_mhc_cd4)

t_tm, lnp_tm = separate_diss_data(tcr_mhc)
t_cm, lnp_cm = separate_diss_data(cd4_mhc)
t_tmc, lnp_tmc = separate_diss_data(tcr_mhc_cd4)
p_tm = lnp2p.(lnp_tm)
p_cm = lnp2p.(lnp_cm)
p_tmc = lnp2p.(lnp_tmc)

plot()
fs_tm = extract_avg_forces(tcr_mhc)
fs_cm = extract_avg_forces(cd4_mhc)
fs_tmc = extract_avg_forces(tcr_mhc_cd4)

scatter(t_tm,lnp_tm)


## TCR-MHC
u₀ = [1.0, 0.]          # Initial condition
tspan = (0.0,10.0)  # Simulation interval
model_opt = (-5.,3.)
model_opt = (0.,100.)
pinit =[11.,0.01,0.14,0.4,3,0.2,20,4]
# pinit = log10.(pinit)
dens_j = fill([],length(t_tm))
model = bimolecular_diss_bell!
# pinit = 10 .^ [0.8934057970002582, 0.6126600385641812, -2.0499405365827617, 0.7136091557046687, -0.6215834038143795, 0.23515760918441408, -0.11706056162578093, -1.5512600092183506]
rates_tm,hes_tm,l_tm,boot_tm = optimization(model,p_tm,
                    dens_j,fs_tm,model_opt,t_tm)


## CD4-MHC


rates_cm,hes_cm,l_cm = optimization(model,p_cm,
                            dens_j,fs_cm,model_opt,t_cm)

## TCR-MHC-CD4

u₀ = [0.8, 0., 0.2]
pinit = [1,1.e-2,1,1.e-2]
bimolecular_cons = 10 .^ rates_tm
model = trimolecular_diss_bell!
rates_tmc,hes_tmc,l_tmc = optimization(model,p_tmc,
                            dens_j,fs_tmc,model_opt,t_tmc)



p = 10 .^ rates
p = rates
p = [1.1385, -1.6272, -2.4367, -0.127488, 0.057953, -3.70439, 0.440775,-0.398936]
# cd4-TCR
p = rates
p = 10 .^ [0.41751687626687767, -0.23866855433293688, 0.06689410353314687, -0.5510289313898404, -0.09606981254493874, -0.1155499030032302, -3.9999999999999964, -1.7908513448187287]
 #tcr Mhc
p = rates_tm
p = 10 .^ rates_tm
p = 10 .^ [1.1854325378992236, -5.237609441955307, -0.2921252076084758, -2.703452812162055, 0.05179350352175943, -3.0727673775338427, -0.9157012289343508, -2.3711064932852954]

## ploting
dens = []

model = bimolecular_diss_bell!

gr()
p = rates_tm
p = 10 .^ rates_tm
p = [11.0,0.1,
        0.14,0.4,
        5.0,0.2,
        20.0,4.0]
show_diss_fit(model,tcr_mhc,fs_tm)
p1 = plot()
show_diss_fit!(p1,M,rates,data;cmap=nothing)
p = rates_cm
p = 10 .^ rates_cm
show_diss_fit(model,cd4_mhc,fs_cm)
p = rates_tmc
p = 10 .^ rates_tmc
p = [5.,1.,1.,0.1]
u₀ = [0.7, 0., 0.3]
model = trimolecular_diss_bell!
show_diss_fit(model,tcr_mhc_cd4,fs_tmc)



plotly()
function show_diss_fit(model, data, forces)
  plot()
  plot_diss_fit(model,forces)
  plot_dissociations(data)
end
function show_diss_fit(model, t, lnp, forces)
  plot()
  plot_diss_fit(model,forces)
  scatter!(t, lnp)
end

Fs = collect(0.:0.5:30.)
plot()
plot_force_lifetime_fit(model,Fs,tcr_mhc)
plot_force_lifetime_fit(model,Fs,cd4_mhc)


plot_F_vs_t(tcr_mhc)
plot_F_vs_t(cd4_mhc)



plot()
# p =[11.,0.,0.14,0.4,3,0.2,20,4]
plot_bell_equations(p)


max_t_tm  = extract_max_lifetimes(tcr_mhc)
scatter(fs_tm,max_t_tm)
max_lifetime_bounds(fs_tm,max_t_tm,0:0.1:30)
