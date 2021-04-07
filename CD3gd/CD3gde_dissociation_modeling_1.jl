## CD3 γϵ and δϵ cooperative binding to TCRαβ ecotodomain
##

# bimolecular dissociations
#       γϵ + TCR ← γϵ-TCR
#       δϵ + TCR ← δϵ-TCR

##
using Revise
include("model_fit_global.jl")
include("models.jl")
include("helper_functions.jl")
include("binning_helper_fs.jl")
using DataFrames, Plots
CD3_data = CSV.read("C:/Users/stravaglino3/Downloads/julia_code/CD3_data/dissociations.csv",DataFrame)

γϵ = delete_missing(CD3_data[:,2:3])
δϵ = delete_missing(CD3_data[:,8:9])
mix = delete_missing(CD3_data[:,14:15])
t_γϵ_thrml = delete_missing(CD3_data[:,5:5])
t_δϵ_thrml = delete_missing(CD3_data[:,11:11])
t_mix_thrml = delete_missing(CD3_data[:,17:17])
p_γϵ_thrml = survival_p.(t_γϵ_thrml)
p_δϵ_thrml = survival_p.(t_δϵ_thrml)
p_mix_thrml = survival_p.(t_mix_thrml)


s1=scatter(γϵ[1],γϵ[2], label="γϵ")
s2=scatter(δϵ[1],δϵ[2], label="δϵ")
s3=scatter(mix[1],mix[2], label="mix")
plot(s1,s2,s3)
plotly()
plot(t_γϵ_thrml,normalize_lnp.(p2lnp(p_γϵ_thrml)), label="γϵ")
plot!(t_δϵ_thrml,normalize_lnp.(p2lnp(p_δϵ_thrml)), label="δϵ")
plot!(t_mix_thrml,normalize_lnp.(p2lnp(p_mix_thrml)), label="mix")

## Thermal fluctuation bimolecular dissociation

tspan = (0.0,10.0)
t = collect(0.0:0.1:10)

gr()

u₀ = [1.0]          # Initial condition
tspan = (0.0,10.0)  # Simulation interval
# model_opt = (-5.,3.)
model_opt = (0.,100.)
pinit =[0.0]
dens = fill([],length(t_γϵ_thrml))
f_thrml = [[0.0]]
model = bimolecular_diss!

# CD3γϵ ------------------------------------------------------------------------
rates_γϵ_thrml,hes_γϵ_thrml,l_γϵ_thrml,boot_γϵ_thrml = optimization(model,p_γϵ_thrml,
                    dens,f_thrml,model_opt,t_γϵ_thrml)

# p = 10 .^ rates_γϵ_thrml
p = rates_γϵ_thrml
# p = [14.1]
k₋₁ᵒ = rates_γϵ_thrml
show_diss_fit(model,t_γϵ_thrml,p2lnp(p_γϵ_thrml),f_thrml)

# CD3δϵ ------------------------------------------------------------------------
rates_δϵ_thrml,hes_δϵ_thrml,l_δϵ_thrml,boot_δϵ_thrml = optimization(model,p_δϵ_thrml,
                    dens,f_thrml,model_opt,t_δϵ_thrml)
p = rates_δϵ_thrml
# p = [11.1]
k₋₂ᵒ = rates_δϵ_thrml
show_diss_fit(model,t_δϵ_thrml,p2lnp(p_δϵ_thrml),f_thrml)

# mix ------------------------------------------------------------------------
cons_mix_thrml = [[k₋₁ᵒ[1],k₋₂ᵒ[1]]]
cons_mix_thrml = [[14.1,11.1]]
u₀ = [0.4,0.4,0.2]          # Initial condition
model = cd3_dissociation!
rates_mix_thrml,hes_mix_thrml,l_mix_thrml,boot_mix_thrml = optimization(model,p_mix_thrml,
                    dens,cons_mix_thrml,model_opt,t_mix_thrml)
p = rates_mix_thrml
k₋₃ᵒ = rates_δϵ_thrml[1]
# f = (du,u,p,t) -> cd3_dissociation!(du,u,p,t,[],cons_mix_thrml[1])
# sol = solve_model(k₋₃ᵒ)
# sumsol =  sum(sol,dims=1)
show_diss_fit(model,t_mix_thrml,p2lnp(p_mix_thrml),cons_mix_thrml)

## Load Force-dependent bimolecular dissociation

f_γϵ, t_γϵ, γϵ_bin_ids = bin_dissociations(γϵ,6;bin_type="force")
f_δϵ, t_δϵ, δϵ_bin_ids = bin_dissociations(δϵ,6)
f_mix, t_mix, mix_bin_ids = bin_dissociations(mix,6)
f_γϵ, t_γϵ = resort(f_γϵ, t_γϵ)
f_δϵ, t_δϵ = resort(f_δϵ, t_δϵ)
f_mix, t_mix = resort(f_mix, t_mix)
p_γϵ = survival_p.(t_γϵ)
p_δϵ = survival_p.(t_δϵ)
p_mix = survival_p.(t_mix)

s1=scatter(f_γϵ, t_γϵ, title="γϵ")
s2=scatter(f_δϵ, t_δϵ, title="δϵ")
s3=scatter(f_mix, t_mix, title="mix")
plot(s1,s2,s3)

s1=scatter(t_γϵ, p2lnp(p_γϵ), title="γϵ")
s2=scatter(t_δϵ, p2lnp(p_δϵ), title="δϵ")
s3=scatter(t_mix, p2lnp(p_mix), title="mix")
plot(s1,s2,s3)


using StatsBase
t_γϵ_avg = mean.(t_γϵ)
f_γϵ_avg = mean.(f_γϵ)
plot(f_γϵ_avg,t_γϵ_avg)



t_γϵ[1]
prob = survival_p_2(t_γϵ[1]) 

survival_p_2(t) = collect((length(t)-1 :-1:0) / length(t))
plot(t_γϵ[1],p2lnp(prob))
mean(t_γϵ[2])
sum(t_γϵ[1] .* prob ./ sum(prob))
sum(t_γϵ[2].*survival_p_2(t_γϵ[2]))

# u₀ = [1.0]          # Initial condition
# tspan = (0.0,15.0)  # Simulation interval
# # model_opt = (-5.,3.)
# model_opt = (0.,100.)
# pinit =[1.0,1.0]
# dens = fill([],length(t_γϵ))
# model = bimolecular_diss_bell_single!
# cons_γϵ = [[k₋₁ᵒ[1], f] for f in f_γϵ_avg]
#
# rates_γϵ,hes_γϵ,l_γϵ,boot_γϵ = optimization(model,p_γϵ,
#                     dens,f_γϵ_avg,model_opt,t_γϵ)
#
# p = rates_γϵ
# show_diss_fit(model, t_γϵ, p2lnp(p_γϵ), f_γϵ_avg)
# plot_F_vs_t(f_γϵ,t_γϵ)


## try bimolecular dissociation bell model with activated state

u₀ = [1.0,0.0]
tspan = (0.0,15.0)
# model_opt = (-5.,3.)  # Simulation interval
model_opt = (0.,100.)
pinit = ones(8)
dens = fill([],length(t_γϵ))
model = bimolecular_diss_bell!
rates_γϵ,hes_γϵ,l_γϵ,boot_γϵ = optimization(model,p_γϵ,
                    dens,f_γϵ_avg,model_opt,t_γϵ)

p = rates_γϵ
# p = 10 .^ rates_γϵ
show_diss_fit(model, t_γϵ, p2lnp(p_γϵ), f_γϵ_avg)
Fs = collect(0.:0.5:30.)
plot()
plot_F_vs_t(f_γϵ,t_γϵ)
plot_force_lifetime_fit(model,Fs,f_γϵ,t_γϵ)

## try to include thermal fluctuation data to bfp data


p_γϵ_all = vcat(p_γϵ_thrml, p_γϵ)
t_γϵ_all = vcat(t_γϵ_thrml, t_γϵ)
f_γϵ_thrml = [zeros(size(t_γϵ_thrml[1]))]
f_γϵ_all = vcat(f_γϵ_thrml,f_γϵ)
f_γϵ_avg_all = vcat(0.0,f_γϵ_avg)
dens = fill([],length(t_γϵ_all))

rates_γϵ,hes_γϵ,l_γϵ,boot_γϵ = optimization(model,p_γϵ_all,
                    dens,f_γϵ_avg_all,model_opt,t_γϵ_all)

p = rates_γϵ
# p = 10 .^ rates_γϵ
show_diss_fit(model, t_γϵ_all, p2lnp(p_γϵ_all), f_γϵ_avg_all)
Fs = collect(0.:0.5:30.)
plot()
plot_F_vs_t(f_γϵ_all,t_γϵ_all)
plot_force_lifetime_fit(model,Fs,f_γϵ,t_γϵ)
