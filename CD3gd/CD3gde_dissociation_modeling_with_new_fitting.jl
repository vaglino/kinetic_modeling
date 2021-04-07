## CD3 γϵ and δϵ cooperative binding to TCRαβ ecotodomain
##

# bimolecular dissociations
#       γϵ + TCR ← γϵ-TCR
#       δϵ + TCR ← δϵ-TCR

##
using Revise
# include("model_fit_global.jl")
include("models.jl")
include("models_dissociation.jl")
include("helper_functions.jl")
include("binning_helper_fs.jl")
include("fitting.jl")
include("helper_plotting.jl")
using DataFrames, Plots

## load thermal data

CD3_data = CSV.read("C:/Users/stravaglino3/Downloads/julia_code/CD3_data/dissociations.csv",DataFrame)

t_γϵ_thrml = delete_missing(CD3_data[:,5:5])
t_δϵ_thrml = delete_missing(CD3_data[:,11:11])
t_mix_thrml = delete_missing(CD3_data[:,17:17])
p_γϵ_thrml = survival_p.(t_γϵ_thrml)
p_δϵ_thrml = survival_p.(t_δϵ_thrml)
p_mix_thrml = survival_p.(t_mix_thrml)


# γϵ_thrml = DissData(p_γϵ_thrml,t_γϵ_thrml,[[]])
# δϵ_thrml = DissData(p_δϵ_thrml,t_δϵ_thrml,[[]])
# mix_thrml = DissData(p_mix_thrml,t_mix_thrml,[[]])

γϵ_thrml = DissData(p_γϵ_thrml,rescale_t.(t_γϵ_thrml),[[]])
δϵ_thrml = DissData(p_δϵ_thrml,rescale_t.(t_δϵ_thrml),[[]])
mix_thrml = DissData(p_mix_thrml,rescale_t.(t_mix_thrml),[[]])

## Thermal fluctuation bimolecular dissociation

gr()

u₀ = [1.0]          # Initial condition
tspan = (0.0,10.0)  # Simulation interval
dens = fill([],length(t_γϵ_thrml))
model = bimolecular_diss!
cons_thrml = [[0.0]]

#pack parameters
γϵM_thrml = M(model,u₀,tspan,dens,cons_thrml)

opts = (0.,100.)
pinit = [0.0]

# CD3γϵ thermal fitting---------------------------------------------------------
rates_γϵ_thrml,l_γϵ_thrml = optimization(γϵM_thrml, γϵ_thrml, opts)
k₋₁ᵒ = rates_γϵ_thrml


# CD3δϵ thermal fitting---------------------------------------------------------
δϵM_thrml = M(model,u₀,tspan,dens,cons_thrml)
rates_δϵ_thrml,l_δϵ_thrml = optimization(δϵM_thrml, δϵ_thrml, opts)
k₋₂ᵒ = rates_δϵ_thrml


# CD3 mix thermal fitting-------------------------------------------------------
u₀ = [0.4,0.4,0.2]
cons_mix_thrml = [[k₋₁ᵒ[1],k₋₂ᵒ[1]]]
cons_mix_thrml = [[14.1,11.1]]
model = cd3_dissociation!

opts = (0.,100.)
pinit = [0.0,0.0]
mixM_thrml = M(model,u₀,tspan,dens,cons_mix_thrml)

rates_mix_thrml,l_mix_thrml = optimization(mixM_thrml, mix_thrml, opts)

pyplot()
plotly()
# thermal fluctuation plotting
# gr()
p1 = plot()
show_diss_fit!(p1,γϵM_thrml,rates_γϵ_thrml, γϵ_thrml; cmap=:blue)
show_diss_fit!(p1,δϵM_thrml,rates_δϵ_thrml, δϵ_thrml; cmap=:red)
show_diss_fit!(p1,mixM_thrml,rates_mix_thrml, mix_thrml; cmap=:green)

show_diss_fit!(p1,γϵM_thrml,[14.1], γϵ_thrml; cmap=:blue)
show_diss_fit!(p1,δϵM_thrml,[11.1], δϵ_thrml; cmap=:red)
u3 = rates_mix_thrml[1]
mixM_thrml.u₀ = [(1-u3)*0.610, (1-u3)*0.390, u3]
show_diss_fit!(p1,mixM_thrml,rates_mix_thrml[2], mix_thrml; cmap=:green)

## load force dependent data

γϵ_data = CSV.read("C:/Users/stravaglino3/Downloads/julia_code/CD3_data/cd3ge-tcr-unbinned.txtbinned bond lifetime.csv",DataFrame)
δϵ_data = CSV.read("C:/Users/stravaglino3/Downloads/julia_code/CD3_data/cd3de-tcr-unbinned.txtbinned bond lifetime.csv",DataFrame)
mix_data = CSV.read("C:/Users/stravaglino3/Downloads/julia_code/CD3_data/cd3mix-tcr-unbinned.txtbinned bond lifetime.csv",DataFrame)

γϵ = clean_dissociation_data(γϵ_data)
δϵ = clean_dissociation_data(δϵ_data)
mix = clean_dissociation_data(mix_data)

t_γϵ, lnp_γϵ = separate_diss_data(γϵ)
t_δϵ, lnp_δϵ = separate_diss_data(δϵ)
t_mix, lnp_mix = separate_diss_data(mix)

p_γϵ = lnp2p.(lnp_γϵ)
p_δϵ = lnp2p.(lnp_δϵ)
p_mix = lnp2p.(lnp_mix)

using StatsBase
f_γϵ_avg = extract_avg_forces(γϵ)
f_δϵ_avg = extract_avg_forces(δϵ)
f_mix_avg = extract_avg_forces(mix)

γϵ = DissData(p_γϵ,t_γϵ,f_γϵ_avg)
δϵ = DissData(p_δϵ,t_δϵ,f_δϵ_avg)
mix = DissData(p_mix,t_mix,f_mix_avg)
# mix_thrml = DissData(p_mix_thrml,t_mix_thrml,[[]])

## force dependent dissociation


u₀ = [1.0, 0.0]          # Initial condition
tspan = (0.0,20.0)       # Simulation interval
dens = fill([],length(t_γϵ))
model = bimolecular_diss_bell_activated_state_w_cons!

# CD3γϵ force-dependent dissociations-------------------------------------------
cons_γϵ = [[k₋₁ᵒ[1], f] for f in f_γϵ_avg]
cons_γϵ = [[14.1, f] for f in f_γϵ_avg]
#pack parameters
γϵM = M(model,u₀,tspan,dens,cons_γϵ)

opts = (0.,100.)
pinit = zeros(7)

# # CD3γϵ force fitting
rates_γϵ,l_γϵ = optimization(γϵM, γϵ, opts)

pyplot()
p2 = plot()
show_diss_fit!(p2,γϵM,rates_γϵ, γϵ)

# CD3δϵ force-dependent dissociations-------------------------------------------

cons_δϵ = [[11.1, f] for f in f_δϵ_avg]
#pack parameters
δϵM = M(model,u₀,tspan,dens,cons_δϵ)

opts = (0.,100.)
pinit = zeros(7)

# # CD3γϵ force fitting
rates_δϵ,l_δϵ = optimization(δϵM, δϵ, opts)

p3 = plot()
show_diss_fit!(p3,δϵM,rates_δϵ, δϵ)

# plot()
# plot_F_vs_t(clean_dissociation_data(δϵ_data))


# CD3mix force-dependent dissociations-------------------------------------------



u₀ = [0.4,0.0, 0.4,0.0, 0.2,0.0]          # Initial condition
tspan = (0.0,20.0)       # Simulation interval
dens = fill([],length(t_mix))
model = cd3_tri_dissociation_bell!

k₋₃ᵒ = rates_mix_thrml[1]
cons_mix = [[vcat(k₋₁ᵒ,rates_γϵ), vcat(k₋₂ᵒ,rates_δϵ), k₋₃ᵒ,  f] for f in f_mix_avg]

#pack parameters
mixM = M(model,u₀,tspan,dens,cons_mix)

opts = (0.,100.)
opts = [(0.,0.1),(0.,100.),(0.,0.1),(0.,100.),(0.,0.1),(0.,100.),(0.,0.1)]
pinit = zeros(7)

# # CD3γϵ force fitting
rates_mix,l_mix = optimization(mixM, mix, opts)
p4 = plot()
show_diss_fit!(p4,mixM,rates_mix, mix)
show_diss_fit!(p4,mixM, zeros(7), mix)



ft1=plot()
plot_F_vs_t(clean_dissociation_data(γϵ_data))
title!("γϵ")
ft2=plot()
plot_F_vs_t(clean_dissociation_data(δϵ_data))
title!("δϵ")
ft3=plot()
plot_F_vs_t(clean_dissociation_data(mix_data))
title!("mix")
plot(ft1,ft2,ft3)
