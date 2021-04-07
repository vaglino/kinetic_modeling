using Revise
using StatsBase, KernelDensity, Distributions
dir = "C:/Users/stravaglino3/Downloads/julia_code/src/"
includet(dir*"models.jl")
includet(dir*"models_dissociation.jl")
includet(dir*"helper_functions.jl")
includet(dir*"binning_helper_fs.jl")
includet(dir*"mle_fitting.jl")
includet(dir*"mle_plotting.jl")
using DataFrames, Plots, ColorSchemes
data_dir = "C:/Users/stravaglino3/Downloads/julia_code/CD3gd/data/"
CD3_data = CSV.read(data_dir*"dissociations_cleaned.csv",DataFrame)

γϵ_binned = CSV.read(data_dir*"data/rebinned/cd3ge-tcr-rebin binned bond lifetime.csv",DataFrame)
δϵ_binned = CSV.read(data_dir*"cd3de-tcr-unbinned.txtbinned bond lifetime.csv",DataFrame)
mix_binned = CSV.read(data_dir*"data/rebinned/cd3mix-tcr-rebin.txtbinned bond lifetime.csv",DataFrame)


γϵ = delete_missing(CD3_data[:,2:3])
δϵ = delete_missing(CD3_data[:,8:9])
mix = delete_missing(CD3_data[:,14:15])

pyplot()
gr()
s1=Plots.scatter(γϵ[1],γϵ[2], label="γϵ");
s2=Plots.scatter(δϵ[1],δϵ[2], label="δϵ");
s3=Plots.scatter(mix[1],mix[2], label="mix");
plot(s1,s2,s3)

γϵ  = F_t_Data(γϵ[1],γϵ[2],   InterpKDE(kde(γϵ[1])))
δϵ  = F_t_Data(δϵ[1],δϵ[2],   InterpKDE(kde(δϵ[1])))
mix = F_t_Data(mix[1],mix[2], InterpKDE(kde(mix[1])))

h1=histogram2d(γϵ.F, γϵ.t, nbins=(40,50),title="γϵ")#, c=:jet);#,aspect_ratio=2)
h2=histogram2d(δϵ.F, δϵ.t, nbins=(40,60),title="δϵ")#, c=:turbo);#,aspect_ratio=2)
h3=histogram2d(mix.F, mix.t, nbins=(40,50),title="mix")#, c=:turbo);#,aspect_ratio=2)
plot(h1,h2,h3)

#________________________________________________________________________________________________
## try fitting with 2 species model
# rates from Dunn's Science paper:
r_ini = [14.0,0.1,
        0.14,0.4,
        5.0,0.2,
        20.0,4.0]

# γϵ.t = γϵ.t .- minimum(γϵ.t)

model = bimolecular_diss_bell!
# u₀ = [1.0,0.0]
u₀ = [] # don't provide initial states
u₀_type = "equilibrium" # calculate initial states based on equilibrium
# u₀ = "detailed_balance"
# u₀ = "force"
# u₀ = []
tspan = (0.0,maximum(γϵ.t))
γϵ_ft_M = M(model,u₀,u₀_type,tspan,[],[],0.0)

pinit = [14.1,0.01,
        0.1,0.6,
        3.,0.2,
        150.0,4.0]

opt =  [(14.05,14.15),(0.,10.),
        (0.,1.),(0.,10.),
        (0.,10.),(0.,10.),
        (0.,500.),(0.,10.)]


# pinit = rand(8)*10
@time mle_loss(pinit,γϵ,γϵ_ft_M)
#=
res = so.minimize(x->mle_loss(10 .^ x,γϵ,γϵ_ft_M), log10.(pinit),method="COBYLA",
                        bounds=opt,
                        options=Dict("maxiter"=>10000,
                                 "disp" => true) )
res["x"]
=#   # the random seed

@time rates_γϵ, l_γϵ, hes_γϵ = optimization(γϵ_ft_M,γϵ,opt)
using LinearAlgebra
errors = sqrt.(abs.( diag(inv(hes_γϵ)) ))
# ks_γϵ = 10 .^ rates_γϵ
ks_γϵ = rates_γϵ
ks_γϵ = pinit

# pyplot()
plotly()
dt = 0.1
ts = 0:dt:maximum(γϵ.t)
Fs = 0:1:maximum(γϵ.F)
# Fs = 0:5:150
ks_γϵ = [14.1,0.0,
        0.1,0.6,
        3.,0.2,
        150.0,4.0]

sols_γϵ = plot_ft_fit(γϵ_ft_M,ks_γϵ,γϵ)

plot(Fs,mean_lifetime.(sols_γϵ),ylabel="<t>",xlabel="F",title="γϵ")
plot_F_vs_t(clean_dissociation_data(γϵ_binned))

us = map(x -> initial_u(ks_γϵ,γϵ_ft_M,x), 0:30)
scatter(0:30,reduce(hcat,us)')

plot_bell_equations([11.0,0.0,
                     0.14,0.4,
                     3.0,0.2,
                     20.0,4.0])
plt = plot_bell_equations(vcat(k₋ᵒf_γϵ, ks_γϵ))
# display(plt)
pyplot()
u₀_γϵ = initial_state(ks_γϵ[1],ks_γϵ[3],ks_γϵ[5],ks_γϵ[7])
initial_state(ks_γϵ[5],ks_γϵ[7])
a
________________________________________________________________________________________________
## δϵ fitting with 2 species model
model = bimolecular_diss_bell!
# u₀ = [1.0,0.0]
u₀ = "equilibrium"
tspan = (0.0,maximum(δϵ.t))
δϵ_ft_M = M(model,u₀,tspan,[],[],0.0)
pinit = [11.1,0.1,
        0.1,0.5,
        0.1,0.5,
        10.0,5.0]
opt =  [(11.05,11.15),(0.,10.),
        (0.,1.),(0.,10.),
        (0.,10.),(0.,10.),
        (0.,100.),(0.,10.)]
@time rates_δϵ, l_δϵ, hes_δϵ = optimization(δϵ_ft_M,δϵ,opt)

ks_δϵ = rates_δϵ

ks_δϵ = [11.087136651210209, 0.0026318767020595732,
        0.02256859725182868, 0.5066206884832165,
        0.21328859839615166, 0.058733617801501145,
        1.5349514553341914, 0.6116955164634462]

plotly()
ts = 0:dt:maximum(δϵ.t)
Fs = 0:1:maximum(δϵ.F)
sols_δϵ = plot_ft_fit(δϵ_ft_M,ks_δϵ,δϵ)
pyplot()
plotly()
plot(Fs,mean_lifetime.(sols_δϵ),ylabel="<t>",xlabel="F",title="δϵ")
plot_F_vs_t(clean_dissociation_data(δϵ_binned))


pyplot()
plot_bell_equations(ks_δϵ)
u₀_δϵ = initial_state(ks_δϵ[1],ks_δϵ[3],ks_δϵ[5],ks_δϵ[7])



________________________________________________________________________________________________
## mix fitting with 2 species model

model = cd3_tri_dissociation_bell!
u₀ = vcat(u₀_γϵ./3, u₀_δϵ./3, 1/3,0.0)
u₀ = vcat(u₀_γϵ./2, u₀_δϵ./2, 0.0,0.0)
tspan = (0.0,maximum(mix.t))
mix_ft_M = M(model,u₀,tspan,[],[])
r1 = ks_γϵ
r2 = ks_δϵ
pinit = [0.6,0.1,
        0.1,0.5,
        0.1,0.5,
        10.0,5.0]
opt =  [(0.59,0.61),(0.,10.),
        (0.,1.),(0.,1.),
        (0.,10.),(0.,1.),
        (0.,100.),(0.,10.)]
@time rates_mix, l_mix = optimization(mix_ft_M,mix,opt)

ks_mix = rates_mix
ks_mix = pinit
ts = 0:dt:maximum(mix.t)
Fs = 0:1:35
h2 = plot_ft_fit(mix_ft_M,ks_mix,mix)







histogram(γϵ.F,bins=10)

h = fit(Histogram, rand(100))
plot(h)
U = kde(γϵ.F;kernel=Distributions.Normal)
ik = InterpKDE(U)

pdf(U, [1,2,3])
plot(U.x,U.density)
scatter!(γϵ.F, zeros(size(γϵ.F)), markersize=7,markeralpha=0.3)

size(γϵ.F)

X = randn(100)
Y = randn(100)

k = kde(X)

pdf(k, k.x) ≈ k.density

ik = InterpKDE(k)

pdf(ik,2)
