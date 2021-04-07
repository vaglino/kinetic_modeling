using Measurements
r_ini = [14.0,0.1,
        0.14,0.4,
        5.0,0.2,
        20.0,4.0]

model = bimolecular_diss_bell_activated_state_w_cons!
# u₀ = [1.0,0.0]
u₀ = [] # don't provide initial states
u₀_type = "equilibrium" # calculate initial states based on equilibrium
# u₀ = "detailed_balance"
# u₀ = "force"
# u₀ = []
tspan = (0.0,maximum(γϵ.t))
# k₋ᵒf_γϵ = 14.1
k₋ᵒf_γϵ = 12.2
γϵ_ft_M = M(model,u₀,u₀_type,tspan,[],[k₋ᵒf_γϵ],0.0,Rodas5)
# γϵ_ft_M = M(model,u₀,u₀_type,tspan,[],[0.5],0.0,Rodas5)

pinit = [   0.01,
        0.1,0.6,
        3.,0.2,
        150.0,4.0]
opt =  [        (0.,10.),
        (0.,1.),(0.,10.),
        (0.,10.),(0.,10.),
        (0.,500.),(0.,10.)]


@time mle_loss(pinit,γϵ,γϵ_ft_M)
# loss_in = (p) -> mle_loss(p,γϵ,γϵ_ft_M)
# hes = ForwardDiff.hessian(loss_in,rates_γϵ)
# σ_γϵ = hessian2σ(hes)

@time rates_γϵ, l_γϵ, hes_γϵ = optimization(γϵ_ft_M,γϵ,opt)
ks_γϵ = rates_γϵ
σ_γϵ = hessian2σ(hes_γϵ)
k_γϵ = rates_γϵ .± σ_γϵ

plotly()
plot_ft_fit(γϵ_ft_M,ks_γϵ,γϵ)
# plot_ft_fit(γϵ_ft_M,ks_γϵ,γϵ,dF=30)
plot(); Color = :blue
plot_mean_ft_fit(γϵ_ft_M,k_γϵ,γϵ,γϵ_binned;species="γϵ")
# sols_γϵ = ftsolve(γϵ_ft_M,k_γϵ) #change solver to Tsit5
# plot(Fs,mean_lifetime.(sols_γϵ),ylabel="<t>",xlabel="F",title="γϵ")
# plot_F_vs_t(clean_dissociation_data(γϵ_binned))
u₀γϵ = initial_u(ks_γϵ,γϵ_ft_M)
________________________________________________________________________________________________
## δϵ fitting with 2 species model
model = bimolecular_diss_bell_activated_state_w_cons!
u₀ = []
u₀_type = "equilibrium"
tspan = (0.0,maximum(δϵ.t))
# k₋ᵒf_δϵ = 11.1
k₋ᵒf_δϵ = 13.8
δϵ_ft_M = M(model,u₀,u₀_type,tspan,[],[k₋ᵒf_δϵ],0.0,Rodas5)
pinit = [   0.1,
        0.1,0.5,
        0.1,0.5,
        10.0,5.0]
opt =  [        (0.,10.),
        (0.,1.),(0.,10.),
        (0.,10.),(0.,10.),
        (0.,100.),(0.,10.)]

@time mle_loss(pinit,δϵ,δϵ_ft_M)

@time rates_δϵ, l_δϵ, hes_δϵ = optimization(δϵ_ft_M,δϵ,opt)
ks_δϵ = rates_δϵ
σ_δϵ = hessian2σ(hes_δϵ)
k_δϵ = rates_δϵ .± σ_δϵ
u₀δϵ = initial_u(ks_δϵ,δϵ_ft_M)

plot_ft_fit(δϵ_ft_M,ks_δϵ,δϵ)
plot();
Color = :red
plot_mean_ft_fit(δϵ_ft_M,k_δϵ,δϵ,δϵ_binned;species="δϵ")
# plot(Fs,mean_lifetime.(sols_δϵ),ylabel="<t>",xlabel="F",title="δϵ")
# plot_F_vs_t(clean_dissociation_data(δϵ_binned))

## mix fitting with 2 species model

model = cd3_tri_dissociation_bell!
# u₀_type = "equilibrium"
u₀_type = []
# u₀ = vcat(u₀γϵ./2, u₀δϵ./2, 0.0,0.0)
u₀ = vcat(0.454*u₀γϵ,0.290*u₀δϵ, 0.256*[1.,0.])
# u₀ = vcat(0.454*u₀γϵ,0.290*u₀δϵ, 0.256*[0.,0.])/(1-0.256)
sum(u₀)
# u₀ = vcat(u₀γϵ/3, u₀δϵ/3)
tspan = (0.0,maximum(mix.t))
k₋ᵒf_mix = 0.6
cons = [k₋ᵒf_γϵ;
        ks_γϵ;
        k₋ᵒf_δϵ;
        ks_δϵ;
        k₋ᵒf_mix ]
mix_ft_M = M(model,u₀,u₀_type,tspan,[],cons,0.0,Rodas5)


pinit = [   0.1,
        0.1,0.5,
        0.1,0.5,
        10.0,5.0]
# pinit = zeros(7)
opt =  [        (0.,10.),
        (0.,1.),(0.,1.),
        (0.,10.),(0.,1.),
        (0.,100.),(0.,10.)]

@time mle_loss(pinit,mix,mix_ft_M)
loss_in = (p) -> mle_loss(p,mix,mix_ft_M)
hes_mix = ForwardDiff.hessian(loss_in,pinit)
# hes = Zygote.hessian(loss_in,pinit)

@time rates_mix, l_mix, hes_mix = optimization(mix_ft_M,mix,opt)
ks_mix = rates_mix
σ_mix = hessian2σ(hes_mix)
k_mix = rates_mix .± σ_mix


ts = 0:dt:maximum(mix.t)
Fs = 0:1:maximum(mix.F)
plot_ft_fit(mix_ft_M,ks_mix,mix)
# sols_mix = plot_ft_fit(mix_ft_M,pinit,mix)


plot_ft_fit(mix_ft_M,ks_mix,mix)
plot(); Color = :green
plot_mean_ft_fit(mix_ft_M,k_mix,mix,mix_binned;species="mix")
# plot(Fs,mean_lifetime.(sols_mix),ylabel="<t>",xlabel="F",title="mix")
# plot_F_vs_t(clean_dissociation_data(mix_binned))


## plot rates

plt = plot_bell_equations(vcat(k₋ᵒf_γϵ, ks_γϵ))
title!("γϵ")
plt = plot_bell_equations(vcat(k₋ᵒf_δϵ, ks_δϵ))
title!("δϵ")
plt = plot_bell_equations(vcat(k₋ᵒf_mix, ks_mix))
title!("mix")

pyplot()
plot(); Color = :blue
plot_mean_ft_fit(γϵ_ft_M,k_γϵ,γϵ,γϵ_binned;species="γϵ")
Color = :red
plot_mean_ft_fit(δϵ_ft_M,k_δϵ,δϵ,δϵ_binned;species="δϵ")
Color = :green
plot_mean_ft_fit(mix_ft_M,k_mix,mix,mix_binned;species="mix")

plot(); Color = :blue
extrapolate_fit(γϵ_ft_M,ks_γϵ,γϵ,150.0;species="γϵ",dF=10.0)
Color = :red
extrapolate_fit(δϵ_ft_M,ks_δϵ,δϵ,150.0;species="δϵ",dF=10.0)
Color = :green
extrapolate_fit(mix_ft_M,ks_mix,mix,150.0;species="mix",dF=10.0)

mean_t_γϵ = mean_ft_fit(γϵ_ft_M,k_γϵ,γϵ,35;species=[],dt=0.1,dF=1.0)
mean_t_δϵ = mean_ft_fit(δϵ_ft_M,k_δϵ,δϵ,35;species=[],dt=0.1,dF=1.0)
mean_t_mix = mean_ft_fit(mix_ft_M,k_mix,mix,35;species=[],dt=0.1,dF=1.0)

mean_t_sum_bi = 0.610*mean_t_γϵ + 0.390*mean_t_δϵ
plot!(mean_t_sum_bi)

pyplot()
lftm_coop = (mean_t_mix -mean_t_sum_bi) ./ mean_t_sum_bi
plot(0:35, 100 .* lftm_coop,errorbar=:ribbon, xlabel="F",ylabel="lifetime cooperativity (%)")

CSV.write("results/γϵ_fit.csv", measurements2Df(0:35,mean_t_γϵ))
CSV.write("results/δϵ_fit.csv", measurements2Df(0:35,mean_t_δϵ))
CSV.write("results/mix_fit.csv", measurements2Df(0:35,mean_t_mix))
CSV.write("results/lftm_coop.csv", measurements2Df(0:35,lftm_coop))
