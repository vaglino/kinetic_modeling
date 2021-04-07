edges = collect(0:5:30)
f_γϵ, t_γϵ, γϵ_bin_ids = bin_dissociations(γϵ,6;edges=edges)
f_δϵ, t_δϵ, δϵ_bin_ids = bin_dissociations(δϵ,6;edges=edges)
f_mix, t_mix, mix_bin_ids = bin_dissociations(mix,6;edges=edges)

s1=scatter(f_γϵ, t_γϵ, title="γϵ")
s2=scatter(f_δϵ, t_δϵ, title="δϵ")
s3=scatter(f_mix, t_mix, title="mix")
plot(s1,s2,s3)

plot(); Color = :blue
plot_F_vs_t(f_γϵ, t_γϵ)
Color = :red
plot_F_vs_t(f_δϵ, t_δϵ)
Color = :green
plot_F_vs_t(f_mix, t_mix)

Γ = Ft_average_and_error(f_γϵ, t_γϵ)
Δ = Ft_average_and_error(f_δϵ, t_δϵ)
Mix = Ft_average_and_error(f_mix, t_mix)

Γ_t = Γ.t_mean .± Γ.t_sem
Δ_t = Δ.t_mean .± Δ.t_sem
Mix_t = Mix.t_mean .± Mix.t_sem

sum_bi_t = 0.610*Γ_t + 0.390*Δ_t
# t_sum_bi = 0.5*Γ.t_mean + 0.5*Δ.t_mean
coop_means = (Mix_t - sum_bi_t) ./ sum_bi_t
plot(Γ.f_mean,coop_means)
plot!(0:35, lftm_coop,errorbar=:ribbon, xlabel="F",ylabel="lifetime cooperativity")


CSV.write("results/ft_mean_std_lftm_coop.csv", measurements2Df(Γ.f_mean,coop_means))

CSV.write("results/ft_mean_std_γϵ.csv",Ft_average_and_error(clean_dissociation_data(γϵ_binned)))
CSV.write("results/ft_mean_std_δϵ.csv", Ft_average_and_error(clean_dissociation_data(δϵ_binned)))
CSV.write("results/ft_mean_std_mix.csv", Ft_average_and_error(clean_dissociation_data(mix_binned)))
