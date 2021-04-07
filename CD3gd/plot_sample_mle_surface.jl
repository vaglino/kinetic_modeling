
M_example = M(model,u₀,u₀_type,tspan,[],[10],15,Rodas5)

rates_example = ks_γϵ
sol = solve_prob(M_example,ks_γϵ;tᵢ=[])


ŷ = sum(sol,dims=1)
mean_t = mean_lifetime_linear(sol)

L = ŷ / mean_t

pyplot()
# plot(sol)
plot(sol.t,ŷ',linewidth=5,label="P")
plot!(sol.t,L',linewidth=5,label="pdf")
xlabel!("t")
ylabel!("distribution")
