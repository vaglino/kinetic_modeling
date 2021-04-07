densities_T = [25.0, 38.0, 1.0] # TCR, pMHC, Ac
densities_TC1 = [15.0, 38.0, 10.0, 1.0] # TCR, pMHC, CD4, Ac
densities_TC2 = [8.0,  38.0, 17.0, 1.0]
densities_TC3 = [3.0,  38.0, 22.0, 1.0]

rates_T = [-3.37638358011457
            -0.302013187183543]
rates_TC1 = [-1.72689171836867
            -0.646695212785203]
rates_TC2 = [-0.644293816490510
            0.617796010276307]
rates_TC3 = [0.891017383250648
            2.16934160248305]

sd_data = CSV.read("C:/Users/stravaglino3/Downloads/julia_code/n_stds_for_julia.csv")
nT_sds = sd_data[:,1]
nTC1_sds = Array(sd_data[:,2:3])
nTC2_sds = Array(sd_data[:,4:5])
nTC3_sds = Array(sd_data[:,6:7])

scatter(n_se_tcr[:,1],n_se_tcr[:,2],yerror=n_se_tcr[:,3],markercolor=:black,markersize=5)
scatter!(n_se_1[:,1],n_se_1[:,2],yerror=n_se_1[:,3],markercolor=:blue,markersize=5)
scatter!(n_se_2[:,1],n_se_2[:,2],yerror=n_se_2[:,3],markercolor=:red,markersize=5)
h = scatter!(n_se_3[:,1],n_se_3[:,2],yerror=n_se_3[:,3],markercolor=:cyan,markersize=5,legend=false)


tspan = (0.0,16.0)
t = collect(0.01:0.01:16)

gr()
# bimolecular TCR-pMHC bonds
u₀ = [0.0]
f = (du,u,p,t) -> bimolecular!(du,u,p,t,densities_T,[])
n_T = Array(solve_model(rates_T))'[2:end]
n_T_sd = n_T .± nT_sds

# TCR-pMHC-CD4 (1)
u₀ = [0.0, 0.0]
f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_TC1,10 .^ rates_T)
n_TC1 = Array(solve_model(rates_TC1))'[2:end,:]
n_TC1_sd = sum(n_TC1 .± nTC1_sds, dims=2)

# TCR-pMHC-CD4 (2)
f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_TC2,10 .^ rates_T)
n_TC2 = Array(solve_model(rates_TC2))'[2:end,:]
n_TC2_sd = sum(n_TC2 .± nTC2_sds, dims=2)

# TCR-pMHC-CD4 (3)
f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_TC3,10 .^ rates_T)
n_TC3 = Array(solve_model(rates_TC3))'[2:end,:]
n_TC3_sd = sum(n_TC3 .± nTC3_sds, dims=2)


plot!(t, n_T_sd, errorbars=:ribbon, alpha=0.1,msc=:black)
plot!(t, n_T, lw=3, lc=:black)
plot!(t,n_TC1_sd, errorbars=:ribbon, alpha=0.1, msc=:blue)
plot!(t, sum(n_TC1,dims=2), lw=3, lc=:blue)
plot!(t,n_TC2_sd, errorbars=:ribbon, alpha=0.1, msc=:red)
plot!(t, sum(n_TC2,dims=2), lw=3, lc=:red)
plot!(t,n_TC3_sd, errorbars=:ribbon, alpha=0.1, msc=:cyan)
plot!(t, sum(n_TC3,dims=2), lw=3, lc=:cyan)

xlabel!("Contact Time tc")
ylabel!("<n>")

#------------------------------------------------------------------------------
# CD4 fitting

cd4_1 = CSV.read("cd4_data_from_muaz.csv")
t_cd41, Pa_cd41, n_se_cd41 = clean_data(cd4_1)
n_cd41 = Pa2n(Pa_cd41)

f = (du,u,p,t) -> bimolecular!(du,u,p,t,dens,cons) # enclose constants

u₀ = [0.0]          # Initial condition
tspan = (0.0,8.0)  # Simulation interval
pinit = [1e-4,0.5]      # equation parameter. p = [k₁, k₋₁]
densities_C = [1700., 850., 1.] #CD4, pMHC, A
cons = []

rates_C,hes_C,l_C,boot_C = optimization(bimolecular!,n_cd41,densities_C',[],model_opt2,t_cd41)

# Setup the ODE problem, then solve
# plotly()
t = collect(0.01:0.01:16)
tspan = (0.0,16.0)
f = (du,u,p,t) -> bimolecular!(du,u,p,t,densities_C,cons)
prob = ODEProblem(f,u₀,tspan,10 .^ rates_C)
sol_C = Array(solve(prob,Vern7(),saveat=t,abstol=1e-10,reltol=1e-10))'
plot(t,sol_C,lc=:orange,lw=3)
scatter!(n_se_cd41[:,1],n_se_cd41[:,2],yerror=n_se_cd41[:,3],markercolor=:orange)


#-----------------------------------------------------------------------------

densities_T =   [25.0, 38.0, 1.0]       # TCR, pMHC, Ac
densities_TC1 = [15.0, 38.0, 10.0, 1.0] # TCR, pMHC, CD4, Ac
densities_TC2 = [8.0,  38.0, 17.0, 1.0]
densities_TC3 = [3.0,  38.0, 22.0, 1.0]
densities_C =   [1700., 850., 1.]       # CD4, pMHC, A


n_T_norm = n_T ./25. .*15.
n_TC1_norm = sum(n_TC1,dims=2)
n_TC2_norm = sum(n_TC2,dims=2) ./8. .*15.# ./17. .*10.
n_TC3_norm = sum(n_TC3,dims=2) ./3. .*15.# ./22. .*10.
n_C_norm = sol_C ./1700. .*10. ./850. .*38

plotly()
t = collect(0.01:0.01:16)
plot(t,n_T_norm,lw=3,lc=:black,label="TCR")
plot!(t,n_TC1_norm,lw=3,lc=:blue,label="TCR 15, CD4 10")
plot!(t,n_TC2_norm,lw=3,lc=:red,label="TCR 8, CD4 17")
plot!(t,n_TC3_norm,lw=3,lc=:cyan,label="TCR 3, CD4 22")
plot!(t,n_C_norm,lw=3,lc=:orange,yscale=:linear,label="CD4",legend=:bottomright)
# yaxis!("log<n>")
yaxis!("<n>")
xlabel!("Contact Time tc")
# savefig("n_log_scaled")


coop_indx_TC1 = cooperativity_index.(n_T_norm,n_C_norm,n_TC1_norm)
coop_indx_TC2 = cooperativity_index.(n_T_norm,n_C_norm,n_TC2_norm)
coop_indx_TC3 = cooperativity_index.(n_T_norm,n_C_norm,n_TC3_norm)
plot(t,coop_indx_TC1 .* 100,lw=3,lc=:blue,label="Coop TCR 15, CD4 10")
plot!(t,coop_indx_TC2 .* 100,lw=3,lc=:red,label="Coop TCR 8, CD4 17")
plot!(t,coop_indx_TC3 .* 100,lw=3,lc=:cyan,label="Coop TCR 3, CD4 22")
xlabel!("Contact Time tc")
ylabel!("Cooperative Index (%)")
