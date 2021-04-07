n1 = 1. # bimolecular 1
n2 = 2. # bimolecular 2
n3 = 3. # trimolecular
function cooperativity_index(n1,n2,n3)
    Σ_bi = n1 + n2
    Σ_tri = n3
    Δn = Σ_tri - Σ_bi
    # cooperativity index formulav (relative change)
    Coop_indx = Δn / Σ_bi
end

cooperativity_index(n1,n2,n3)

n1 = [1. ± 0.1, 1. ± 0.1]  # bimolecular 1
n2 = [2. ± 0.1, 2. ± 0.1]# bimolecular 2
n3 = [3. ± 0.1, 4. ± 0.1]

cooperativity_index.(n1,n2,n3)



function solve_with_sds(rates_sds)
    prob = ODEProblem(f,u₀,tspan,10 .^ rates_sds)
    # prob = remake(prob, p=10 .^ rates)
    sol = solve(prob,Vern7(),saveat=0.01, abs_tol=1e-10,rel_tol=1e-10)
    n_tot = reduce(vcat,sum(sol,dims=1))
    return n_tot
end
function solve_model(rates_sds)
    prob = ODEProblem(f,u₀,tspan,10 .^ rates_sds)
    # prob = remake(prob, p=10 .^ rates)
    sol = solve(prob,Vern7(),saveat=0.01, abs_tol=1e-10,rel_tol=1e-10)
    return sol
end

#-------------------------------------------------------------------------------

densities_tcr = [25.0, 38.0, 1.0] # TCR, pMHC, Ac
densities_tc1 = [15.0, 38.0, 10.0, 1.0] # TCR, pMHC, CD4, Ac
densities_tc2 = [8.0,  38.0, 17.0, 1.0]
densities_tc3 = [3.0,  38.0, 22.0, 1.0]

# rates_T = [-3.37638358011457 ± 0.0483562400944489
#             -0.302013187183543 ± 0.0607189115063286]
# rates_TC1 = [-1.72689171836867 ± 0.469108114325548
#             -0.646695212785203 ± 0.685834531823618]
# rates_TC2 = [-0.644293816490510 ± 0.408735939559103
#             0.617796010276307 ± 0.467542299161307]
# rates_TC3 = [0.891017383250648 ± 0.559684851912230
#             2.16934160248305 ± 0.568940078989975]

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

tspan = (0.0,16.0)
t = collect(0.01:0.01:16)

gr()
# bimolecular TCR-pMHC bonds
u₀ = [0.0]
f = (du,u,p,t) -> bimolecular!(du,u,p,t,densities_tcr,[])
n_tcr = Array(solve_model(rates_T))'[2:end]
n_tcr_sd = n_tcr .± nT_sds

# TCR-pMHC-CD4 (1)
u₀ = [0.0, 0.0]
f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_tc1,10 .^ rates_tcr)
n_tc1 = Array(solve_model(rates_TC1))'[2:end,:]
n_tc1_sd = sum(n_tc1 .± nTC1_sds, dims=2)


# TCR-pMHC-CD4 (2)
f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_tc2,10 .^ rates_tcr)
n_tc2 = Array(solve_model(rates_TC2))'[2:end,:]
n_tc2_sd = sum(n_tc2 .± nTC2_sds, dims=2)
plot!(t,n_tc2_sd, errorbars=:ribbon, alpha=.1,c=:red)

# TCR-pMHC-CD4 (3)
f = (du,u,p,t) -> trimolecular!(du,u,p,t,densities_tc3,10 .^ rates_tcr)
n_tc3 = Array(solve_model(rates_TC3))'[2:end,:]
n_tc3_sd = sum(n_tc3 .± nTC3_sds, dims=2)
plot!(t,n_tc3_sd, errorbars=:ribbon, alpha=0.1,ribboncolor=:cyan)


plot(t, n_tcr_sd, errorbars=:ribbon, alpha=0.1,msc=:black)
plot!(t, n_tcr, lw=3, lc=:black)
plot!(t,n_tc1_sd, errorbars=:ribbon, alpha=0.1, msc=:blue)
plot!(t, sum(n_tc1,dims=2), lw=3, lc=:blue)
plot!(t,n_tc2_sd, errorbars=:ribbon, alpha=0.1, msc=:red)
plot!(t, sum(n_tc2,dims=2), lw=3, lc=:red)
plot!(t,n_tc3_sd, errorbars=:ribbon, alpha=0.1, msc=:cyan)
plot!(t, sum(n_tc3,dims=2), lw=3, lc=:cyan)

scatter!(n_se_tcr[:,1],n_se_tcr[:,2],yerror=n_se_tcr[:,3],markercolor=:black,markersize=5)
scatter!(n_se_1[:,1],n_se_1[:,2],yerror=n_se_1[:,3],markercolor=:blue,markersize=5)
scatter!(n_se_2[:,1],n_se_2[:,2],yerror=n_se_2[:,3],markercolor=:red,markersize=5)
h = scatter!(n_se_3[:,1],n_se_3[:,2],yerror=n_se_3[:,3],markercolor=:cyan,markersize=5,legend=false)

xlabel!("Contact Time tc")
ylabel!("<n>")
# savefig("all_n_with_error_bands")

#-------------------------------------------------------------------------------
# FIT CD4 only data


cd4_1 = CSV.read("cd4_data_from_muaz.csv")
t_cd41, Pa_cd41, n_se_cd41 = clean_data(cd4_1)
n_cd41 = Pa2n(Pa_cd41)

f = (du,u,p,t) -> bimolecular!(du,u,p,t,dens,cons) # enclose constants

# u0 = [0.0;0.0]
u₀ = [0.0]          # Initial condition
tspan = (0.0,8.0)  # Simulation interval
pinit = [1e-4,0.5]      # equation parameter. p = [k₁, k₋₁]
densities_cd41 = [1700., 850., 1.] #CD4, pMHC, A
cons = []

rates_cd41,hes_cd41,l_cd41,boot_cd41 = optimization(bimolecular!,n_cd41,densities_cd41',[],model_opt2,t_cd41)

# Setup the ODE problem, then solve
plotly()
t = collect(0.01:0.01:16)
tspan = (0.0,16.0)
f = (du,u,p,t) -> bimolecular!(du,u,p,t,densities_cd41,cons)
prob = ODEProblem(f,u₀,tspan,10 .^ rates_cd41)
sol_cd41 = Array(solve(prob,Vern7(),saveat=t,abstol=1e-10,reltol=1e-10))'
plot(t,sol_cd41,lc=:orange,lw=3)
scatter!(n_se_cd41[:,1],n_se_cd41[:,2],yerror=n_se_cd41[:,3],markercolor=:orange)


densities_tcr = [25.0, 38.0, 1.0] # TCR, pMHC, Ac
densities_tc1 = [15.0, 38.0, 10.0, 1.0] # TCR, pMHC, CD4, Ac
densities_tc2 = [8.0,  38.0, 17.0, 1.0]
densities_tc3 = [3.0,  38.0, 22.0, 1.0]
densities_cd41 = [1700., 850., 1.] #CD4, pMHC, A
cooperativity_index(n1,n2,n3)
#---------------------
#normalize everything to tc2

n_tcr_norm = n_tcr./ 25 .* 8
n_tc1_norm = sum(n_tc1,dims=2) ./15 .* 8 ./10 .*17
n_tc2_norm = sum(n_tc2,dims=2)
n_tc3_norm = sum(n_tc3,dims=2) ./3 .* 8 ./22 .* 17
n_cd41_norm = sol_cd41 ./1700 .*17 ./850 .*38
n_tcr_norm = n_tcr ./25. .*15.
n_tc1_norm = sum(n_tc1,dims=2)
n_tc2_norm = sum(n_tc2,dims=2) ./8. .*15. ./17. .*10.
n_tc3_norm = sum(n_tc3,dims=2) ./3. .*15. ./22. .*10.
n_cd41_norm = sol_cd41 ./1700. .*10. ./850. .*38

plotly()
t = collect(0.01:0.01:16)
plot(t,n_tcr_norm,lw=3,lc=:black)
plot!(t,n_tc1_norm,lw=3,lc=:blue)
plot!(t,n_tc2_norm,lw=3,lc=:red)
plot!(t,n_tc3_norm,lw=3,lc=:cyan)
plot!(t,n_cd41_norm,lw=3,lc=:orange,yscale=:linear,legend=false)
yaxis!("log<n>")
xlabel!("Contact Time tc")
savefig("n_log_scaled")

coop_indx_tc1 = cooperativity_index.(n_tcr_norm,n_cd41_norm,n_tc1_norm)
coop_indx_tc2 = cooperativity_index.(n_tcr_norm,n_cd41_norm,n_tc2_norm)
coop_indx_tc3 = cooperativity_index.(n_tcr_norm,n_cd41_norm,n_tc3_norm)
plot(t,coop_indx_tc1 .* 100,lw=3,lc=:blue)
plot!(t,coop_indx_tc2 .* 100,lw=3,lc=:red)
plot!(t,coop_indx_tc3 .* 100,lw=3,lc=:cyan,legend=false)
xlabel!("Contact Time tc")
ylabel!("Cooperative Index (%)")
# savefig("Cooperative Index.png")
