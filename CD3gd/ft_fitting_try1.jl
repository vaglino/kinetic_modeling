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



s1=scatter(γϵ[1],γϵ[2], label="γϵ")
s2=scatter(δϵ[1],δϵ[2], label="δϵ")
s3=scatter(mix[1],mix[2], label="mix")
plot(s1,s2,s3)

mutable struct F_t_Data
    F
    t
end

mutable struct M
    model   # model function
    u₀      # initial conditions
    tspan   # time span
    dens    # site densities
    cons    # constant parameters
end
γϵ = F_t_Data(γϵ[1],γϵ[2])

function solve_model(M,p,ts)

    f = (du,u,p,t) -> M.model(du,u,10 .^ p,t,M.dens,M.cons)
    prob = ODEProblem(f,M.u₀,M.tspan,p)
    sol = solve(prob,Rodas4P(),tstops=ts, cb=PositiveDomain(),abstol=1e-7,reltol=1e-7,verbose=false)
end

model = bimolecular_diss_bell_single!
# model = bimolecular_diss!
u₀ = [1.0]
tspan = (0.0,maximum(γϵ.t))

γϵ_ft = M(model,u₀,tspan,[],[])

rates = [0.05,0.5]
# rates = [1.0]

Fs = 0.0:1.:30.0

function ftsolve(model,rates,ts)
    sols = []
    for f in Fs
        model.cons = f
        sol = solve_model(model,rates,ts)
        push!(sols,sol)
    end
    return sols
end


sols = ftsolve(γϵ_ft,rates)
dt = 1
ts = 0:dt:10
p = zeros(length(ts),length(Fs))
pdf = zeros(length(ts),length(Fs))
function construct_surf()
    for (i,f) in enumerate(Fs)
        pi = sum(Array(sols[i](ts))',dims=2)
        p[:,i] = pi
        pdfi = p2pdf(pi,ts)
        pdf[:,i] = pdfi
    end
end
function p2pdf(p,t)
    @show find_mean_lifetime(p,t)
    pdf = p ./ find_mean_lifetime(p,t)
end

function find_mean_lifetime(p,t)
    dt =  t[2]-t[1]
    sum(p)*dt
end


construct_surf()

# pyplot()
gr()
surface(Fs,ts,p)
surface(Fs,ts,pdf,camera=(0,90))


# gr()
# sol = sols[5]
# plot(sol)
# plot!(sol.t,Array(sol)')
#
# sol[7]

model = bimolecular_diss_bell!
# model = bimolecular_diss!
u₀ = [1.0,0.0]
tspan = (0.0,maximum(γϵ.t))

γϵ_ft = M(model,u₀,tspan,[],[])

rates = [11.0,0.1,
        0.14,0.4,
        5.0,0.2,
        20.0,4.0]

sols = ftsolve(γϵ_ft,rates)
dt = 1
ts = 0:dt:10
p = zeros(length(ts),length(Fs))
pdf = zeros(length(ts),length(Fs))

construct_surf()
p
pyplot()
gr()
using ColorSchemes
plot(Fs,ts,p,st=:contourf,fill=(true,cgrad(:inferno)),camera=(0,90))
plot(Fs,ts,pdf,camera=(0,90),st=:contourf,fill=(true,cgrad(:turbo)))

sol = sols[15]
sum(sol,dims=2)*dt
ps = sum(sol(ts),dims=1)
mean_lifetime(ps,ts)
plot(sol)
plot(sum(sol(0:0.1:10)',dims=2))
plot!(sol.t,Array(sol)')
