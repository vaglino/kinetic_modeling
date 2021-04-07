include("C:/Users/stravaglino3/Downloads/julia_code/julia_code/model_fit_global.jl")
include("C:/Users/stravaglino3/Downloads/julia_code/models.jl")
include("C:/Users/stravaglino3/Downloads/julia_code/helper_functions.jl")

f = (du,u,p,t) -> bimolecular_diss!(du,u,p,t,dens,cons) # enclose constants

# u0 = [0.0;0.0]
u₀ = [1.0]          # Initial condition
tspan = (0.0,10.0)  # Simulation interval
p = 0.5      # equation parameter. p = [k₁, k₋₁]
dens = []
cons = []

# Setup the ODE problem, then solve
prob = ODEProblem(f,u₀,tspan,p)
sol = solve(prob,Tsit5())
# sol = solve(prob,Tsit5())
plot(sol)
dataset1 = generate_data(sol)
# plot(sol.t,log.(Array(sol))')
#
# fs = collect(0:20)
# ks = bell_diss.(0.5,10,fs)
# plot(fs,ks)
model_opt2 = (1.0e-6,10.0)
pinit = 0.0
t = collect(range(0,stop=10,length=101))
rates,hes,l = optimization(bimolecular_diss!,[vec(dataset1)],
                            [[]],[[]],model_opt2,[t])

Juno.@enter optimization(bimolecular_diss!,[vec(dataset1)],
            [[]],[[]],model_opt2,[t])
