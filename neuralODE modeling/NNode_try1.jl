using DifferentialEquations
using ModelingToolkit
# using DataDrivenDiffEq
using LinearAlgebra, DiffEqSensitivity, Optim
using DiffEqFlux, Flux
using Plots, Revise
gr()


function triNN!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ,nₓ,n₃ = u
    k₁,k₋₁,k₂,k₋₂ = cons
    mᵣ,mₗ,mₓ,A = dens
    z = L(u,p)
    # model
    du[1] = dnᵣ = A*k₁*mᵣ*mₗ - k₋₁*nᵣ - z[1] + z[2]
    du[2] = dnₓ = A*k₂*mₓ*mₗ - k₋₂*nₓ - z[3] + z[4]
    du[3] = dn₃ = z[1] + z[3] - z[2] - z[4]

end
function triNN!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ,nₓ,n₃ = u
    k₁,k₋₁,k₂,k₋₂ = cons
    mᵣ,mₗ,mₓ,A = dens
    z =  L(u,p)
    # model
    du[1] = dnᵣ = A*k₁*mᵣ*mₗ - k₋₁*nᵣ - z[1] + z[2]
    du[2] = dnₓ = A*k₂*mₓ*mₗ - k₋₂*nₓ - z[3] + z[4]
    du[3] = dn₃ = z[1] + z[3] - z[2] - z[4]

end
function triNN!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ,nₓ,n₃ = u
    k₁,k₋₁,k₂,k₋₂ = cons
    mᵣ,mₗ,mₓ,A = dens
    z = L(u,p)
    # model
    du[1] = dnᵣ = A*k₁*mᵣ*mₗ - k₋₁*nᵣ + z[1]
    du[2] = dnₓ = A*k₂*mₓ*mₗ - k₋₂*nₓ + z[2]
    du[3] = dn₃ = z[1] + z[2]

end

L = FastChain(FastDense(3, 32, tanh),FastDense(32, 32, tanh), FastDense(32, 2))
L = FastChain(FastDense(3, 32, tanh),FastDense(32, 50, tanh),FastDense(50, 32, tanh), FastDense(32, 3))

# Define the experimental parameter
tspan = (0.0,16.1)
u0 = Float32[0.0,0.0,0.0]

# Xₙ = n_tc1 #load data
Xₙ = [0.02020270731751947
 0.0
 0.06187540371808745
 0.10536051565782635
 0.083381608939051
 0.19845093872383823
 0.3011050927839216
 0.3011050927839216
 0.3566749439387324
 0.38566248081198473
 0.5108256237659907
 0.6931471805599453
 0.8209805520698303
 0.7339691750802005
 0.8209805520698303
 0.7339691750802005]
 t = [0.25
  0.25
  0.5
  0.5
  1.0
  1.0
  2.0
  2.0
  4.0
  4.0
  6.0
  6.0
  8.0
  8.0
 16.0
 16.0]
densities_tc1 = [15.0,38.0,10.0,1.0]
cons_tc = [ 0.0006692541890287495,0.8662361534770547,1.169965568192585e-6,0.4]



function predict(θ)

    tmp_prob = remake(prob_nn,u0=u0,p=θ)
    tmp_sol =  solve(tmp_prob, VCABM(), saveat = t,
                  abstol=1e-5, reltol=1e-5)
                  # sensealg=BacksolveAdjoint())
                  # backsolve=true)

    Σ_sol = sum(Array(tmp_sol),dims=1)
end

function loss(θ)
    pred = predict(θ)
    # L = sum(abs2, (Xₙ .- pred').*t), pred
    L = sum(abs2, Xₙ .- pred'), pred
    # prediction = pred
    # return L, prediction
end

callback(θ,l,pred) = begin
    push!(losses, l)
    @show l
    if length(losses)%10==0
        println("Current loss after $(length(losses)) iterations: $(losses[end])")
        fig = plot(t,pred')
        scatter!(fig,t,Xₙ)
        display(fig)
    end
    false
end

# function run_model(learning_rate)
p = initial_params(L)
f = (du,u,p,t) -> triNN!(du,u, p,t,densities_tc1,cons_tc)
prob_nn = ODEProblem(f,u0, tspan, p)
sol_nn = solve(prob_nn, Rodas4(), saveat = t)

# plot(solution)
plot(sol_nn)
summ = reduce(vcat,sum(sol_nn,dims=1))
h = plot!(sol_nn.t,summ,linecolor=:black)
# Test
@time loss(p)
# L,predition = loss(p)

# p = initial_params(L)
losses = []
# First train with ADAM for better convergence
res1 = DiffEqFlux.sciml_train(loss, p, ADAM(0.03), cb=callback, maxiters = 300)

# Train with BFGS
res2 = DiffEqFlux.sciml_train(loss, res1.minimizer, BFGS(initial_stepnorm=0.001),
                                cb=callback, maxiters = 100)

println("Final training loss after $(length(losses)) iterations: $(losses[end])")

# Plot the losses
plot(losses, yaxis = :log, xaxis = :log, xlabel = "Iterations", ylabel = "Loss")
end

# run_model(0.005)
# Plot the data and the approximation
NNsolution = predict(res2.minimizer)
# Trained on noisy data vs real solution
plot(t, NNsolution')
scatter!(t, Xₙ)

tmp_prob = remake(prob_nn,u0=u0,p=res2.minimizer)
tmp_sol =  solve(tmp_prob, Tsit5(), saveat = 0.1,
              abstol=1e-5, reltol=1e-5,
              sensealg=BacksolveAdjoint())
              # backsolve=true)

Σ_sol = sum(Array(tmp_sol),dims=1)

plot(tmp_sol.t,Σ_sol')
