using OrdinaryDiffEq
using ModelingToolkit
using DataDrivenDiffEq
using LinearAlgebra, DiffEqSensitivity, Optim
using DiffEqFlux, Flux
using Plots
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
    du[3] = dn₃ = z[1] + z[3] - z[2] - z[4] + z[5]

end
function triNN!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ,nₓ,n₃ = u
    k₁,k₋₁,k₂,k₋₂ = cons
    mᵣ,mₗ,mₓ,A = dens
    z = L(u,p)
    # model
    du[1] = dnᵣ = A*k₁*mᵣ*mₗ - k₋₁*nᵣ + z[1] - z[2]
    du[2] = dnₓ = A*k₂*mₓ*mₗ - k₋₂*nₓ + z[3] - z[4]
    du[3] = dn₃ = z[5] - z[6]

end

L = FastChain(FastDense(3, 50, tanh),
                FastDense(50, 100, tanh),
                FastDense(100, 100, tanh),
                FastDense(100, 50, tanh),
                FastDense(50, 6))
#L = FastChain(FastDense(3, 20, tanh),FastDense(20, 20, tanh), FastDense(20, 2))
p = initial_params(L)

# Define the experimental parameter
tspan = (0.0,16.1) # A little bit longer because there are two values to save at final t?
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

f = (du,u,p,t) -> triNN!(du,u, p,t,densities_tc1,cons_tc)
prob_nn = ODEProblem(f,u0, tspan, p)
sol_nn = solve(prob_nn, Tsit5(), saveat = t)

# plot(solution)
plot(sol_nn)
# summ = reduce(vcat,sum(sol_nn,dims=1))
# h = plot!(sol_nn.t,summ,linecolor=:black)

function predict(θ)

    tmp_prob = remake(prob_nn,u0=u0,p=θ)
    tmp_sol =  solve(tmp_prob, VCABM(), saveat = t,
                  abstol=1e-5, reltol=1e-5)
                  # backsolve=true)
    Σ_sol = sum(Array(tmp_sol),dims=1) # Note: this returns a row vector!
end

function loss(θ)
    pred = predict(θ)
    # sum(abs2, (Xₙ .- pred') .* t), pred
    sum(abs2, (Xₙ .- pred')), pred
end

# Test
@time loss(p)

losses = []

callback(θ,l,pred) = begin
    push!(losses, l)
    @show l
    if length(losses)%50==0
        println("Current loss after $(length(losses)) iterations: $(losses[end])")
    end
    p = plot(t, pred')
    scatter!(p, t, Xₙ)
    display(p)
    false
end

# First train with ADAM for better convergence
@time res1 = DiffEqFlux.sciml_train(loss, p, ADAM(0.005), cb=callback, maxiters = 10)

# Train with BFGS
res2 = DiffEqFlux.sciml_train(loss, res1.minimizer, BFGS(initial_stepnorm=0.01),
                                cb=callback, maxiters = 1000)

res2 = DiffEqFlux.sciml_train(loss, res1.minimizer, NelderMead(),
                                cb=callback)

println("Final training loss after $(length(losses)) iterations: $(losses[end])")

# Plot the losses
plot(losses, yaxis = :log, xaxis = :log, xlabel = "Iterations", ylabel = "Loss")

# Plot the data and the approximation
NNsolution = predict(res2.minimizer)
# Trained on noisy data vs real solution
plot(t, NNsolution')
scatter!(t, Xₙ)
