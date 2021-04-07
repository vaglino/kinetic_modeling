using DifferentialEquations
using ModelingToolkit
using DataDrivenDiffEq
using LinearAlgebra, DiffEqSensitivity, Optim
using DiffEqFlux, Flux
using Plots
gr()

# function triNN!(du,u,p,t,dens,cons)
#     # unpack rates and constants
#     nᵣ,nₓ,n₃ = u
#     k₁,k₋₁,k₂,k₋₂ = cons
#     mᵣ,mₗ,mₓ,A = dens
#     z = L(u,p)
#     # model
#     du[1] = dnᵣ = A*k₁*mᵣ*mₗ - k₋₁*nᵣ + z[1]
#     du[2] = dnₓ = A*k₂*mₓ*mₗ - k₋₂*nₓ + z[2]
#     du[3] = dn₃ = z[1] + z[2]
#
# end
function bimolecular!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ = u[1]
    k₁,k₋₁  = p
    mᵣ,mₗ,A = dens
    # model
    du[1] = dnᵣ = A*k₁*mᵣ*mₗ - k₋₁*nᵣ

end

dens = [100.,10.,1.]
cons = Float64[]
f = (du,u,p,t) -> bimolecular!(du,u,p,t,dens,cons)
# Define the experimental parameter
tspan = (0.0f0,10.0f0)
u0 = Float32[0.]
p_ = Float32[1e-4, 0.5]
prob = ODEProblem(f, u0,tspan, p_)
solution = solve(prob, Tsit5(), saveat = 0.1)

scatter(solution, alpha = 0.25)
plot!(solution, alpha = 0.5)

# Ideal data
X = Array(solution)
# Add noise to the data
println("Generate noisy data")
Xₙ = X + Float32(1e-2)*randn(eltype(X), size(X))

# function dudt_(u, p,t)
#     x = u
#     z = L(u,p)
#     [p_[1]*x + z[1],
#     -p_[4]*y + z[2]]
# end
function biNN!(du,u,p,t,dens,cons)
    # unpack rates and constants
    nᵣ = u[1]
    k₁,k₋₁  = p
    mᵣ,mₗ,A = dens
    z = L(u,p)
    # model
    du[1] = dnᵣ = z[1]

end

L = FastChain(FastDense(1, 32, tanh),FastDense(32, 32, tanh), FastDense(32, 1))
#L = FastChain(FastDense(3, 20, tanh),FastDense(20, 20, tanh), FastDense(20, 2))
p = initial_params(L)

#
fNN = (du,u,p,t) -> biNN!(du,u, p,t,dens,cons)
prob_nn = ODEProblem(fNN,u0, tspan, p)
sol_nn = solve(prob_nn, Tsit5(), saveat =  solution.t)

plot(solution)
plot!(sol_nn)
# summ = reduce(vcat,sum(sol_nn,dims=1))
# h = plot!(sol_nn.t,summ,linecolor=:black)

function predict(θ)

    tmp_prob = remake(prob_nn,u0=u0,p=θ)
    tmp_sol =  solve(tmp_prob, VCABM(), saveat = solution.t,
                  abstol=1e-5, reltol=1e-5)
                  # backsolve=true)
    # Σ_sol = sum(Array(tmp_sol),dims=1) # Note: this returns a row vector!
    Array(tmp_sol)
end

function loss(θ)
    pred = predict(θ)
    sum(abs2, Xₙ .- pred), pred
end

# Test
@time loss(p)

const losses = []

callback(θ,l,pred) = begin
    push!(losses, l)
    @show l
    if length(losses)%50==0
        println("Current loss after $(length(losses)) iterations: $(losses[end])")
    end
    h = plot(solution.t, pred')
    scatter!(h, solution.t, Xₙ')
    display(h)
    false
end

p = initial_params(L)
# First train with ADAM for better convergence
@time res1 = DiffEqFlux.sciml_train(loss, p, ADAM(0.01), cb=callback, maxiters = 40)

# Train with BFGS
res2 = DiffEqFlux.sciml_train(loss, res1.minimizer, BFGS(),
                                cb=callback, maxiters = 100)

# @time res1 = DiffEqFlux.sciml_train(loss, res2.minimizer, ADAM(0.001), cb=callback, maxiters = 20)

println("Final training loss after $(length(losses)) iterations: $(losses[end])")

# Plot the losses
imgloss= plot(losses, yaxis = :log, xaxis = :log, xlabel = "Iterations", ylabel = "Loss",fmt = :svg)
savefig(imgloss,"training_loss")
# Plot the data and the approximation
NNsolution = predict(res2.minimizer)
# Trained on noisy data vs real solution
plot(solution.t, NNsolution',linewidth=3,xlabel = "t (s)", ylabel = "n",fmt = :svg)
scatter!(solution.t, Xₙ')
savefig("fit")


L̂ = L(Xₙ,res2.minimizer)
## Sparse Identification
#
# # Create a Basis
# @variables u[1:2]
# # Lots of polynomials
# polys = Operation[1]
#
# for i ∈ 1:2
#     push!(polys, u[1]^i)
#     push!(polys, u[2]^i)
#     for j ∈ i:2
#         if i != j
#             push!(polys, (u[1]^i)*(u[2]^j))
#             push!(polys, u[2]^i*u[1]^i)
#         end
#     end
# end
#
# # And some other stuff
# h = [cos.(u)...; sin.(u)...; polys...]
# h = [ polys...]
# basis = Basis(h, u)
#
# # Create an optimizer for the SINDy problem
# opt = SR3()
# # Create the thresholds which should be used in the search process
# λ = exp10.(-7:0.1:3)
# # Target function to choose the results from; x = L0 of coefficients and L2-Error of the model
# g(x) = x[1] < 1 ? Inf : norm(x, 2)
# # Test on original data and without further knowledge
# # println("SINDy on full ideal, unavailable data")
# # Ψ = SINDy(Xₙ[:, :], DX[:, :], basis, λ, opt, g = g, maxiter = 10000) # Fail
# # println(Ψ)
# # print_equations(Ψ)
# #
# # # Test on ideal derivative data for unknown function ( not available )
# # println("SINDy on partial ideal, unavailable data")
# # Ψ = SINDy(Xₙ[:, 1:end], L̄[:, 1:end], basis, λ,opt, g = g, maxiter = 10000) # Succeed
# # println(Ψ)
# # print_equations(Ψ)
#
# # Test on uode derivative data
# println("SINDy on learned, partial, available data")
# Ψ = SINDy(Xₙ[:, 2:end], L̂[:, 2:end], basis, λ,  opt, g = g, maxiter = 100000, normalize = true, denoise = true) # Succeed
# println(Ψ)
# print_equations(Ψ)
#
# # Extract the parameter
# p̂ = parameters(Ψ)
# println("First parameter guess : $(p̂)")

@variables u[1:1]
h = Operation[u; u.^2; 1]
basis = Basis(h, u)

# Create the thresholds which should be used in the search process
λ = exp10.(-7:0.1:3)
# Target function to choose the results from; x = L0 of coefficients and L2-Error of the model
g(x) = x[1] < 1 ? Inf : norm(x, 2)

opt = SR3(3e-1, 1.0)
# Ψ = SInDy(X[:, 1:1000], DX[:, 1:1000], basis, maxiter = 10000, opt = opt, normalize = true)
Ψ = SINDy(Xₙ[:, 2:end], L̂[:, 2:end], basis, λ,  opt, g = g, maxiter = 10000, normalize = true) # Succeed
# Ψ = SINDy(Xₙ[:, 2:end], L̂[:, 2:end], basis, maxiter = 10000, opt = opt, normalize = true,denoise = true)
println(Ψ)
print_equations(Ψ)

sys = ODESystem(Ψ)
p = parameters(Ψ)

dudt = ODEFunction(sys)

estimator = ODEProblem(dudt, u0, tspan, p)
estimation = solve(estimator, Tsit5(), saveat = solution.t)

plot!(estimation,linecolor=:red,linewidth=3)
savefig("sindy_raw_estim")

using DiffEqParamEstim
using Optim
# cost_function = build_loss_objective(estimator,t,Xₙ,Tsit5(),maxiters=10000)

cost_function = build_loss_objective(estimator,Tsit5(),L2Loss(solution.t,Xₙ),
                                     maxiters=10000)#,verbose=false)
result = optimize(cost_function, p, BFGS())
tmp_estimator = remake(estimator, p=result.minimizer)
better_estimation = solve(tmp_estimator, Tsit5(), saveat = solution.t)

plot!(better_estimation,linecolor=:black,linestyle=:dash,linewidth=4)
savefig("sindy_better_estim")
#
# # The parameters are a bit off, but the equations are recovered
# # Start another SINDy run to get closer to the ground truth
# # Create function
# unknown_sys = ODESystem(Ψ)
# unknown_eq = ODEFunction(unknown_sys)
#
# # Just the equations
# b = Basis((u, p, t)->unknown_eq(u, [1.; 1.], t), u)
#
# # Retune for better parameters -> we could also use DiffEqFlux or other parameter estimation tools here.
# Ψf = SINDy(Xₙ[:, 2:end], L̂[:, 2:end], b, STRRidge(0.01), maxiter = 100, convergence_error = 1e-18) # Succeed
# println(Ψf)
# p̂ = parameters(Ψf)
# println("Second parameter guess : $(p̂)")
#
# # Create function
# recovered_sys = ODESystem(Ψf)
# recovered_eq = ODEFunction(recovered_sys)
# estimator = ODEProblem(dudt, u0, tspan, p)
# estimation = solve(estimator, Tsit5(), saveat = solution.t)
#
# plot!(estimation)
#
#
#
