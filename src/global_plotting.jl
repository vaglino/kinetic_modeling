function plot_sum(rates,Color =nothing)
    prob = ODEProblem(f,u₀,tspan,rates)
    # prob = remake(prob, p=10 .^ rates)
    sol = solve(prob,Rosenbrock23())
    summ = reduce(vcat,sum(sol,dims=1))
    if Color == nothing
        plot!(sol.t,summ)
    else
        plot!(sol.t,summ,c=Color)
    end
end

function display_rates(rates)
    for rate in rates
        @show rate
    end
end
