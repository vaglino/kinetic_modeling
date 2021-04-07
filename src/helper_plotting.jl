function solve_model(M,rates,dens,cons)
    f = (du,u,p,t) -> M.model(du,u,p,t,dens,cons)
    prob = ODEProblem(f,M.u₀,M.tspan,rates)
    sol = solve(prob,AutoVern7(KenCarp4()), abstol=1e-8,reltol=1e-8)
end


function plot_diss_fit!(plt,M,rates,data;cmap=nothing)

    for i in 1:length(M.cons)
        # Setup the ODE problem, then solve
        sol = solve_model(M,rates,M.dens[i],M.cons[i])

        Σ = sum(Array(sol),dims=1)
        Σ = Im2Re(Σ)
        Σ = normalize_lnp_fit(Σ,data.y[i])
        # h = plot!(sol.t,real.(Σ)',ylims=(-5,0))

        if cmap==nothing
            color = i
        else
            color = cmap
        end
        plot!(plt,sol.t,Σ,
            ylims=(0.0,1.0),
            # xlims=(0,15),
            color=color,
            label=string(data.f[i]))
    end
    if cmap==nothing
        colors = collect(1:length(M.cons))'
    else
        colors = cmap
    end
    plot_dissociations!(plt,data,colors)
    # return plt
end

function Im2Re(S)  # ln of the solutions (gets rid of integration instabilities)
    Σ = Complex.(S)
    Σ = real.(log.(Σ)')
end

function normalize_lnp(lnp)
    normlnpv = lnp ./ log(length(lnp)) .+ 1.0
end

function normalize_lnp_fit(lnp_fit,lnp_data)
    normlnpv_fit = lnp_fit ./ log(length(lnp_data)) .+ 1.0
end


function plot_dissociations!(plt,data,colors)
    normlnp = normalize_lnp.(p2lnp.(data.y))
    ε=0.01
    scatter!(plt,data.t,normlnp,
            xlims = (0.0-ε, maximum(maximum.(data.t))+ε),
            ylims = (0.0-ε, 1.0+ε),
            markercolors = colors, #, markerstrokecolor = colors,
            # markersize = 5, markershape=:square,
            # markeralpha = 0., markerstrokealpha = 1.,
            label=false)
end

function show_diss_fit!(plt,M,rates,data;cmap=nothing)

    plot_diss_fit!(plt,M,rates,data;cmap)
    # return plt
end
# function Im2Re(S)
#     Σ = Complex.(S)
#     Σ = real.(Σ)
#     Σ = log.(Σ')
# end
#


# function plot_diss_fit(rates,dens,cons)
#     for (i,F) in enumerate(cons)
#         # Setup the ODE problem, then solve
#         solve_model(rates,dens,cons)
#
#         Σ = Complex.(sum(Array(sol),dims=1))
#         # h = plot!(sol.t,real.(Σ)',ylims=(-5,0))
#         h = plot!(sol.t,real.(log.(Σ)'),
#             ylims=(-6.,0),
#             xlims=(0,15),
#             color=i,
#             label=string(F),
#             palette = :Dark2_5)
#
#         if generate_data
#             push!(dataset,vec(generate_data(sol;logplot=true)))
#         end
#
#         # display(h)
#     end
#     return h
# end
#
