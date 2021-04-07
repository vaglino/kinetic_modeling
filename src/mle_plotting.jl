using ColorSchemes

# solve model for given rates and force interval
function ftsolve(M,rates,Fs,ts)
    sols = []
    for f in Fs
        # M.cons = f
        M.f = f
        sol = solve_prob(M,rates;tᵢ=ts)
        push!(sols,sol)
    end
    return sols
end

# take model solution and construct P and PDF surfaces for plotting
function construct_surf(sols,Fs,ts)
    p = []
    pdf = []
    for (i,f) in enumerate(Fs)

        pᵢ = sum(Array(sols[i](ts))',dims=2)
        pdfᵢ = p2pdf(pᵢ,sols[i])

        p = push!(p,pᵢ)
        pdf = push!(pdf,pdfᵢ)
    end
    p = reduce(hcat,p)
    pdf = reduce(hcat,pdf)
    return p, pdf
end

# compute pdf from p distribution
function p2pdf(p,sol)
    mean_lifetime(sol)
    pdf = p ./ mean_lifetime(sol)
end


function plot_ft_fit(M,rates,data;dt=0.1,dF=1.0)

    ts = 0:dt:maximum(data.t)
    Fs = 0:dF:maximum(data.F)
    # Fs = 0:dF:20
    sols = ftsolve(M,rates,Fs,ts)

    p, pdf = construct_surf(sols,Fs,ts)

    # plot contours

    p1 = plot(Fs,ts,log.(p),st=:contour,
                            fill=(true),#cgrad(:plasma)),
                            camera=(0,90),
                            clims=(-7,0),
                            ylabel="ln(p)",xlabel="F")
    scatter!(data.F,data.t,c=:white,markeralpha=0.7)
    display(p1)

    pex = surface(Fs,ts,(p),#,st=:contour,
                            fill=(true),#cgrad(:plasma)),
                            # camera=(0,90),
                            # clims=(-7,0),
                            zlabel="p",xlabel="F")
    display(pex)

    p2 = plot(Fs,ts,p,st=:contour,
                        fill=(true),#,cgrad(:plasma)),
                        camera=(0,90),
                        ylabel="p",xlabel="F")
    display(p2)

    p3 = plot(Fs,ts,log.(pdf),st=:contour,
                        fill=(true,cgrad(:plasma)),
                        camera=(0,90),
                        ylabel="pdf",xlabel="F")

    # plot lnp curves

    p4=plot()
    lnp = p2lnp.(p)
    for i in 1:1:size(lnp,2)
            plot!(ts,lnp[:,i],c = i,
                                xlims=(0,5),
                                ylims=(-7,0),
                                legend=false,
                                ylabel="ln(p)",xlabel="t")
    end


    # plot mean lifetime distribution
    t_means = mean_lifetime.(sols)
    p5 = plot(Fs,t_means,ylabel="<t>",xlabel="F")

    display(plot(p1,p2,p3,p4,p5))
    # display(p2)
    # display(p3)
    # display(p4)
    # display(p5)
    return sols
end


function plot_mean_ft_fit(M,k,data,binned;species=[],dt=0.1,dF=1.0)

    ts = 0:dt:maximum(data.t)
    Fs = 0:dF:maximum(data.F)

    M.solver = Tsit5  #change solver to Tsit5 (can handle ± std)
    sols = ftsolve(M,k,Fs,ts)
    mean_t = mean_lifetime.(sols)
    h = plot!(Fs,mean_lifetime.(sols),c=Color,
                    # ribbon=0.1,
                    markerstrokecolor=:auto,
                    linealpha = 0.3,
                    ylabel="<t>",
                    xlabel="F",
                    label=species)
    plot_F_vs_t(clean_dissociation_data(binned))
end

function extrapolate_fit(M,k,data,F_max;species=[],dt=0.1,dF=1.0)

    ts = 0:dt:maximum(data.t)
    Fs = 0:dF:maximum(F_max)

    # M.solver = Tsit5  #change solver to Tsit5 (can handle ± std)
    sols = ftsolve(M,k,Fs,ts)
    mean_t = mean_lifetime.(sols)
    h = plot!(Fs,log10.(mean_lifetime.(sols)),c=Color,
                    # ribbon=0.1,
                    markerstrokecolor=:auto,
                    linealpha = 1,
                    ylabel="<t>",
                    xlabel="F",
                    label=species)
end


function mean_ft_fit(M,k,data,F_max;species=[],dt=0.1,dF=1.0)

    ts = 0:dt:maximum(data.t)
    Fs = 0:dF:F_max

    M.solver = Tsit5  #change solver to Tsit5 (can handle ± std)
    sols = ftsolve(M,k,Fs,ts)
    mean_t = mean_lifetime.(sols)
end
