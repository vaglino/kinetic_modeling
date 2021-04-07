
#--------------------------helper functions------------------------------------
using Statistics
function clean_data(data)

    t = data[:,1]
    Pa1 = data[:,2]
    Pa2 = data[:,3]

    Pa_se = std(hcat(Pa1,Pa2),dims=2)
    n1 = Pa2n(Pa1)
    n2 = Pa2n(Pa2)
    n_mean = mean(hcat(n1,n2),dims=2)
    n_se = std(hcat(n1,n2),dims=2)

    tstacked = vcat(t,t)
    Pa = vcat(Pa1,Pa2)

    p = sortperm(tstacked)
    tstacked = tstacked[p]
    Pa = Pa[p]
    return tstacked, Pa, hcat(t,n_mean,n_se)
end


function generate_data(sol)
  t = collect(range(0,stop=sol.t[end],length=101))

  solution = sum(sol(t),dims=1)
  # randomized = VectorOfArray([(sol(t[i]) + .02randn(size(sol,1))) for i in 1:length(t)])
  # dataset = convert(Array,randomized)
  randomized = solution .+ .01randn(size(solution))
  dataset = randomized
  display(scatter!(t,dataset'))
  return dataset
end

function generate_data(sol,t)

  solution = sum(sol(t),dims=1)
  # randomized = VectorOfArray([(sol(t[i]) + .02randn(size(sol,1))) for i in 1:length(t)])
  # dataset = convert(Array,randomized)
  randomized = solution .+ .01randn(size(solution))
  dataset = randomized
  display(scatter!(t,dataset'))
  return dataset
end


function n2Pa(n) # convert adhesion frequency to <n>
    Pa = 1 .- exp.(-n)
end

function Pa2n(Pa) # convert <n> to adhesion frequency
    n = log.(1 ./ (1 .- Pa))
end

using CSV
function clean_dissociation_data(data)
    n = size(data,2) ÷ 4 # 4 columns per force bin
    clean_data = []
    for i in 0:n-1
        ind_bin = i*4+1:i*4+4
        data_bin = data[:,ind_bin]
        mask_non_zeros = data_bin[:,1] .!= 0
        if sum(mask_non_zeros) > 1
             # clean_data[i+1] = data_bin[mask_non_zeros,:]
             push!(clean_data,data_bin[mask_non_zeros,:])
        end
    end
    return clean_data
end

function separate_diss_data(data)
    t = [bin[:,2] for bin in data]
    lnp = [bin[:,4] for bin in data]
    return t, lnp
end


function rescale_t(data)
    rescaled = data .- minimum(data)
end

function plot_dissociations(data)
    lb = 0
    for i in 1:length(data)
        F = mean(data[i][:,1])
        h = scatter!(Array(data[i][:,2]),Array(data[i][:,4]),
                        color=i,
                        markeralpha=0.5,
                        markerstrokealpha=1,
                        markersize=5,
                        label = string(F),
                        palette = :Dark2_5)
        if data[i][end,4] < lb
            lb = data[i][end,4]
        end
        # ylims! = (-10,0.)
        # display(h)
    end
end

function plot_F_vs_t(data)
    D = Ft_average_and_error(data)
    h = scatter!(D.f_mean,D.t_mean,
                    yerror=D.t_sem,
                    xerror=D.f_sem,
                    c=Color,
                    label="",
                    markerstrokecolor=:auto,
                    markersize=7)
    display(h)
end

function plot_F_vs_t(f,t)
    D = Ft_average_and_error(f,t)
    h = scatter!(D.f_mean,D.t_mean,
                    yerror=D.t_sem,
                    xerror=D.f_sem,
                    c=Color,
                    label="",
                    markerstrokecolor=:auto,
                    markersize=7)
end

function Ft_average_and_error(f,t)
    t_mean = []
    f_mean = []
    t_sem = []
    f_sem = []
    for (fi,ti) = zip(f,t)
        n = length(fi[:,1])
        push!(t_mean, mean(ti))
        push!(f_mean, mean(fi))
        push!(t_sem, std(ti)/sqrt(n) )
        push!(f_sem, std(fi)/sqrt(n) )
    end
    df = DataFrame(f_mean=f_mean,f_sem=f_sem,
                    t_mean=t_mean,t_sem=t_sem)
end
function Ft_average_and_error(data)
    t_mean = []
    f_mean = []
    t_sem = []
    f_sem = []
    for bin in data
        n = length(bin[:,1])
        push!(t_mean, mean(bin[:,2]))
        push!(f_mean, mean(bin[:,1]))
        push!(t_sem, std( bin[:,2])/sqrt(n) )
        push!(f_sem, std(bin[:,1])/sqrt(n) )
    end
    df = DataFrame(f_mean=f_mean,f_sem=f_sem,
                    t_mean=t_mean,t_sem=t_sem)
end


using Statistics
function extract_avg_forces(data)
    avg_fs = [mean(bin[:,1]) for bin in data]
end


function extract_max_lifetimes(data)
    max_t = [maximum(bin[:,2]) for bin in data]
end


# using Interpolations, Loess, ImageFiltering
#
# function max_lifetime_bounds(f,t,xs)
#     max_f = maximum([xs[end], f[end]])
#     fs = vcat(0., f, max_f)
#     ts = vcat(t[1], t, t[end])
#     # fs = vcat(f, max_f)
#     # ts = vcat(t, t[end])
#     interp_linear = LinearInterpolation(fs, ts)
#     ys = interp_linear(xs)
#     # h = plot!(xs,ys)
#
#     ker = ImageFiltering.Kernel.gaussian((3,))
#     smoothed_ys = imfilter(ys, ker)
#     # smoothed_ys = moving_average(ys,20)
#     # model = loess(collect(xs), ys)
#     # smoothed_ys = predict(model, xs)
#     # h = plot!(xs,smoothed_ys)
#     # display(h)
#     return smoothed_ys
# end

function lnp2p(lnp)
    p = exp.(lnp)
end

function plot_bell_equations(p)
    f = collect(0.0:1:35)
    plt = []
    if length(p) == 8
        y1 = log10.(bell_diss_kT.(p[1],p[2],f))
        y2 = log10.(bell_diss_kT.(p[3],p[4],f))
        y3 = log10.(bell_diss_kT.(p[5],p[6],f))
        y4 = log10.(bell_ass_kT.(p[7],p[8],f))
        # ytot = y1+y2
        rates = [y1,y2,y3,y4]
        plt = plot(f,rates,
                    ylims=(minimum(minimum.(rates)),maximum(maximum.(rates))),
                    labels=["k₋f" "k₋s" "ka" "k₋a"])
    elseif length(p) == 4
        y1 = bell_ass_kT.(p[1],p[2],f)
        y2 = bell_diss_kT.(p[3],p[4],f)
        # ytot = y1+y2
        plt = plot(f,[y1,y2],ylims=(0,30))
    end
    plt
end


function plot_diss_fit(model,Fs; generate_data = false)
    for (i,F) in enumerate(Fs)
        # Setup the ODE problem, then solve
        f = (du,u,p,t) -> model(du,u,p,t,dens,F)
        prob = ODEProblem(f,u₀,tspan,p)
        sol = solve(prob,AutoVern7(KenCarp4()),abstol=1e-8,reltol=1e-8)

        Σ = Complex.(sum(Array(sol),dims=1))
        # h = plot!(sol.t,real.(Σ)',ylims=(-5,0))
        h = plot!(sol.t,real.(log.(Σ)'),
            ylims=(-5.5,0),
            xlims=(0,6),
            color=i,
            label=string(F),
            palette = :Dark2_5)

        if generate_data
            push!(dataset,vec(generate_data(sol;logplot=true)))
        end

        # display(h)
    end
end

function show_diss_fit(model, data, forces)
  h = plot()
  plot_diss_fit(model,forces)
  plot_dissociations(data)
  return h
end
function show_diss_fit(model, t, lnp, forces)
    h = plot()
    plot_diss_fit(model,forces)
    scatter!(t, lnp)
    return h
end

using NumericalIntegration
function plot_force_lifetime_fit(model,Fs,data)

    mean_f  = extract_avg_forces(data)
    max_t  = extract_max_lifetimes(data)
    # mean_f = mean.(f)
    # max_t = maximum.(t)
    # scatter(mean_f,max_t)
    t_bounds = max_lifetime_bounds(mean_f,max_t,Fs)

    avg_lifetimes = zeros(length(Fs))
    for (i,F) in enumerate(Fs)
        # Setup the ODE problem, then solve
        f = (du,u,p,t) -> model(du,u,p,t,dens,F)
        prob = ODEProblem(f,u₀,tspan,p)
        sol = solve(prob,AutoVern7(KenCarp4()),abstol=1e-8,reltol=1e-8)

        t_query = collect(range(0.,stop=t_bounds[i],length=100))
        sol_at_bound_t = sol(t_query)

        Σ = real.(Complex.(sum(Array(sol_at_bound_t),dims=1)))
        # Σ_sol = real.(Complex.(sum(sol_at_bound_t)))
        # h = plot!(sol.t,real.(Σ)',ylims=(-5,0))
        # Σ = vec(Σ_sol)
        # t = sol.t[sol.t .<= t_bounds[i]]
        # Σ = Σ[sol.t .<= t_bounds[i]]
        avg_lifetime = integrate(t_query,vec(Σ))
        avg_lifetimes[i] = avg_lifetime
    end
    plot!(Fs,avg_lifetimes)
                # ylims=(0,Inf),
                # xlims=(0,Inf))
end


# given vector of lifetimes > compute lnp
function calculate_lnp(t)
    lnp = log.( (length(t):-1:1) ./ length(t) )
end


function measurements2Df(t,y)
    y_value = Measurements.value.(y)
    y_error = Measurements.uncertainty.(y)
    df = DataFrame(x = t, y = y_value, error = y_error)
end


function combine_bins(data)
    F = []
    t = []
    for bin in data
        F = vcat(F,bin[:,1])
        t = vcat(t,bin[:,2])
    end
    perm_id = sortperm(F)
    F = F[perm_id]
    t = t[perm_id]
    return Float64.(F), Float64.(t)
end
