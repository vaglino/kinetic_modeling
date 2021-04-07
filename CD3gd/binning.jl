include("C:/Users/stravaglino3/Downloads/julia_code/julia_code/model_fit_global.jl")
include("C:/Users/stravaglino3/Downloads/julia_code/julia_code/models.jl")
include("C:/Users/stravaglino3/Downloads/julia_code/julia_code/helper_functions.jl")


using DataFrames

CD3_data = CSV.read("C:/Users/stravaglino3/Downloads/julia_code/CD3_data/dissociations.csv",DataFrame)

delete_missing(x) = [collect(skipmissing(x[:,i])) for i in 1:ncol(x)]

ge = delete_missing(CD3_data[:,2:3])
de = delete_missing(CD3_data[:,8:9])
mix = delete_missing(CD3_data[:,14:15])
gethrml = delete_missing(CD3_data[:,5:5])


scatter(ge[1],ge[2])
scatter(de[1],de[2])
scatter(mix[1],mix[2])

# for each experiment: [forces, lifetimes]


function bin_dissociations(data, n_bins; bin_type = "force")
    F = data[1]
    t = data[2]
    # even force bins
    if bin_type == "force"
        f_min = minimum(F)
        f_max = maximum(F)

        edges = range(f_min,stop=f_max,length=n_bins+1)
        # h = fit(Histogram, F)
        bins,bin_ids = partition(F,edges)
        binned_t = [similar(t,0) for i in bins]
        for (i,bin) in enumerate(bin_ids)
            push!(binned_t[bin],t[i])
        end

    end
    return bins, bin_ids,binned_t
end


function partition(data, intervals)
    ranges = intervals[1:end-1] .=> intervals[2:end]
    bins = [similar(data, 0) for _ in 1:length(ranges)]
    ids = similar(data,0)
    for x in data
        for (i, (a, b)) in pairs(ranges)
            if a <= x < b
                push!(bins[i], x)
                push!(ids,i)
                break
            end
        end
    end
    return bins, Int64.(ids)
end

f_binned, bin_ids, binned_t = bin_dissociations(ge,8)
using StatsBase
t_avg = mean.(binned_t)
f_avg = mean.(f_binned)
plot!(f_avg,t_avg)
histogram(binned_t[5])
    # even events number bins

binned_t_sorted = sort.(binned_t)
sort_id = sortperm.(binned_t)

[f_binned[i][sort_id[i]] for i in 1:length(sort_id)]

survival_p(t) = collect((length(t) :-1:1) / length(t))

ps = (survival_p.(binned_t_sorted))

p2lnp(p) = map((x) -> log.(x), p)
# lnps = map((x) -> log.(x), ps)
plot(binned_t_sorted,lnps)


gethrml_sorted = sort(gethrml)
p_gethrml = survival_p.(gethrml_sorted)
lnp_gethrml = p2lnp(p_gethrml)
plot(gethrml_sorted,lnp_gethrml)





using StatsBase, LinearAlgebra
x = randn(1000)*10
h = fit(Histogram, x)
h.edges
h.weights

normalize(h, mode=:pdf)

plot(h)
