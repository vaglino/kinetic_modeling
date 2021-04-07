
delete_missing(x) = [collect(skipmissing(x[:,i])) for i in 1:ncol(x)]

function bin_dissociations(data, n_bins; bin_type = "force", edges=nothing)
    F = data[1]
    t = data[2]
    # even force bins
    if edges != nothing
        bins,bin_ids = partition(F,edges)
        binned_t = [similar(t,0) for i in bins]

    elseif bin_type == "force"
        f_min = minimum(F)
        f_max = maximum(F)

        edges = range(f_min,stop=f_max,length=n_bins+1)
        bins,bin_ids = partition(F,edges)
        binned_t = [similar(t,0) for i in bins]

    elseif bin_type == "equal number"
        edges_ind = Int64.(round.(range(1,stop=length(F),length=n_bins+1)))
        edges = F[edges_ind]
        # edges = [0,9,12,15,19,30]
        bins,bin_ids = partition(F,edges)
        binned_t = [similar(t,0) for i in bins]
    end
    # h = fit(Histogram, F)

    for (i,bin) in enumerate(bin_ids)
        push!(binned_t[bin],t[i])
    end

    return bins, binned_t, bin_ids
end


function partition(data, intervals)
    ranges = intervals[1:end-1] .=> intervals[2:end]
    bins = [similar(data, 0) for _ in 1:length(ranges)]
    ids = similar(data,0)
    for x in data
        for (i, (a, b)) in pairs(ranges)
            if i != length(bins) && a <= x < b
                push!(bins[i], x)
                push!(ids,i)
                break

            elseif i == length(bins) && a <= x <= b
                push!(bins[i], x)
                push!(ids,i)
                break
        end
        end
    end
    return bins, Int64.(ids)
end

function resort(f,t)
    t_sorted = sort.(t)
    sort_id = sortperm.(t)
    f_sorted = [f[i][sort_id[i]] for i in 1:length(sort_id)]
    return f_sorted, t_sorted
end

survival_p(t) = collect((length(t) :-1:1) / length(t))

p2lnp(p) = map((x) -> log.(x), p)
