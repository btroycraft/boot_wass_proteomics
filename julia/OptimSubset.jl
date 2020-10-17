module OptimSubset

import Random, StatsBase, IterTools

export max_subset, min_subset, max_subset_iter, min_subset_iter, dist_subset, dist_subset_rand

max_subset = function(func, length_sub, length_range; keep = 1)

    # Gets global maximum by exhaustive search through subsets of range_arg. Requires a function that takes as input a set of indices. keep gives the number of top elements to return

    val_max = typemin(Float64)
    keep_out = 0

    out_list = Vector{Tuple{Float64, Vector{Int}}}(undef, keep)
    for ind in 1:keep
        out_list[ind] = (val_max, Vector{Int}(undef, length_sub))
    end

    for sub in IterTools.subsets(1:length_range, length_sub)

        val = func(sub)

        if val > val_max

            keep_out += keep_out < keep

            out_list[1] = (val, out_list[1][2])
            out_list[1][2] .= sub

            sort!(out_list;
                by = x->x[1])

            val_max = out_list[1][1]
        end
    end

    for ind in (keep-keep_out+1):keep
        sort!(out_list[ind][2])
    end

    out_mat = Matrix{Int}(undef, length_sub, keep_out)
    for ind in 1:keep_out
        out_mat[:, ind] = out_list[keep-ind+1][2]
    end

    out_vec = Vector{Float64}(undef, keep_out)
    for ind in 1:keep_out
        out_vec[ind] = out_list[keep-ind+1][1]
    end

    return out_vec, out_mat
end

min_subset = function(func, length_sub, length_range; keep = 1)

    # Gets global minimum by exhaustive search through subsets of range_arg. Requires a function that takes as input a set of indices. keep gives the number of bottom elements to return

    val_min = typemax(Float64)
    keep_out = 0

    out_list = Vector{Tuple{Float64, Vector{Int}}}(undef, keep)
    for ind in 1:keep
        out_list[ind] = (val_min, Vector{Int}(undef, length_sub))
    end

    for sub in IterTools.subsets(1:length_range, length_sub)

        val = func(sub)

        if val < val_min

            keep_out += keep_out < keep

            out_list[1] = (val, out_list[1][2])
            out_list[1][2] .= sub

            sort!(out_list;
                by = x->x[1],
                rev = true)

            val_min = out_list[1][1]
        end
    end

    for ind in (keep-keep_out+1):keep
        sort!(out_list[ind][2])
    end

    out_mat = Matrix{Int}(undef, length_sub, keep_out)
    for ind in 1:keep_out
        out_mat[:, ind] = out_list[keep-ind+1][2]
    end

    out_vec = Vector{Float64}(undef, keep_out)
    for ind in 1:keep_out
        out_vec[ind] = out_list[keep-ind+1][1]
    end

    return out_vec, out_mat
end

optim_subset_iter = function(func, length_sub, length_range; reps = 1, keep = 1, max_iter = 1000, type = "max", worker_ids = workers())

    @everywhere worker_ids begin
        optim_subset_iter_ = optim_subset_iter_close($func, $length_sub, $length_range, $keep, $max_iter, $type)
    end

    keep_list = pmap(optim_subset_iter_, worker_ids, 1:reps)

    keep_out = min(keep, sum(map(length, keep_list)))
    out_list = vcat(keep_out...)[1:keep_out]

    return out_list
end

optim_subset_iter_close = function(func, length_sub, length_range, keep, max_iter, type)

    # Gets estimated maximum by fiterative search through subsets of range_arg. Requires a function that takes as input a set of indices. keep gives the number of top elements to return. reps gives the number of randomly selected starting points

    length_opt = length_range - length_sub
    range_perm = collect(1:length_range)
    sub = Vector{Int}(undef, length_sub)
    opt = Vector{Int}(undef, length_range - length_sub)

    if type == "max"
        return function()

            out_list = Vector{Tuple{Float64, Vector{Int}}}(undef, 0)

            Random.shuffle!(range_perm)

            for ind in 1:length_sub
                sub[ind] = range_perm[ind]
            end
            for ind in 1:length_opt
                opt[ind] = range_perm[length_sub + ind]
            end

            for _ in 1:max_iter
                osi_inner_max_(func, sub, opt, out_list. keep)
            end

            return out_list
        end
    elseif type == "min"
        return function()

            out_list = Vector{Tuple{Float64, Vector{Int}}}(undef, 0)

            Random.shuffle!(range_perm)

            for ind in 1:length_sub
                sub[ind] = range_perm[ind]
            end
            for ind in 1:length_opt
                opt[ind] = range_perm[length_sub + ind]
            end

            for _ in 1:max_iter
                osi_inner_min_(func, sub, opt, out_list. keep)
            end

            return out_list
        end
    end
end

osi_inner_max_ = function(func, sub, opt, out_list, keep)

    Random.shuffle!(sub)
    Random.shuffle!(opt)

    sub_temp = copy(sub)

    for ind_old in 1:length(sub)
        for ind_new in 1:length(opt)

            sub_temp[ind_old] = opt[ind_new]

            val = func(sub_temp)

            if length(out_list) == 0 || val > out_list[1][1]

                sort!(sub_temp)

                for ind in 1:length(out_list)
                    if sub_temp == out_list[ind][2]
                        @goto INNER_LOOP_END
                    end
                end

                if length(out_list) < keep
                    push!(out_list, (val, sub_temp))
                else
                    out_list[1] = (val, sub_temp)
                end

                sort!(out_list;
                    by = x->x[1])

                opt[ind_new] = sub[ind_old]
                sub[ind_old] = sub_temp[ind_old]

                return
            end

            @label INNER_LOOP_END
        end

        sub_temp[ind_old] = sub[ind_old]
    end

    return
end

osi_inner_min_ = function(func, sub, opt, out_list, keep)

    Random.shuffle!(sub)
    Random.shuffle!(opt)

    sub_temp = copy(sub)

    for ind_old in 1:length(sub)
        for ind_new in 1:length(opt)

            sub_temp[ind_old] = opt[ind_new]

            val = func(sub_temp)

            if length(out_list) == 0 || val < out_list[1][1]

                sort!(sub_temp)

                for ind in 1:length(out_list)
                    if sub_temp == out_list[ind][2]
                        @goto INNER_LOOP_END
                    end
                end

                if length(out_list) < keep
                    push!(out_list, (val, sub_temp))
                else
                    out_list[1] = (val, sub_temp)
                end

                sort!(out_list;
                    by = x->x[1],
                    rev = true)

                opt[ind_new] = sub[ind_old]
                sub[ind_old] = sub_temp[ind_old]

                return
            end

            @label INNER_LOOP_END
        end

        sub_temp[ind_old] = sub[ind_old]
    end

    return
end

dist_subset = function(func, length_sub, length_range; lower, upper, num_bins)

    out_bin = range(lower, upper;
        length = num_bins+1)

    out_dist = zeros(Int, num_bins)

    for sub in IterTools.subsets(1:length_range, length_sub)
        val = func(sub)
        out_dist[floor(Int, (val-lower)/(upper-lower)*num_bins)+1] += val >= lower && val < upper
    end

    return out_dist, out_bin
end

dist_subset_rand = function(func, length_sub, length_range; lower, upper, num_bins, reps)

    out_bin = range(lower, upper;
        length = num_bins+1)

    out_dist = zeros(Int, num_bins)

    sub = Vector{Int}(undef, length_sub)

    for ind in 1:reps

        StatsBase.sample!(1:length_range, sub)
        val = func(sub)
        val >= lower && val < upper && (out_dist[floor(Int, (val-lower)/(upper-lower)*num_bins)+1] += 1)
    end

    return out_dist, out_bin
end



end
