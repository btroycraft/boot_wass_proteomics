module OptimSubset

import Random, StatsBase, IterTools

export max_subset, min_subset, max_subset_iter, min_subset_iter, dist_subset, dist_subset_rand

max_subset = function(func, length_sub, length_range; keep = 1, show = false)

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

            show == true && @show val_max
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

min_subset = function(func, length_sub, length_range; keep = 1, show = false)

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

            show == true && @show val_min
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

max_subset_iter = function(func, length_sub, length_range; reps = 1, keep = 1, iter = 1000, show = false)

    # Gets estimated maximum by iterative search through subsets of range_arg. Requires a function that takes as input a set of indices. keep gives the number of top elements to return. reps gives the number of randomly selected starting points

    range_perm = vcat(1:length_range)

    opt_size = length_range - length_sub

    opt = Vector{Int}(undef, opt_size)
    sub = Vector{Int}(undef, length_sub)
    sub_temp = Vector{Int}(undef, length_sub)

    val_max_total = val_max = val_max_sub = typemin(Float64)

    ind_max_opt = ind_max_sub = 1

    keep_out = 0

    out_list = Vector{Tuple{Float64, Vector{Int}}}(undef, keep)
    for ind in 1:keep
        out_list[ind] = (val_max_total, zeros(Int, length_sub))
    end

    for _ in 1:reps

        Random.shuffle!(range_perm)

        for ind in 1:length_sub
            sub[ind] = range_perm[ind]
        end
        for ind in 1:opt_size
            opt[ind] = range_perm[ind + length_sub]
        end

        val_max = func(sub)

        for _ in 1:iter

            val_max_sub = typemin(Float64)

            for ind_old in 1:length_sub

                ind_swap = sub[ind_old]

                for ind_new in 1:opt_size

                    sub[ind_old] = opt[ind_new]

                    val = func(sub)

                    if val > val_max_sub
                        val_max_sub = val
                        ind_max_sub = ind_old
                        ind_max_opt = ind_new
                    end

                    if val > val_max_total

                        sub_temp .= sub
                        sort!(sub_temp)

                        for ind in 1:keep
                            if sub_temp == out_list[ind][2]
                                @goto LOOP_EXIT
                            end
                        end

                        keep_out += keep_out < keep
                        out_list[1] = (val, out_list[1][2])
                        out_list[1][2] .= sub_temp

                        sort!(out_list;
                            by = x->x[1])

                        val_max_total = out_list[1][1]

                        show == true && @show val_max_total

                        @label LOOP_EXIT
                    end
                end

                sub[ind_old] = ind_swap
            end

            val_max_sub <= val_max && break

            val_max = val_max_sub

            (sub[ind_max_sub], opt[ind_max_opt]) = (opt[ind_max_opt], sub[ind_max_sub])
        end
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

min_subset_iter = function(func, length_sub, length_range; reps = 1, keep = 1, iter = 1000, show = false)

    # Gets estimated minimum by iterative search through subsets of range_arg. Requires a function that takes as input a set of indices. keep gives the number of top elements to return. reps gives the number of randomly selected starting points

    range_perm = vcat(1:length_range)

    opt_size = length_range - length_sub

    opt = Vector{Int}(undef, opt_size)
    sub = Vector{Int}(undef, length_sub)
    sub_temp = Vector{Int}(undef, length_sub)

    val_min_total = val_min = val_min_sub = Inf

    ind_min_opt = ind_min_sub = 1

    keep_out = 0

    out_list = Vector{Tuple{Float64, Vector{Int}}}(undef, keep)
    for ind in 1:keep
        out_list[ind] = (val_min_total, Vector{Int}(undef, length_sub))
    end

    for _ in 1:reps

        Random.shuffle!(range_perm)

        for ind in 1:length_sub
            sub[ind] = range_perm[ind]
        end
        for ind in 1:opt_size
            opt[ind] = range_perm[ind + length_sub]
        end

        val_min = func(sub)

        for _ in 1:iter

            val_min_sub = -Inf

            for ind_old in 1:length_sub

                ind_swap = sub[ind_old]

                for ind_new in 1:opt_size

                    sub[ind_old] = opt[ind_new]

                    val = func(sub)

                    if val < val_min_sub
                        val_min_sub = val
                        ind_min_sub = ind_old
                        ind_min_opt = ind_new
                    end

                    if val < val_min_total

                        sub_temp .= sub
                        sort!(sub_temp)

                        for ind in 1:keep
                            if sub_temp == out_list[ind][2]
                                @goto LOOP_EXIT
                            end
                        end

                        keep_out += keep_out < keep
                        out_list[1] = (val, out_list[1][2])
                        out_list[1][2] .= sub_temp

                        sort!(out_list;
                            by = x->x[1],
                            rev = true)

                        val_min_total = out_list[1][1]

                        show == true && @show val_min_total

                        @label LOOP_EXIT
                    end
                end

                sub[ind_old] = ind_swap
            end

            val_min_sub >= val_min && break

            val_min = val_min_sub

            (sub[ind_min_sub], opt[ind_min_opt]) = (opt[ind_min_opt], sub[ind_min_sub])
        end
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
