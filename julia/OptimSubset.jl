module OptimSubset

export max_subset, min_subset, max_subset_iter, min_subset_iter, dist_subset, dist_subset_rand

max_subset = function(func, length_sub, range_arg; keep = 1)

    # Gets global maximum by exhaustive search through subsets of range_arg. Requires a function that takes as input a set of indices. keep gives the number of top elements to return

    out_type = typeof(func(range_arg[1:length_sub]))
    val_max = typemin(out_type)
    keep_out = 0

    out_list = Vector{Tuple{out_type, Vector{Int}}}(undef, keep)
    for ind in 1:keep
        out_list[ind] = (val_max, Vector{Int}(undef, length_sub))
    end

    for sub in IterTools.subsets(range_arg, length_sub)

        val = func(sub)

        if val > val_max

            keep_out += keep_out < keep

            out_list[1] = (val, out_list[1][2])
            out_list[1][2] = sub

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

    out_vec = Vector{out_type}(undef, keep_out)
    for ind in 1:keep_out
        out_vec[ind] = out_list[keep-ind+1][1]
    end

    return out_vec, out_mat
end

min_subset = function(func, length_sub, range_arg; keep = 1)

    # Gets global minimum by exhaustive search through subsets of range_arg. Requires a function that takes as input a set of indices. keep gives the number of bottom elements to return

    out_type = typeof(func(range_arg[1:length_sub]))
    val_min = typemax(out_type)
    keep_out = 0

    out_list = Vector{Tuple{out_type, Vector{Int}}}(undef, keep)
    for ind in 1:keep
        out_list[ind] = (val_min, Vector{Int}(undef, length_sub))
    end

    for sub in IterTools.subsets(range_arg, length_sub)

        val = func(sub)

        if val < val_min

            keep_out += keep_out < keep

            out_list[1] = (val, out_list[1][2])
            out_list[1][2] = sub

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

    out_vec = Vector{out_type}(undef, keep_out)
    for ind in 1:keep_out
        out_vec[ind] = out_list[keep-ind+1][1]
    end

    return out_vec, out_mat
end

max_subset_iter = function(func, length_sub, range_arg; reps = 1, keep = 1, iter = 1000)

    # Gets estimated maximum by iterative search through subsets of range_arg. Requires a function that takes as input a set of indices. keep gives the number of top elements to return. reps gives the number of randomly selected starting points

    length_range = length(range_arg)
    opt_size = length_range - length_sub

    out_type = typeof(func(range_arg[1:length_sub]))

    range_perm = vcat(range_arg)
    opt = Vector{Int}(opt_size)
    sub = Vector{Int}(length_sub)
    sub_temp = Vector{Int}(length_sub)

    val_max_total = val_max = val_max_sub = typemin(out_type)

    ind_max_opt = ind_max_sub = 1

    keep_out = 0

    out_list = Vector{Tuple{out_type, Vector{Int}}}(undef, keep)
    for ind in 1:keep
        out_list[ind] = (val_max_total, Vector{Int}(undef, length_sub))
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

            val_max_sub = typemin(out_type)

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
                            val > out_list[ind][1] && break
                            if val == out_list[ind][1] && sub_temp == out_list[ind][2]
                                @goto LOOP_EXIT
                            end
                        end

                        keep_out += keep_out < keep
                        out_list[1] = (val, out_list[1][2])
                        out_list[1][2] .= sub_temp

                        sort!(out_list;
                            by = x->x[1])

                        val_min_total = out_list[1][1]

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

    out_vec = Vector{out_type}(undef, keep_out)
    for ind in 1:keep_out
        out_vec[ind] = out_list[keep-ind+1][1]
    end

    return out_vec, out_mat
end

min_subset_iter = function(func, length_sub, range_arg; reps = 1, keep = 1, iter = 1000)

    # Gets estimated minimum by iterative search through subsets of range_arg. Requires a function that takes as input a set of indices. keep gives the number of top elements to return. reps gives the number of randomly selected starting points

    length_range = length(range_arg)
    opt_size = length_range - length_sub

    out_type = typeof(func(range_arg[1:length_sub]))

    range_perm = vcat(range_arg)
    opt = Vector{Int}(opt_size)
    sub = Vector{Int}(length_sub)
    sub_temp = Vector{Int}(length_sub)

    val_min_total = val_min = val_min_sub = typemax(out_type)

    ind_min_opt = ind_min_sub = 1

    keep_out = 0

    out_list = Vector{Tuple{out_type, Vector{Int}}}(undef, keep)
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

            val_min_sub = typemin(out_type)

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
                            val < out_list[ind][1] && break
                            if val == out_list[ind][1] && sub_temp == out_list[ind][2]
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

    out_vec = Vector{out_type}(undef, keep_out)
    for ind in 1:keep_out
        out_vec[ind] = out_list[keep-ind+1][1]
    end

    return out_vec, out_mat
end

dist_subset = function(func, length_sub, range_arg; bound, length_bin)

    out_type = float(typeof(func(range_arg[1:length_subset])))

    lower = convert(out_type, bound[1])
    upper = convert(out_type, bound[2])

    out_bin = range_arg(lower, upper;
        length = length_bin+1)

    out_dist = zeros(Int, length_bin)

    for sub in IterTools.subsets(range_arg, length_sub)
        val = func(sub)
        out_dist[floor(Int, (val-lower)/(upper-lower))] += val >= lower && val < upper
    end

    return out_dist, out_bin
end

dist_subset_rand = function(func, length_sub, range_arg; bound, length_bin, reps)

    out_type = float(typeof(func(range_arg[1:length_subset])))

    lower = convert(out_type, bound[1])
    upper = convert(out_type, bound[2])

    out_bin = range_arg(lower, upper;
        length = length_bin+1)

    out_dist = zeros(Int, length_bin)

    sub = Vector{Int}(undef, length_sub)

    for _ in 1:reps

        StatsBase.sample!(range_arg, sub)
        val = func(sub)

        out_dist[floor(Int, (val-lower)/(upper-lower))] += val >= lower && val < upper
    end

    return out_dist, out_bin
end

end
