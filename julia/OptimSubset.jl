module OptimSubset

using Random, StatsBase, IterTools, Distributed, ParallelDataTransfer

export max_subset, min_subset, optim_subset_iter, dist_subset, dist_subset_rand

function max_subset(func, length_sub, length_range; keep = 1, show = false)

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

function min_subset(func, length_sub, length_range; keep = 1, show = false)

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

function optim_subset_iter(func, length_sub, length_range; reps = 1, type = "max", keep = 1, max_iter = 100, time = false)

    (range_perm, sub, sub_temp, opt) = osi_alloc_(length_sub, length_range, type == "max", keep)

    temp = Vector{Vector{Tuple{Float64, Vector{Int}}}}(undef, reps)

    for ind in 1:reps
        if time == true
            @time temp[ind] = osi_!(func, range_perm, sub, sub_temp, opt, type == "max", keep, max_iter)
        else
            temp[ind] = osi_!(func, range_perm, sub, sub_temp, opt, type == "max", keep, max_iter)
        end
    end

    return sort!(vcat(temp...); by = x->x[1], rev = type == "max")[1:keep]
end

function osi_alloc_(length_sub, length_range, max_, keep)

    range_perm = collect(1:length_range)
    sub = Vector{Int}(undef, length_sub)
    sub_temp = Vector{Int}(undef, length_sub)
    opt = Vector{Int}(undef, length_range - length_sub)

    return range_perm, sub, sub_temp, opt
end

function osi_!(func, range_perm, sub, sub_temp, opt, max_, keep, max_iter)

    length_sub = length(sub)

    out_list = [(max_ == true ? -Inf : Inf, zeros(Int, length_sub)) for _ in 1:keep]

    shuffle!(range_perm)

    for ind in 1:length_sub
        sub[ind] = range_perm[ind]
    end
    for ind in 1:length(opt)
        opt[ind] = range_perm[length_sub + ind]
    end

    for ind in 1:max_iter
        if osi_inner_!(func, sub, sub_temp, opt, out_list, max_) == true
            break
        end
    end

    keep_out = sum(x[2][1] > 0 for x in out_list)

    return out_list[(keep-keep_out+1):keep]
end

function osi_inner_!(func, sub, sub_temp, opt, out_list, max_)

    shuffle!(sub)
    shuffle!(opt)

    length_sub = length(sub)
    length_opt = length(opt)

    keep = length(out_list)

    for ind in 1:length_sub
        sub_temp[ind] = sub[ind]
    end

    val_ext = out_list[1][1]

    for ind_old in 1:length_sub
        for ind_new in 1:length_opt

            sub_temp[ind_old] = opt[ind_new]

            val = func(sub_temp)

            if max_ == true ? val > val_ext : val < val_ext

                sort!(sub_temp)

                for ind in 1:keep
                    if max_ == true ? val < out_list[ind][1] : val > out_list[ind][1]
                        break
                    elseif val == out_list[ind][1] && sub_temp == out_list[ind][2]
                        for ind in 1:length_sub
                            sub_temp[ind] = sub[ind]
                        end
                        @goto LOOP_END
                    end
                end

                let temp = out_list[1][2]
                    out_list[1] = (val, temp)
                    for ind in 1:length_sub
                        temp[ind] = sub_temp[ind]
                    end
                end

                sort!(out_list;
                    by = x->x[1], rev = !max_)

                val_ext = out_list[1][1]

                let temp = sub[ind_old]
                    sub[ind_old] = opt[ind_new]
                    opt[ind_new] = temp
                end

                return false
            else
                sub_temp[ind_old] = sub[ind_old]
            end
            
            @label LOOP_END
        end
    end

    return true
end

function dist_subset(func, length_sub, length_range; lower, upper, num_bins)

    out_bin = range(lower, upper;
        length = num_bins+1)

    out_dist = zeros(Int, num_bins)

    for sub in IterTools.subsets(1:length_range, length_sub)
        val = func(sub)
        out_dist[floor(Int, (val-lower)/(upper-lower)*num_bins)+1] += val >= lower && val < upper
    end

    return out_dist, out_bin
end

function dist_subset_rand(func, length_sub, length_range; lower, upper, num_bins, reps)

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
