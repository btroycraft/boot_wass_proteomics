module OptimSubset

using Random, StatsBase, IterTools, Distributed, ParallelDataTransfer

export optim_subset_iter, dist_subset, dist_subset_rand

function optim_subset_iter(func, length_sub, length_range; reps = 1, type = "max", keep = 1, max_iter = 10^2, worker_ids = workers())

    if length(worker_ids) == 1 && worker_ids[1] == 1

        out_list = vcat(map(1:reps) do _
            optim_subset_iter_!(func, length_sub, length_range, type == "max", keep, max_iter)
        end...)
    else

        for id in worker_ids
            sendto(id, func = func, length_sub = length_sub, length_range = length_range, max_ = type == "max", keep = keep, max_iter = max_iter)
        end

        out_list = vcat(pmap(CachingPool(worker_ids), 1:reps) do _
            optim_subset_iter_!(func, length_sub, length_range, type == "max", keep, max_iter)
        end...)
    end

    sort!(out_list;
        by = x -> x[2])
    for ind in length(out_list):-1:2
        if(out_list[ind][2] == out_list[ind-1][2])
            out_list[ind] = (type == "max" ? -Inf : Inf, out_list[ind][2])
        end
    end
    sort!(out_list;
        by = x -> x[1],
        rev = type == "max")
    
    out_list = out_list[1:keep]
    keep = sum(map(out_list) do x
        type == "max" ? x[1] > -Inf : x[1] < Inf
    end)
    out_list = out_list[1:keep]
    
    return out_list
end

function optim_subset_iter_!(func, length_sub, length_range, max_, keep, max_iter)

    range_perm = collect(1:length_range)
    sub = Vector{Int}(undef, length_sub)
    sub_temp = Vector{Int}(undef, length_sub)
    opt = Vector{Int}(undef, length_range - length_sub)

    out_list = [(max_ == true ? -Inf : Inf, zeros(Int, length_sub)) for _ in 1:keep]

    shuffle!(range_perm)

    for ind in 1:length_sub
        sub[ind] = range_perm[ind]
    end
    for ind in 1:length(opt)
        opt[ind] = range_perm[length_sub + ind]
    end

    for ind in 1:max_iter
        if optim_subset_iter_inner_!(func, sub, sub_temp, opt, out_list, max_) == true
            break
        end
    end

    keep_out = sum(x[2][1] > 0 for x in out_list)

    return out_list[(keep-keep_out+1):keep]
end

function optim_subset_iter_inner_!(func, sub, sub_temp, opt, out_list, max_)

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

                val = func(sort!(sub_temp))

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
                    by = x->x[1],
                    rev = !max_)

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
