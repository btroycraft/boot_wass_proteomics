module OptimSubset

using IterTools: subsets
using Random: shuffle!

export max_subset, min_subset, max_subset_iter, min_subset_iter

max_subset = function(func; subset_size, index_range, keep = 1)

    # Gets global maximum by exhaustive search through subsets of index_range. Requires a function that takes as input a set of indices. keep gives the number of top elements to return

    inds_total = vcat(index_range)
    val_max = typemin(func(inds_total[1:subset_size]))

    num_out = zero(keep)

    max_list = Vector{Tuple{typeof(val_max), Vector{eltype(inds_total)}}}(undef, keep)
    for ind in 1:keep
        max_list[ind] = (val_max, zeros(eltype(inds_total), subset_size))
    end

    for inds_subset in subsets(inds_total, subset_size)

        val = func(inds_subset)

        if val > val_max

            num_out += num_out < keep

            max_list[1] = (val, max_list[1][2])
            max_list[1][2] .= inds_subset

            sort!(max_list, by = x->x[1])

            val_max = max_list[1][1]
        end
    end

    for ind in (keep-num_out+1):keep
        sort!(max_list[ind][2])
    end

    out_mat = Matrix{eltype(inds_total)}(undef, subset_size, num_out)
    for ind in 1:num_out
        out_mat[:, ind] = max_list[keep-ind+1][2]
    end

    out_vec = Vector{typeof(val_max)}(undef, num_out)
    for ind in 1:num_out
        out_vec[ind] = max_list[keep-ind+1][1]
    end

    return out_vec, out_mat
end

min_subset = function(func; subset_size, index_range, keep = 1)

    # Gets global minimum by exhaustive search through subsets of index_range. Requires a function that takes as input a set of indices. keep gives the number of bottom elements to return

    inds_total = vcat(index_range)
    val_min = typemin(func(inds_total[1:subset_size]))

    num_out = zero(eltype(inds_total))

    min_list = Vector{Tuple{typeof(val_min), Vector{eltype(inds_total)}}}(undef, keep)
    for ind in 1:keep
        min_list[ind] = (val_min, zeros(eltype(inds_total), subset_size))
    end

    for inds_subset in subsets(inds_total, subset_size)

        val = func(inds_subset)

        if val > val_min

            num_out += num_out < keep

            min_list[1] = (val, min_list[1][2])
            min_list[1][2] .= inds_subset

            sort!(min_list, by = x->x[1])

            val_min = min_list[1][1]
        end
    end

    for ind in (keep-num_out+1):keep
        sort!(min_list[ind][2])
    end

    out_mat = Matrix{eltype(inds_total)}(undef, subset_size, num_out)
    for ind in 1:num_out
        out_mat[:, ind] = min_list[keep-ind+1][2]
    end

    out_vec = Vector{typeof(val_min)}(undef, num_out)
    for ind in 1:num_out
        out_vec[ind] = min_list[keep-ind+1][1]
    end

    return out_vec, out_mat
end


max_subset_iter = function(func; subset_size, index_range, n_reps = 1, keep = 1, max_iters = 1000)

    # Gets estimated maximum by iterative search through subsets of index_range. Requires a function that takes as input a set of indices. keep gives the number of top elements to return. n_reps gives the number of randomly selected starting points

    range_size = length(index_range)
    option_size = range_size - subset_size

    inds_total = vcat(index_range)
    inds_option = similar(inds_total, option_size)
    inds_subset = similar(inds_total, subset_size)
    inds_subset_temp = similar(inds_total, subset_size)

    val_max_total = val_max = val_max_subset = typemin(func(inds_total[1:subset_size]))

    ind_max_option = ind_max_subset = one(eltype(inds_total))

    num_out = zero(eltype(inds_total))

    max_list = Vector{Tuple{typeof(val_max_total), Vector{eltype(inds_total)}}}(undef, keep)
    for ind in 1:keep
        max_list[ind] = (val_max_total, zeros(eltype(inds_total), subset_size))
    end

    for _ in 1:n_reps

        shuffle!(inds_total)

        for ind in 1:subset_size
            inds_subset[ind] = inds_total[ind]
        end
        for ind in 1:option_size
            inds_option[ind] = inds_total[ind + subset_size]
        end

        val_max = func(inds_subset)

        for _ in 1:max_iters

            val_max_subset = typemin(val_max_total)

            for ind_old in 1:subset_size

                ind_swap = inds_subset[ind_old]

                for ind_new in 1:option_size

                    inds_subset[ind_old] = inds_option[ind_new]

                    val = func(inds_subset)

                    if val > val_max_subset
                        val_max_subset = val
                        ind_max_subset = ind_old
                        ind_max_option = ind_new
                    end

                    if val > val_max_total

                        inds_subset_temp .= inds_subset
                        sort!(inds_subset_temp)

                        flag = true
                        for ind in 1:keep
                            val > max_list[ind][1] && break
                            if val == max_list[ind][1] && inds_subset_temp == max_list[ind][2]
                                flag = false
                                break
                            end
                        end

                        if flag == true
                            num_out < keep && (num_out += 1)
                            max_list[1] = (val, max_list[1][2])
                            max_list[1][2] .= inds_subset_temp

                            sort!(max_list, by = x->x[1])

                            val_max_total = max_list[1][1]
                        end
                    end
                end

                inds_subset[ind_old] = ind_swap
            end

            val_max_subset <= val_max && break

            val_max = val_max_subset

            (inds_subset[ind_max_subset], inds_option[ind_max_option]) = (inds_option[ind_max_option], inds_subset[ind_max_subset])
        end
    end

    out_mat = Matrix{eltype(inds_total)}(undef, subset_size, num_out)
    for ind in 1:num_out
        out_mat[:, ind] = max_list[keep-ind+1][2]
    end

    out_vec = Vector{typeof(val_max_total)}(undef, num_out)
    for ind in 1:num_out
        out_vec[ind] = max_list[keep-ind+1][1]
    end

    return out_vec, out_mat
end

min_subset_iter = function(func; subset_size, index_range, n_reps = 1, keep = 1, max_iters = 1000)

    # Gets estimated minimum by iterative search through subsets of index_range. Requires a function that takes as input a set of indices. keep gives the number of bottom elements to return. n_reps gives the number of randomly selected starting points

    range_size = length(index_range)
    option_size = range_size - subset_size

    inds_total = vcat(index_range)
    inds_option = similar(inds_total, option_size)
    inds_subset = similar(inds_total, subset_size)
    inds_subset_temp = similar(inds_total, subset_size)

    val_min_total = val_min = val_min_subset = typemax(func(inds_total[1:subset_size]))

    ind_min_option = ind_min_subset = one(eltype(inds_total))

    num_out = zero(eltype(inds_total))

    min_list = Vector{Tuple{typeof(val_min_total), Vector{eltype(inds_total)}}}(undef, keep)
    for ind in 1:keep
        min_list[ind] = (val_min_total, zeros(eltype(inds_total), subset_size))
    end

    for _ in 1:n_reps

        shuffle!(inds_total)

        for ind in 1:subset_size
            inds_subset[ind] = inds_total[ind]
        end
        for ind in 1:option_size
            inds_option[ind] = inds_total[ind + subset_size]
        end

        val_min = func(inds_subset)

        for _ in 1:max_iters

            val_min_subset = typemax(val_min_total)

            for ind_old in 1:subset_size

                ind_swap = inds_subset[ind_old]

                for ind_new in 1:option_size

                    inds_subset[ind_old] = inds_option[ind_new]

                    val = func(inds_subset)

                    if val < val_min_subset
                        val_min_subset = val
                        ind_min_subset = ind_old
                        ind_min_option = ind_new
                    end

                    if val < val_min_total

                        inds_subset_temp .= inds_subset
                        sort!(inds_subset_temp)

                        flag = true
                        for ind in 1:keep
                            val < min_list[ind][1] && break
                            if val == min_list[ind][1] && inds_subset_temp == min_list[ind][2]
                                flag = false
                                break
                            end
                        end

                        if flag == true
                            num_out < keep && (num_out += 1)
                            min_list[1] = (val, min_list[1][2])
                            min_list[1][2] .= inds_subset_temp

                            sort!(min_list, by = x->x[1], rev = true)

                            val_min_total = min_list[1][1]
                        end
                    end
                end

                inds_subset[ind_old] = ind_swap
            end

            val_min_subset >= val_min && break

            val_min = val_min_subset

            (inds_subset[ind_min_subset], inds_option[ind_min_option]) = (inds_option[ind_min_option], inds_subset[ind_min_subset])
        end
    end

    out_mat = Matrix{eltype(inds_total)}(undef, subset_size, num_out)
    for ind in 1:num_out
        out_mat[:, ind] = min_list[keep-ind+1][2]
    end

    out_vec = Vector{typeof(val_min_total)}(undef, num_out)
    for ind in 1:num_out
        out_vec[ind] = min_list[keep-ind+1][1]
    end

    return out_vec, out_mat
end

end
