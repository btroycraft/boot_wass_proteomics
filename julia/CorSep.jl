module CorSep

using Statistics

export cor_sep_alloc, cor_sep!, cor_sep_group_alloc, cor_sep_group!

function cor_sep_alloc(sub_size, group_list, group_vec, data_mat, range; trans = false)

    num_groups = length(group_list)
    num_pairs = binomial(sub_size, 2)

    cor_mat_list = [trans == true ? atanh.(cor(data_mat[group_vec[group], range])) : cor(data_mat[group_vec[group], range]) for group in group_list]

    vec_cor_list = [Vector{Float64}(undef, num_pairs) for _ in 1:num_groups]

    return cor_mat_list, vec_cor_list
end

function cor_sep!(sub, group_list, group_vec, data_mat, range, cor_mat_list, vec_cor_list; trans = false)

    num_groups = length(group_list)
    sub_size = length(sub)
    num_pairs = binomial(sub_size, 2)

    for ind_group in 1:num_groups

        cor_mat = cor_mat_list[ind_group]
        vec_cor = vec_cor_list[ind_group]

        ind_cor = 1
        for ind1 in 1:sub_size
            for ind2 in (ind1+1):sub_size
                vec_cor[ind_cor] = cor_mat[sub[ind1], sub[ind2]]
                ind_cor += 1
            end
        end
    end

    out = Inf

    for ind1 in 1:num_groups
        vec_cor1 = vec_cor_list[ind1]
        for ind2 in (ind1+1):num_groups
            vec_cor2 = vec_cor_list[ind2]
            dist = 0.
            for ind_cor in 1:num_pairs
                dist += (vec_cor1[ind_cor] - vec_cor2[ind_cor])^2
            end
            out = min(out, dist)
        end
    end

    return sqrt(out)::Float64
end

function cor_sep_group_alloc(sub_size, group_num, group_list, group_vec, data_mat, range; trans = false)

    num_groups = length(group_list)
    num_pairs = binomial(sub_size, 2)

    not_group = vcat(1:(group_num-1), (group_num+1):(num_groups))

    cor_mat_list = [trans == true ? atanh.(cor(data_mat[group_vec[group], range])) : cor(data_mat[group_vec[group], range]) for group in group_list]

    vec_cor_list = [Vector{Float64}(undef, num_pairs) for _ in 1:num_groups]

    return not_group, cor_mat_list, vec_cor_list
end

function cor_sep_group!(sub, group_list, group_vec, data_mat, range, not_group, cor_mat_list, vec_cor_list; trans = false)

    sub_size = length(sub)
    num_groups = length(group_list)
    num_pairs = binomial(sub_size, 2)

    for ind_group in 1:num_groups

        cor_mat = cor_mat_list[ind_group]
        vec_cor = vec_cor_list[ind_group]

        ind_cor = 1
        for ind1 in 1:sub_size
            for ind2 in (ind1+1):sub_size
                vec_cor[ind_cor] = cor_mat[sub[ind1], sub[ind2]]
                ind_cor += 1
            end
        end
    end

    out = Inf

    vec_cor1 = vec_cor_list[group_num]
    for ind2 in not_group
        vec_cor2 = vec_cor_list[ind2]
        dist = 0.
        for ind_cor in 1:num_pairs
            dist += (vec_cor1[ind_cor] - vec_cor2[ind_cor])^2
        end
        out = min(out, dist)
    end

    return sqrt(out)::Float64
end

end
