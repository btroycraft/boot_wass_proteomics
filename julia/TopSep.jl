module TopSep

using Ripser, Statistics

export top_sep_alloc, TopSepCall, top_sep_group_alloc, TopSepGroupCall

function top_sep_alloc(sub_size, group_list, group_vec, data_mat, range = 1:size(data_mat, 2))

    dist_mat_list = [1 .- cor(data_mat[group_vec[group], range]) for group in group_list]
    dist_sub_mat_list = [zeros(sub_size, sub_size) for group in group_list]

    return dist_mat_list, dist_sub_mat_list
end

struct TopSepCall

    dim::Int
    dist_mat_list::Vector{Matrix{Float64}}
    dist_sub_mat_list::Vector{Matrix{Float64}}
end

function (TopSepCall_::TopSepCall)(sub)

    sub_size = length(sub)
    num_groups = length(TopSepCall_.dist_mat_list)

    for ind_group in 1:num_groups

        dist_mat = TopSepCall_.dist_mat_list[ind_group]
        dist_sub_mat = TopSepCall_.dist_sub_mat_list[ind_group]

        for ind2 in 1:sub_size
            for ind1 in 1:sub_size
                dist_sub_mat[ind1, ind2] = dist_mat[sub[ind1], sub[ind2]]
            end
        end
    end

    diag_list = map(TopSepCall_.dist_sub_mat_list) do dist_sub_mat

        diag = ripser(dist_sub_mat;
            dim_max = TopSepCall_.dim)[TopSepCall_.dim+1]
        birth = map(diag) do x
            return x[1]
        end
        death = map(diag) do x
            if x[2] == Inf
                return prevfloat(Inf)
            else
                return x[2]
            end
        end

        return birth, death
    end

    out = Inf

    for ind1 in 1:(num_groups-1)

        (birth1, death1) = diag_list[ind1]
        for ind2 in (ind1+1):num_groups

            (birth2, death2) = diag_list[ind2]

            breaks = sort!(unique(vcat(birth1, death1, birth2, death2)))

            betti1 = map(breaks) do x
                sum((birth1 .<= x) .* (death1 .> x))
            end

            betti2 = map(breaks) do x
                sum((birth2 .<= x) .* (death2 .> x))
            end

            dist = 0.
            for ind in 1:(length(breaks)-1)
                dist += abs(betti1[ind] - betti2[ind]) * (breaks[ind+1] - breaks[ind])
            end

            out = min(out, dist)
        end
    end

    return out
end

function top_sep_group_alloc(sub_size, group_num, group_list, group_vec, data_mat, range = 1:size(data_mat, 2))

    dist_mat_list = [1 .- cor(data_mat[group_vec[group], range]) for group in group_list]
    dist_sub_mat_list = [zeros(sub_size, sub_size) for group in group_list]

    temp = dist_mat_list[group_num]
    for ind in 1:(group_num-1)
        dist_mat_list[ind+1] = dist_mat_list[ind]
    end
    dist_mat_list[1] = temp

    return dist_mat_list, dist_sub_mat_list
end

struct TopSepGroupCall

    dim::Int
    dist_mat_list::Vector{Matrix{Float64}}
    dist_sub_mat_list::Vector{Matrix{Float64}}
end

function (TopSepGroupCall_::TopSepGroupCall)(sub)

    sub_size = length(sub)
    num_groups = length(TopSepGroupCall_.dist_mat_list)

    for ind_group in 1:num_groups

        dist_mat = TopSepGroupCall_.dist_mat_list[ind_group]
        dist_sub_mat = TopSepGroupCall_.dist_sub_mat_list[ind_group]

        for ind2 in 1:sub_size
            for ind1 in 1:sub_size
                dist_sub_mat[ind1, ind2] = dist_mat[sub[ind1], sub[ind2]]
            end
        end
    end

    diag_list = map(TopSepGroupCall_.dist_sub_mat_list) do dist_sub_mat

        diag = ripser(dist_sub_mat;
            dim_max = TopSepGroupCall_.dim)[TopSepGroupCall_.dim+1]
        birth = map(diag) do x
            return x[1]
        end
        death = map(diag) do x
            if x[2] == Inf
                return prevfloat(Inf)
            else
                return x[2]
            end
        end

        return birth, death
    end

    out = Inf

    (birth1, death1) = diag_list[1]
    for ind2 in 2:num_groups

        (birth2, death2) = diag_list[ind2]

        breaks = sort!(unique(vcat(birth1, death1, birth2, death2)))

        betti1 = map(breaks) do x
            sum((birth1 .<= x) .* (death1 .> x))
        end

        betti2 = map(breaks) do x
            sum((birth2 .<= x) .* (death2 .> x))
        end

        dist = 0.
        for ind in 1:(length(breaks)-1)
            dist += abs(betti1[ind] - betti2[ind]) * (breaks[ind+1] - breaks[ind])
        end

        out = min(out, dist)
    end

    return out
end

end
