using StatsBase: sample!
using Random: shuffle!

split_groups = function(group_lab_vec)
  group_enum = collect(enumerate(group_lab_vec))
  sort!(group_enum; by = x -> x[2])
  out_vec = map(x -> x[1], group_enum)

  out_list = map(unique(group_lab_vec)) do lab
    lower = findfirst(x -> x[2] == lab, group_enum)
    upper = findlast(x -> x[2] == lab, group_enum)
    return lower:upper
  end

  return out_list, out_vec
end

boot_reps_cor = function(boot_reps, group_list, group_vec, data_mat, range)

  out = Matrix{Int}(undef, length(group_vec), boot_reps)

  vec_X = Vector{Float64}(undef, length(group_vec))
  vec_samp = Vector{Int}(undef, length(group_vec))

  for ind_boot in 1:boot_reps
    for ind_group in 1:length(group_list)

      group = group_list[ind_group]
      length_samp = length(group)

      @label LOOP_START
      sample!(group, vec_samp;
        replace = true)

      for ind_range in 1:length(range)
        for ind_samp in 1:length(group)
          vec_X[ind_samp] = data_mat[group_vec[vec_samp[group[ind_samp]]], range[ind_range]]
        end

        all(vec_X[ind] == vec_X[1] for ind in 2:length(group)) && @goto LOOP_START
      end

      for ind_samp in 1:length(group)
        out[group[ind_samp], ind_boot] = vec_samp[ind_samp]
      end
    end
  end

  return out
end

perm_cor = function(perm_reps, group_list, group_vec, data_mat, range)

  out = Matrix{Int}(undef, length(group_vec), perm_reps)

  vec_X = Vector{Float64}(undef, length(group_vec))
  vec_samp = collect(1:length(group_vec))

  for ind_perm in 1:perm_reps

    @label LOOP_START
    shuffle!(vec_samp)

    for ind_group in 1:length(group_list)

      group = group_list[ind_group]

      for ind_range in 1:length(range)
        for ind_samp in 1:length(group)
          vec_X[ind_samp] = data_mat[group_vec[vec_samp[group[ind_samp]]], range[ind_range]]
        end

        all(vec_X[ind] == vec_X[1] for ind in 2:length(group)) && @goto LOOP_START
      end
    end

    for ind_samp in 1:length(group_vec)
      out[ind_samp, ind_perm] = vec_samp[ind_samp]
    end
  end

  return out
end
