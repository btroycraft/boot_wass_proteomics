module WassBootCor

using StatsBase: sample!
using Statistics: cor, var, mean, median

export wass_boot_cor_sep, cor_sep

wass_boot_cor_sep = function(X; groups, boot_reps = 100, index_range = 1:size(X, 2), wp = "w2", transformed = true)

  inds_total = vcat(index_range)
  range_size = length(inds_total)

  # Find group labels from input

  inds_group_list = [findall(x -> x == ind, groups) for ind in unique(groups)]

  n_groups = length(inds_group_list)

  # preallocate vectors

  inds_boot_list = Vector{Matrix{Int}}(undef, n_groups)
  x_temp_list = Vector{Vector{eltype(X)}}(undef, n_groups)
  y_temp_list = similar(x_temp_list)

  # setup bootstrap samples, need var > 0

  for (group_size, ind) in ((length(inds_group_list[ind]), ind) for ind in 1:n_groups)

    inds_cand = Vector{Int}(undef, group_size)
    vec_cand = Vector{Float64}(undef, group_size)

    inds_boot_list[ind] = Matrix{Int}(undef, group_size, boot_reps)

    for ind1 in 1:boot_reps
      while true
        flag = true
        sample!(inds_group_list[ind], inds_cand; replace = true)
        for ind2 in inds_total
          for ind3 in 1:group_size
            vec_cand[ind3] = X[inds_cand[ind3], ind2]
          end

          if var(vec_cand) == 0
            flag = false
            break
          end
        end

        if flag
          inds_boot_list[ind][:, ind1] = inds_cand
          break
        end
      end
    end

    x_temp_list[ind] = Vector{Float64}(undef, group_size)
    y_temp_list[ind] = Vector{Float64}(undef, group_size)
  end

  cor_vec_list = Vector{Vector{Float64}}(undef, n_groups)
  for ind in 1:n_groups
    cor_vec_list[ind] = Vector{Float64}(undef, boot_reps)
  end

  temp_vec = Vector{Float64}(undef, n_groups)
  var_vec = Vector{Float64}(undef, boot_reps)

  wass_sep_out = zeros(Float64, range_size, range_size)

  # start wass calculations

  for ind_x in 1:range_size
    for ind_y in (ind_x+1):range_size

      for ind_group in 1:n_groups

        group_size = length(inds_group_list[ind_group])

        inds_boot = inds_boot_list[ind_group]

        cor_vec = cor_vec_list[ind_group]

        x_temp = x_temp_list[ind_group]
        y_temp = y_temp_list[ind_group]

        # for (i, j) pair and group number, get correlation from bootstrap indices

        for ind_boot in 1:boot_reps
          for ind in 1:group_size
            x_temp[ind] = X[inds_boot[ind, ind_boot], inds_total[ind_x]]
          end
          for ind in 1:group_size
            y_temp[ind] = X[inds_boot[ind, ind_boot], inds_total[ind_y]]
          end

          cor_vec[ind_boot] = transformed ? atanh(cor(x_temp, y_temp)) : cor(x_temp, y_temp)
        end

        sort!(cor_vec)
      end

      wass_sep = 0.
      for ind1 in 1:boot_reps

        for ind2 in 1:n_groups
          temp_vec[ind2] = cor_vec_list[ind2][ind1]
        end

        # get contribution to Wp distance

        if wp == "w1"
          med_cor = median(temp_vec)
          @. temp_vec = temp_vec - med_cor
          @. temp_vec = abs(temp_vec)
          sum_cor = sum(temp_vec)

          wass_sep += sum_cor
        elseif wp == "w2"
          mean_cor = mean(temp_vec)
          @. temp_vec = temp_vec - mean_cor
          @. temp_vec = temp_vec^2
          sum_cor = sum(temp_vec)

          wass_sep += sum_cor
        end

      end

      wass_sep_out[ind_x, ind_y] = wass_sep/boot_reps
      wass_sep_out[ind_y, ind_x] = wass_sep/boot_reps
    end
  end

  return wass_sep_out
end

cor_sep = function(X; groups, boot_reps = 100, index_range = 1:size(X, 2), lp = "l2", transformed = true)

  inds_total = vcat(index_range)
  range_size = length(inds_total)

  inds_group_list = [findall(x -> x == ind, groups) for ind in unique(groups)]

  n_groups = length(inds_group_list)

  # preallocate vectors

  x_temp_list = Vector{Vector{eltype(X)}}(undef, n_groups)
  y_temp_list = similar(x_temp_list)

  for (group_size, ind) in ((length(inds_group_list[ind]), ind) for ind in 1:n_groups)
    x_temp_list[ind] = Vector{Float64}(undef, group_size)
    y_temp_list[ind] = Vector{Float64}(undef, group_size)
  end

  temp_vec = Vector{Float64}(undef, n_groups)

  cor_sep_out = zeros(Float64, range_size, range_size)

  for ind_x in 1:range_size
    for ind_y in (ind_x+1):range_size

      for ind_group in 1:n_groups

        inds_group = inds_group_list[ind_group]
        group_size = length(inds_group)

        x_temp = x_temp_list[ind_group]
        y_temp = y_temp_list[ind_group]

        # for (i, j) pair and group number, get correlation

        for ind in 1:group_size
          x_temp[ind] = X[inds_group[ind], inds_total[ind_x]]
        end
        for ind in 1:group_size
          y_temp[ind] = X[inds_group[ind], inds_total[ind_y]]
        end

        temp_vec[ind_group] = transformed ? atanh(cor(x_temp, y_temp)) : cor(x_temp, y_temp)
      end

      # get correlation separation

      if lp == "l1"
        med_cor = median(temp_vec)
        @. temp_vec = temp_vec - med_cor
        @. temp_vec = abs(temp_vec)
      elseif lp == "l2"
        mean_cor = mean(temp_vec)
        @. temp_vec = temp_vec - mean_cor
        @. temp_vec = temp_vec^2
      end

      cor_sep_out[ind_x, ind_y] = cor_sep_out[ind_y, ind_x] = sum(temp_vec)
    end
  end

  return cor_sep_out
end

end
