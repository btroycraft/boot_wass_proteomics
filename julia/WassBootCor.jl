using StatsBase: sample!
using Random: shuffle!
using Statistics: mean, cor, median

wbc_mat = function(group_list, group_vec, data_mat, range, boot_mat; wp, trans)

  out = zeros(Float64, length(range), length(range))

  vec_x1_list = [Vector{Float64}(undef, length(group)) for group in group_list]
  vec_x2_list = deepcopy(vec_x1_list)

  vec_cor = Vector{Float64}(undef, size(boot_mat, 2))
  vec_cor_mat = Matrix{Float64}(undef, length(group_list), size(boot_mat, 2))

  vec_calc = Vector{Float64}(undef, length(group_list))

  for ind1 in 1:length(range)
    for ind2 in (ind1+1):length(range)

      val = wbc_pw(ind1, ind2, group_list, group_vec, data_mat, range, boot_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)
      out[ind1, ind2] = val
      out[ind2, ind1] = val
    end
  end

  return out
end

cor_mat = function(group_list, group_vec, data_mat, range; lp, trans)

  out = zeros(Float64, length(range), length(range))

  vec_x1 = [Vector{Float64}(undef, length(group)) for group in group_list]
  vec_x2 = deepcopy(vec_x1)

  for ind1 in 1:length(range)
    for ind2 in (ind1+1):length(range)
      val = cor_pw(ind1, ind2, group_list, group_vec, data_mat, range, vec_x1_list, vec_x2_list, vec_calc; lp, trans)
      out[ind1, ind2] = val
      out[ind2, ind1] = val
    end
  end

  return out
end

wbc_pw = @inline function(ind1, ind2, group_list, group_vec, data_mat, range, boot_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)

  for ind_group in 1:length(group_list)

    group = group_list[ind_group]
    vec_x1 = vec_x1_list[ind_group]
    vec_x2 = vec_x2_list[ind_group]

    for ind_boot in 1:size(boot_mat, 2)
      for ind_samp in 1:length(group)
        vec_x1[ind_samp] = data_mat[group_vec[boot_mat[group[ind_samp], ind_boot]], range[ind1]]
        vec_x2[ind_samp] = data_mat[group_vec[boot_mat[group[ind_samp], ind_boot]], range[ind2]]
      end

      vec_cor[ind_boot] = trans == true ? atanh(cor(vec_x1, vec_x2)) : cor(vec_x1, vec_x2)
    end

    sort!(vec_cor)
    vec_cor_mat[ind_group, 1:size(boot_mat, 2)] = vec_cor
  end

  out = 0.
  for ind_boot in 1:size(boot_mat, 2)
    for ind_group in 1:length(group_list)
      vec_calc[ind_group] = vec_cor_mat[ind_group, ind_boot]
    end

    if wp == 1
      vec_calc .-= median(vec_calc)
      vec_calc .= abs.(vec_calc)
    elseif wp == 2
      vec_calc .-= mean(vec_calc)
      vec_calc .= vec_calc.^2
    end
    out += sum(vec_calc)
  end

  out /= size(boot_mat, 2)*length(group_list)

  return out
end

cor_pw = @inline function(ind1, ind2, group_list, group_vec, data_mat, range, vec_x1, vec_x2, vec_calc; lp, trans)

  for ind_group in 1:length(group_list)

    group = group_list[ind_group]
    vec_x1 = vec_x1_list[ind_group]
    vec_x2 = vec_x2_list[ind_group]

    for ind_samp in 1:length(group)
      vec_x1[ind_samp] = data_mat[group_vec[group[ind_samp]], range[ind1]]
      vec_x2[ind_samp] = data_mat[group_vec[group[ind_samp]], range[ind2]]
    end

    vec_calc[ind_group] = trans == true ? atanh(cor(vec_x1, vec_x2)) : cor(vec_x1, vec_x2)
  end

  if lp == 1
    vec_calc .-= median(vec_calc)
    vec_calc .= abs.(vec_calc)
  elseif lp == 2
    vec_calc .-= mean(vec_calc)
    vec_calc .= vec_calc.^2
  end
  out = sum(vec_calc) / length(group_list)

  return out
end

wbc_perm_close = function(group_list, group_vec, data_mat, range, boot_mat, perm_mat; wp, trans)

  group_vec_perm = copy(group_vec)

  vec_x1_list = [Vector{Float64}(undef, size(group, 1)) for group in 1:length(group_list)]
  vec_x2_list = deepcopy(vec_x1_list)

  vec_cor = Vector{Float64}(undef, size(boot_mat, 2))
  vec_cor_mat = Matrix{Float64}(undef, size(boot_mat, 2), length(group_list))

  vec_calc = Vector{Float64}(undef, length(group_list))

  return function(sub)
    return wbc_perm(sub, group_list, group_vec, data_mat, range, boot_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc, group_vec_perm; wp, trans)
  end
end

wbc_perm = @inline function(sub, group_list, group_vec, data_mat, range, boot_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc, group_vec_perm; wp, trans)

  sum_orig = wbc_sum_subset(sub, group_list, group_vec, data_mat, range, boot_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)

  below = 0
  above = 1

  for ind_perm in 1:size(perm_mat, 2)
    for ind_samp in 1:length(group_vec)
      group_vec_perm[ind_samp] = group_vec[perm_mat[ind_samp, ind_perm]]
    end

    val = wbc_sum_subset(sub, group_list, group_vec_perm, data_mat, range, boot_list, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)

    below += val < sum_orig
    above += val >= sum_orig
  end

  return convert(Float64, above) / (convert(Float64, above) + convert(Float64, below))
end

wbc_sum_subset = @inline function(sub, group_list, group_vec, data_mat, range, boot_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)

  out = 0.
  for ind1 in 1:length(sub)
    for ind2 in (ind1+1):length(sub)
      val = wbc_pw(sub[ind1], sub[ind2], group_list, group_vec, data_mat, range, boot_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)
      out += val
    end
  end

  return out / binomial(length(sub), 2)
end
