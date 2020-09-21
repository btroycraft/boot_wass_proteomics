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
      if isnan(vec_cor[ind_boot])
        @show(vec_x1, vec_x2)
      end
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
    return wbc_perm(sub, group_list, group_vec, data_mat, range, boot_mat, perm_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc, group_vec_perm; wp, trans)
  end
end

wbc_perm = @inline function(sub, group_list, group_vec, data_mat, range, boot_mat, perm_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc, group_vec_perm; wp, trans)

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

cor_sum_subset = @inline function(sub, group_list, group_vec, data_mat, range, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; lp, trans)

  out = 0.
  for ind1 in 1:length(sub)
    for ind2 in (ind1+1):length(sub)
      val = cor_pw(sub[ind1], sub[ind2], group_list, group_vec, data_mat, range, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; lp, trans)
      out += val
    end
  end

  return out / binomial(length(sub), 2)
end

wbc_cv = function(sub_mat, train_list, test, group_vec_mat, data_mat, range, boot_mat; wp, trans)

  group_vec_perm = Vector{Int}(undef, size(group_vec_mat, 1))
  succ = Vector{Int}(undef, size(sub_mat, 2))

  vec_x1_test = Vector{Float64}(undef, length(test))
  vec_x2_test = copy(vec_x1_test)

  vec_x1_train = [Vector{Float64}(undef, length(train)) for train in train_list]
  vec_x2_train = deepcopy(vec_x1_train)

  cor_vec_test = Vector{Float64}(undef, size(boot_mat, 2))
  cor_vec_train = copy(cor_vec_test)

  wass = Vector{Float64}(undef, length(train_list))
  vote = Vector{Int}(undef, binomial(size(sub_mat, 1), 2))

  succ_list = map(1:size(group_vec_mat, 2)) do ind_rep

    succ = zeros(Int, size(sub_mat, 2))

    for ind in 1:size(group_vec_mat, 1)
      group_vec_perm[ind] = group_vec_mat[ind, ind_rep]
    end

    for ind_sub in 1:size(sub_mat, 2)
        for ind in 1:size(sub_mat, 1)
          sub[ind] = sub_mat[ind, ind_sub]
        end

        succ[ind_sub] = lab_vec[ind_rep] == wbc_cv_(sub, train_list, test, group_vec_perm, data_mat, range, boot_mat, vec_x1_test, vec_x2_test, vec_x1_list_train, vec_x2_list_train, cor_vec_test, cor_vec_train, wass, vote; wp, trans)
    end

    return succ
  end

  return map(1:size(sub_mat, 2)) do ind_sub
    mean([succ[ind_sub] for succ in succ_list])
  end
end

wbc_cv_ = function(sub, train_list, test, group_vec, data_mat, range, boot_mat, vec_x1_test, vec_x2_test, vec_x1_list_train, vec_x2_list_train, cor_vec_test, cor_vec_train, wass, vote; wp, trans)

  ind_vote = 1
  for ind1 in 1:length(sub)
    for ind2 in (ind1+1):length(sub)

      for ind_boot in 1:size(boot_mat, 2)
        for ind_samp in 1:length(test)
          vec_x1_test[ind_samp] = data_mat[group_vec[boot_mat[test[ind_samp], ind_boot]], range[ind_1]]
          vec_x2_test[ind_samp] = data_mat[group_vec[boot_mat[test[ind_samp], ind_boot]], range[ind_2]]
        end

        cor_vec_test[ind_boot] = trans == true ? atanh(cor(vec_x1_test, vec_x2_test)) : cor(vec_x1_test, vec_x2_test)
      end

      sort!(core_vec_test)

      for ind_group in 1:length(train_list)
        train = train_list[ind_group]
        vec_x1_train = vec_x1_train_list[ind_group]
        vec_x2_train = vec_x2_train_list[ind_group]

        for ind_samp in 1:length(train)
          vec_x1_train[ind_samp] = data_mat[group_vec[boot_mat[train[ind_samp], ind_boot]], range[ind_1]]
          vec_x2_train[ind_samp] = data_mat[group_vec[boot_mat[train[ind_samp], ind_boot]], range[ind_2]]
        end

        cor_vec_train[ind_boot] = trans == true ? atanh(cor(vec_x1_train, vec_x2_train)) : cor(vec_x1_train, vec_x2_train)

        sort!(cor_vec_train)

        total = 0.
        for ind_boot in 1:size(boot_mat, 2)
          if wp == 1
            total += abs(cor_vec_train[ind_boot] - cor_vec_test[ind_boot])
          elseif wp == 2
            total += (cor_vec_train[ind_boot] - cor_vec_test[ind_boot])^2
          end
        end

        wass[ind_group] = total
      end

      vote[ind_vote] = findmin(wass)
      ind_vote += 1
    end
  end

  return findmax(map(x -> sum(vote == x), 1:length(list_train)))
end
