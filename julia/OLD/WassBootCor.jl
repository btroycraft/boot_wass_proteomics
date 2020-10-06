module WassBootCor

using Statistics, LIBSVM

export wbc_mat, cor_mat, wbc_cv, cor_cv, svm_close

wbc_mat = function(boot_reps, group_list, group_vec, data_mat, range, boot_mat; wp, trans)

  out = zeros(Float64, length(range), length(range))

  vec_x1_list = [Vector{Float64}(undef, length(group)) for group in group_list]
  vec_x2_list = deepcopy(vec_x1_list)

  vec_cor = Vector{Float64}(undef, boot_reps)
  vec_cor_mat = Matrix{Float64}(undef, length(group_list), boot_reps)

  vec_calc = Vector{Float64}(undef, length(group_list))

  for ind1 in 1:length(range)
    for ind2 in (ind1+1):length(range)

      val = wbc_pw(ind1, ind2, boot_reps, group_list, group_vec, data_mat, range, boot_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)
      out[ind1, ind2] = val
      out[ind2, ind1] = val
    end
  end

  return out
end

wbc_pw =  function(ind1, ind2, boot_reps, group_list, group_vec, data_mat, range, boot_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)

  for ind_group in 1:length(group_list)

    group = group_list[ind_group]
    vec_x1 = vec_x1_list[ind_group]
    vec_x2 = vec_x2_list[ind_group]

    ind_boot_col = 0
    for ind_boot in 1:boot_reps

      @label LOOP_START

      ind_boot_col += 1

      ind_boot_col == size(boot_mat, 2) && return NaN

      for ind_samp in 1:length(group)
        vec_x1[ind_samp] = data_mat[group_vec[boot_mat[group[ind_samp], ind_boot_col]], range[ind1]]
        vec_x2[ind_samp] = data_mat[group_vec[boot_mat[group[ind_samp], ind_boot_col]], range[ind2]]
      end

      val = trans == true ? atanh(cor(vec_x1, vec_x2)) : cor(vec_x1, vec_x2)

      isnan(val) || isinf(val) && @goto LOOP_START

      vec_cor[ind_boot] = val
    end

    sort!(vec_cor)
    vec_cor_mat[ind_group, 1:boot_reps] = vec_cor
  end

  out = 0.
  for ind_boot in 1:boot_reps
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

  out /= boot_reps*length(group_list)

  return out
end

cor_mat = function(group_list, group_vec, data_mat, range; lp, trans)

  out = zeros(Float64, length(range), length(range))

  vec_x1_list = [Vector{Float64}(undef, length(group)) for group in group_list]
  vec_x2_list = deepcopy(vec_x1_list)

  vec_calc = Vector{Float64}(undef, length(group_list))

  for ind1 in 1:length(range)
    for ind2 in (ind1+1):length(range)
      val = cor_pw(ind1, ind2, group_list, group_vec, data_mat, range, vec_x1_list, vec_x2_list, vec_calc; lp, trans)
      out[ind1, ind2] = val
      out[ind2, ind1] = val
    end
  end

  return out
end

cor_pw =  function(ind1, ind2, group_list, group_vec, data_mat, range, vec_x1_list, vec_x2_list, vec_calc; lp, trans)

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

  vec_x1_list = [Vector{Float64}(undef, length(group)) for group in group_list]
  vec_x2_list = deepcopy(vec_x1_list)

  vec_cor = Vector{Float64}(undef, size(boot_mat, 2))
  vec_cor_mat = Matrix{Float64}(undef, size(boot_mat, 2), length(group_list))

  vec_calc = Vector{Float64}(undef, length(group_list))

  return function(sub)
    return wbc_perm_(sub, group_list, group_vec, data_mat, range, boot_mat, perm_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc, group_vec_perm; wp, trans)
  end
end

wbc_perm_ =  function(sub, group_list, group_vec, data_mat, range, boot_mat, perm_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc, group_vec_perm; wp, trans)

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

wbc_sum_subset =  function(sub, group_list, group_vec, data_mat, range, boot_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)

  out = 0.
  for ind1 in 1:length(sub)
    for ind2 in (ind1+1):length(sub)
      val = wbc_pw(sub[ind1], sub[ind2], group_list, group_vec, data_mat, range, boot_mat, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)
      out += val
    end
  end

  return out / binomial(length(sub), 2)
end

cor_sum_subset =  function(sub, group_list, group_vec, data_mat, range, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; lp, trans)

  out = 0.
  for ind1 in 1:length(sub)
    for ind2 in (ind1+1):length(sub)
      val = cor_pw(sub[ind1], sub[ind2], group_list, group_vec, data_mat, range, vec_x1_list, vec_x2_list, vec_cor, vec_cor_mat, vec_calc; lp, trans)
      out += val
    end
  end

  return out / binomial(length(sub), 2)
end

wbc_cv = function(sub, boot_reps, train_list, test, test_lab, train_mat, group_vec, data_mat, range, boot_mat; wp, trans)

  train_vec = Vector{Int}(undef, size(train_mat, 1))
  vec_succ = Vector{Bool}(undef, size(train_mat, 2))

  vec_x1_test = Vector{Float64}(undef, length(test))
  vec_x2_test = copy(vec_x1_test)

  vec_x1_train_list = [Vector{Float64}(undef, length(train)) for train in train_list]
  vec_x2_train_list = deepcopy(vec_x1_train_list)

  cor_vec_test = Vector{Float64}(undef, boot_reps)
  cor_vec_train = copy(cor_vec_test)

  vec_sep = Vector{Float64}(undef, length(train_list))
  vec_vote = Vector{Int}(undef, binomial(length(sub), 2))

  for ind_train in 1:size(train_mat, 2)

    for ind_samp in 1:size(train_mat, 1)
      train_vec[ind_samp] = train_mat[ind_samp, ind_train]
    end

    val = wbc_cv_(sub, boot_reps, train_list, test, test_lab, train_vec, group_vec, data_mat, range, boot_mat, vec_x1_test, vec_x2_test, vec_x1_train_list, vec_x2_train_list, cor_vec_test, cor_vec_train, vec_sep, vec_vote; wp, trans)
    vec_succ[ind_train] = test_lab[ind_train] == val

  end

  return mean(vec_succ)
end

wbc_cv_ = function(sub, boot_reps, train_list, test, test_lab, train_vec, group_vec, data_mat, range, boot_mat, vec_x1_test, vec_x2_test, vec_x1_train_list, vec_x2_train_list, cor_vec_test, cor_vec_train, vec_sep, vec_vote; wp, trans)

  ind_vote = 1
  for ind_1 in 1:length(sub)
    for ind_2 in (ind_1+1):length(sub)

      ind_boot_col = 0
      for ind_boot in 1:boot_reps

        @label LOOP_START_TEST

        ind_boot_col += 1
        ind_boot_col == size(boot_mat, 2) && return 0

        for ind_samp in 1:length(test)
          vec_x1_test[ind_samp] = data_mat[group_vec[train_vec[boot_mat[test[ind_samp], ind_boot_col]]], range[sub[ind_1]]]
          vec_x2_test[ind_samp] = data_mat[group_vec[train_vec[boot_mat[test[ind_samp], ind_boot_col]]], range[sub[ind_2]]]
        end

        cor_test = trans == true ? atanh(cor(vec_x1_test, vec_x2_test)) : cor(vec_x1_test, vec_x2_test)

        (isnan(cor_test) || isinf(cor_test)) && @goto LOOP_START_TEST

        cor_vec_test[ind_boot] = cor_test
      end

      sort!(cor_vec_test)

      for ind_group in 1:length(train_list)
        train = train_list[ind_group]
        vec_x1_train = vec_x1_train_list[ind_group]
        vec_x2_train = vec_x2_train_list[ind_group]

        ind_boot_col = 0
        for ind_boot in 1:boot_reps

          @label LOOP_START_TRAIN

          ind_boot_col += 1
          ind_boot_col == size(boot_mat, 2) && return 0

          for ind_samp in 1:length(train)
            vec_x1_train[ind_samp] = data_mat[group_vec[train_vec[boot_mat[train[ind_samp], ind_boot_col]]], range[sub[ind_1]]]
            vec_x2_train[ind_samp] = data_mat[group_vec[train_vec[boot_mat[train[ind_samp], ind_boot_col]]], range[sub[ind_2]]]
          end

          cor_train = trans == true ? atanh(cor(vec_x1_train, vec_x2_train)) : cor(vec_x1_train, vec_x2_train)

          (isnan(cor_train) || isinf(cor_train)) && @goto LOOP_START_TRAIN

          cor_vec_train[ind_boot] = cor_train
        end

        sort!(cor_vec_train)

        total = 0.
        for ind_boot in 1:boot_reps
          if wp == 1
            total += abs(cor_vec_train[ind_boot] - cor_vec_test[ind_boot])
          elseif wp == 2
            total += (cor_vec_train[ind_boot] - cor_vec_test[ind_boot])^2
          end
        end

        vec_sep[ind_group] = total
      end

      vec_vote[ind_vote] = findmin(vec_sep)[2]
      ind_vote += 1
    end
  end
  return findmax([sum([vote == x for vote in vec_vote]) for x in 1:length(train_list)])[2]
end

cor_cv = function(sub, train_list, test, test_lab, train_mat, group_vec, data_mat, range; lp, trans)

  train_vec = Vector{Int}(undef, size(train_mat, 1))
  vec_succ = Vector{Bool}(undef, size(train_mat, 2))

  vec_x1_test = Vector{Float64}(undef, length(test))
  vec_x2_test = copy(vec_x1_test)

  vec_x1_train_list = [Vector{Float64}(undef, length(train)) for train in train_list]
  vec_x2_train_list = deepcopy(vec_x1_train_list)

  vec_sep = Vector{Float64}(undef, length(train_list))
  vec_vote = Vector{Int}(undef, binomial(length(sub), 2))

  for ind_train in 1:size(train_mat, 2)

    for ind_samp in 1:size(train_mat, 1)
      train_vec[ind_samp] = train_mat[ind_samp, ind_train]
    end

    val = cor_cv_(sub, train_list, test, test_lab, train_vec, group_vec, data_mat, range, vec_x1_test, vec_x2_test, vec_x1_train_list, vec_x2_train_list, vec_sep, vec_vote; lp, trans)
    vec_succ[ind_train] = test_lab[ind_train] == val
  end

  return mean(vec_succ)
end

cor_cv_ = function(sub, train_list, test, test_lab, train_vec, group_vec, data_mat, range, vec_x1_test, vec_x2_test, vec_x1_train_list, vec_x2_train_list, vec_sep, vec_vote; lp, trans)

  ind_vote = 1
  for ind_1 in 1:length(sub)
    for ind_2 in (ind_1+1):length(sub)

      for ind_samp in 1:length(test)
        vec_x1_test[ind_samp] = data_mat[group_vec[train_vec[test[ind_samp]]], range[sub[ind_1]]]
        vec_x2_test[ind_samp] = data_mat[group_vec[train_vec[test[ind_samp]]], range[sub[ind_2]]]
      end

      cor_test = trans == true ? atanh(cor(vec_x1_test, vec_x2_test)) : cor(vec_x1_test, vec_x2_test)

      for ind_group in 1:length(train_list)
        train = train_list[ind_group]
        vec_x1_train = vec_x1_train_list[ind_group]
        vec_x2_train = vec_x2_train_list[ind_group]

        for ind_samp in 1:length(train)
          vec_x1_train[ind_samp] = data_mat[group_vec[train_vec[train[ind_samp]]], range[sub[ind_1]]]
          vec_x2_train[ind_samp] = data_mat[group_vec[train_vec[train[ind_samp]]], range[sub[ind_2]]]
        end

        cor_train = trans == true ? atanh(cor(vec_x1_train, vec_x2_train)) : cor(vec_x1_train, vec_x2_train)

        vec_sep[ind_group] = abs(cor_train - cor_test)
      end

      vec_vote[ind_vote] = findmin(vec_sep)[2]
      ind_vote += 1
    end
  end

  return findmax([sum([vote == x for vote in vec_vote]) for x in 1:length(train_list)])[2]
end

svm_close = function(group_list, group_vec, data_mat, range, boot_mat, sub_size)

  vec_x1_list = [Vector{Float64}(undef, length(group)) for group in group_list]
  vec_x2_list = deepcopy(vec_x1_list)

  vec_cor_mat = Matrix{Float64}(undef, binomial(sub_size, 2), length(group_list)*size(boot_mat, 2))

  labels = repeat(1:length(group_list); inner = size(boot_mat, 2))

  return function(sub)

    for ind_group in 1:length(group_list)

      group = group_list[ind_group]
      vec_x1 = vec_x1_list[ind_group]
      vec_x2 = vec_x2_list[ind_group]

      ind_count = 1
      for ind1 in 1:length(sub)
        for ind2 in (ind1+1):length(sub)

          for ind_boot in 1:size(boot_mat, 2)
            for ind_samp in 1:length(group)
              vec_x1[ind_samp] = data_mat[group_vec[boot_mat[group[ind_samp], ind_boot]], range[ind1]]
              vec_x2[ind_samp] = data_mat[group_vec[boot_mat[group[ind_samp], ind_boot]], range[ind2]]
            end

            vec_cor_mat[ind_count, (ind_group-1)*size(boot_mat, 2) + ind_boot] = cor(vec_x1, vec_x2)
          end

          ind_count += 1
        end
      end
    end

    model = svmtrain(vec_cor_mat, labels)
    (p_lab, dec) = svmpredict(model, vec_cor_mat)

    return mean(p_lab .== labels)
  end
end

end
