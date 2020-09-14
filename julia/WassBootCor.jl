module WassBootCor

export wbc_mat, cor_mat, split_groups, boot_reps_cor, wbc_pw_close, cor_pw_close,

wbc_mat = function(X, group_list, boot_list, range_arg = 1:size(X, 2); wp, trans)

  wbc_pw = wbc_pw_close(X, boot_list, range_arg; wp, trans)

  out = func_mat(wbc_pw, length(range_arg))

  return out
end

cor_mat = function(X, group_list, range_arg = 1:size(X, 2); wp, trans)

  cor_pw = cor_pw_close(X, group_list, range_arg; wp, trans)

  out = func_mat(cor_pw, length(range_arg))

  return out
end

split_groups = function(group_vec)

  group_enum = [x for x in enumerate(group_vec)]
  sort!(group_enum;
    by = x -> x[2])
  out_group = map(group_enum) do x
    x[1]
  end

  out = map(unique(group_vec)) do lab1

    lower = findfirst(group_enum) do x
      x[2] == lab1
    end

    upper = findlast(group_enum) do x
      x[2] == lab1
    end

    return view(out_group, lower:upper)
  end

  return out
end

func_mat = function(func, length_range)

  out_type = typeof(func(1, 1))

  out = Matrix{out_type}(undef, length_range, length_range)

  for ind1 in 1:length_range, ind2 in (ind2+1):length_range
      out[ind1, ind2] = out[ind2, ind1] = func(ind1, ind2)
  end

  return out
end

boot_reps_cor = function(length_boot, group_list, X, range_arg = 1:size(X, 2))

  length_range = length(range_arg)
  length_group = length(group_list)

  out_mat = Matrix{Int}(undef, sum(map(x -> length(x), group_list)), length_boot)

  lower = 0
  out_list = map(group_list) do x
    length_samp = length(x)
    out = view(out_mat, (lower+1):(lower+length_samp), 1:length_boot)
    lower += length_samp

    return out
  end

  for ind_group in 1:length_group

    group = group_list[ind_group]
    length_samp = length(group)
    boot = out_list[ind_group]
    vec_X = Vector{eltype(X)}(undef, length_samp)
    vec_samp = Vector{Int}(undef, length_samp)

    for ind_boot in 1:length_boot

      @label LOOP
      begin

        StatsBase.sample!(group, vec_samp; replace=true)

        for ind_range in 1:length_range
          for ind_samp in 1:length_samp
            vec_X[ind_samp] = X[vec_samp[ind_samp], range_arg[ind_range]]
          end

          if all(vec_X[ind_samp] == vec_X[1] for ind_samp in 1:length_samp)
            @goto LOOP
          end
        end
      end

      boot[1:length_samp, ind_boot] = vec_samp
    end
  end

  return out_list
end

wbc_pw_close = function(X, boot_list, range_arg = 1:size(X, 2); wp, trans)

  out_type = float(eltype(X))

  length_group = length(boot_list)
  length_range = length(range_arg)
  length_boot = size(boot_list[1], 2)

  vec_X1_list = map(boot_list) do boot
    Vector{out_type}(undef, size(boot, 1))
  end
  vec_X2_list = deepcopy(vec_X2_list)

  vec_cor = Vector{out_type}(undef, length_boot)
  vec_cor_mat = Matrix{out_type}(undef, length_boot, length_group)

  vec_calc = Vector{out_type}(undef, length_group)

  return function(ind1, ind2)

    for ind_group in 1:length_group

      boot = boot_list[ind_group]

      length_samp = size(boot, 1)

      vec_X1 = vec_X1_list[ind_group]
      vec_X2 = vec_X2_list[ind_group]

      for ind_boot in 1:length_boot

        for ind_samp in 1:length_samp
          vec_X1[ind_samp] = X[boot[ind_samp, ind_boot], range_arg[ind1]]
          vec_X2[ind_samp] = X[boot[ind_samp, ind_boot], range_arg[ind2]]
        end

        if trans == true
          vec_cor[ind_boot] = atanh(Statistics.cor(vec_X1, vec_X2))
        else
          vec_cor[ind_boot] = Statistics.cor(vec_X1, vec_X2)
        end
      end

      sort!(vec_cor)
      vec_cor_mat[1:length_boot, ind_group] = vec_cor
    end

    out = zero(out_type)
    for ind_boot in 1:length_boot

      for ind_group in 1:length_group
        vec_calc[ind_group] = vec_cor_mat[ind_boot, ind_group]
      end

      if wp == 1
        vec_calc .-= Statistics.median(vec_calc)
        vec_calc .= abs.(vec_calc)
      elseif wp == 2

        vec_calc .-= Statistics.mean(vec_calc)
        vec_calc .= vec_calc.^2
      end

      out += sum(vec_calc)
    end

    out /= length_boot*length_group

    return out
  end
end

cor_pw_close = function(X, group_list, range_arg = 1:size(X, 2); lp, trans)

  out_type = float(eltype(X))

  length_group = length(group_list)
  length_range = length(range_arg)

  vec_X1_list = map(group_list) do group
    Vector{out_type}(undef, size(group, 1))
  end
  vec_X2_list = deepcopy(vec_X2_list)

  vec_calc = Vector{out_type}(undef, length_group)

  return function(ind1, ind2)

    for ind_group in 1:length_group

      group = group_list[ind_group]
      length_samp = length(group)

      vec_X1 = vec_X1_list[ind_group]
      vec_X2 = vec_X2_list[ind_group]

      for ind_samp in 1:length_samp
        vec_X1[ind_samp] = X[ind_samp, range_arg[ind1]]
        vec_X2[ind_samp] = X[ind_samp, range_arg[ind2]]
      end

      if trans == true
        vec_calc[ind_group] = atanh(Statistics.cor(vec_X1, vec_X2))
      else
        vec_calc[ind_group] = Statistics.cor(vec_X1, vec_X2)
      end
    end

    if lp == 1
      vec_calc .-= Statistics.median(vec_calc)
      vec_calc .= abs.(vec_calc)
    elseif lp == 2
      vec_calc .-= Statistics.mean(vec_calc)
      vec_calc .= vec_calc.^2
    end

    out = sum(vec_calc) / length_group

    return out
  end
end

perm_cor = function(length_boot, group_list, X, range_arg = 1:size(X, 2))

  length_range = length(range_arg)
  length_group = length(group_list)

  out_mat = Matrix{Int}(undef, sum(map(x -> length(x), group_list)), length_boot)

  lower = 0
  out_list = map(group_list) do x
    length_samp = length(x)
    out = view(out_mat, (lower+1):(lower+length_samp), 1:length_boot)
    lower += length_samp

    return out
  end

  for ind_group in 1:length_group

    group = group_list[ind_group]
    length_samp = length(group)
    boot = out_list[ind_group]
    vec_X = Vector{eltype(X)}(undef, length_samp)
    vec_samp = Vector{Int}(undef, length_samp)

    for ind_boot in 1:length_boot

      @label LOOP
      begin

        StatsBase.sample!(group, vec_samp; replace=true)

        for ind_range in 1:length_range
          for ind_samp in 1:length_samp
            vec_X[ind_samp] = X[vec_samp[ind_samp], range_arg[ind_range]]
          end

          if all(vec_X[ind_samp] == vec_X[1] for ind_samp in 1:length_samp)
            @goto LOOP
          end
        end
      end

      boot[1:length_samp, ind_boot] = vec_samp
    end
  end
