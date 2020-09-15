using StatsBase: sample!
using Random: shuffle!
using Statistics: mean, cor, median

wbc_mat = function(X, group_list, boot_list, range_arg = 1:size(X, 2); wp, trans)

  length_group = length(group_list)
  length_range = length(range_arg)
  length_boot = size(boot_list[1], 2)

  vec_X1_list = map(boot_list) do boot
    Vector{Float64}(undef, size(boot, 1))
  end
  vec_X2_list = deepcopy(vec_X1_list)

  vec_cor = Vector{Float64}(undef, length_boot)
  vec_cor_mat = Matrix{Float64}(undef, length_boot, length_group)

  vec_calc = Vector{Float64}(undef, length_group)

  out = zeros(Float64, length_range, length_range)

  for ind1 in 1:length_range
    for ind2 in (ind1+1):length_range
      val = wbc_pw(ind1, ind2, X, group_list, boot_list, range_arg, vec_X1_list, vec_X2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)
      @inbounds out[ind1, ind2] = val
      @inbounds out[ind2, ind1] = val
    end
  end

  return out
end

cor_mat = function(X, group_list, range_arg = 1:size(X, 2); lp, trans)

  length_group = length(group_list)
  length_range = length(range_arg)

  vec_X1_list = map(group_list) do group
    Vector{Float64}(undef, length(group))
  end
  vec_X2_list = deepcopy(vec_X1_list)

  vec_calc = Vector{Float64}(undef, length_group)

  out = zeros(Float64, length_range, length_range)

  for ind1 in 1:length_range
    for ind2 in (ind1+1):length_range
      val = cor_pw(ind1, ind2, X, group_list, range_arg, vec_X1_list, vec_X2_list, vec_calc; lp, trans)
      @inbounds out[ind1, ind2] = val
      @inbounds out[ind2, ind1] = val
    end
  end

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
    vec_X = Vector{Float64}(undef, length_samp)
    vec_samp = Vector{Int}(undef, length_samp)

    for ind_boot in 1:length_boot

      @label LOOP_START
      begin

        sample!(1:length_samp, vec_samp; replace=true)

        for ind_range in 1:length_range
          for ind_samp in 1:length_samp
            vec_X[ind_samp] = X[group[vec_samp[ind_samp]], range_arg[ind_range]]
          end

          if all(vec_X[ind_samp] == vec_X[1] for ind_samp in 1:length_samp)
            @goto LOOP_START
          end
        end
      end

      boot[1:length_samp, ind_boot] = vec_samp
    end
  end

  return out_list
end

wbc_pw = function(ind1, ind2, X, group_list, boot_list, range_arg, vec_X1_list, vec_X2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)

  length_group = length(group_list)
  length_range = length(range_arg)
  length_boot = size(boot_list[1], 2)

  for ind_group in 1:length_group

    @inbounds group = group_list[ind_group]
    @inbounds boot = boot_list[ind_group]

    length_samp = length(group)

    @inbounds vec_X1 = vec_X1_list[ind_group]
    @inbounds vec_X2 = vec_X2_list[ind_group]

    for ind_boot in 1:length_boot

      for ind_samp in 1:length_samp
        @inbounds vec_X1[ind_samp] = X[group[boot[ind_samp, ind_boot]], range_arg[ind1]]
        @inbounds vec_X2[ind_samp] = X[group[boot[ind_samp, ind_boot]], range_arg[ind2]]
      end

      @inbounds vec_cor[ind_boot] = trans == true ? atanh(cor(vec_X1, vec_X2)) : cor(vec_X1, vec_X2)
    end

    sort!(vec_cor)
    @inbounds vec_cor_mat[1:length_boot, ind_group] = vec_cor
  end

  out = 0.
  for ind_boot in 1:length_boot
    for ind_group in 1:length_group
      @inbounds vec_calc[ind_group] = vec_cor_mat[ind_boot, ind_group]
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

  out /= length_boot*length_group

  return out::Float64
end

cor_pw = function(ind1, ind2, X, group_list, range_arg, vec_X1_list, vec_X2_list, vec_calc; lp, trans)

  length_group = length(group_list)

  for ind_group in 1:length_group

    @inbounds group = group_list[ind_group]
    length_samp = length(group)

    @inbounds vec_X1 = vec_X1_list[ind_group]
    @inbounds vec_X2 = vec_X2_list[ind_group]

    for ind_samp in 1:length_samp
      @inbounds vec_X1[ind_samp] = X[group[ind_samp], range_arg[ind1]]
      @inbounds vec_X2[ind_samp] = X[group[ind_samp], range_arg[ind2]]
    end

    @inbounds vec_calc[ind_group] = trans == true ? atanh(cor(vec_X1, vec_X2)) : cor(vec_X1, vec_X2)
  end

  if lp == 1
    vec_calc .-= median(vec_calc)
    vec_calc .= abs.(vec_calc)
  elseif lp == 2
    vec_calc .-= mean(vec_calc)
    vec_calc .= vec_calc.^2
  end
  out = sum(vec_calc) / length_group

  return out::Float64
end

perm_cor = function(length_perm, group_list_arg, X, range_arg = 1:size(X, 2))

  length_range = length(range_arg)
  length_group = length(group_list_arg)

  out_mat = Matrix{Int}(undef, sum(map(x -> length(x), group_list_arg)), length_perm)

  lower = 0
  out_list = map(group_list_arg) do x
    length_samp = length(x)
    out = view(out_mat, (lower+1):(lower+length_samp), 1:length_perm)
    lower += length_samp

    return out
  end

  vec_X_list = map(group_list_arg) do group
    Vector{eltype(X)}(undef, length(group))
  end

  group_list = deepcopy(group_list_arg)
  group_vec = parent(group_list[1])

  for ind_perm in 1:length_perm

    @label LOOP_START
    begin

      shuffle!(group_vec)

      for ind_group in 1:length_group

        group = group_list[ind_group]
        length_samp = length(group)
        vec_X = vec_X_list[ind_group]

        for ind_range in 1:length_range
          for ind_samp in 1:length_samp
            vec_X[ind_samp] = X[group[ind_samp], range_arg[ind_range]]
          end

          if all(vec_X[ind_samp] == vec_X[1] for ind_samp in 1:length_samp)
            @goto LOOP_START
          end
        end
      end
    end

    out_mat[:, ind_perm] = group_vec
  end

  return out_list
end

wbc_perm = function(X, group_list_arg, boot_list, perm_list, range_arg; wp, trans)

  group_vec_arg = parent(group_list_arg[1])
  group_list = deepcopy(group_list_arg)
  group_vec = parent(group_list[1])

  perm_mat = parent(perm_list[1])

  length_perm = size(perm_mat, 2)
  length_total = size(perm_mat, 1)

  length_boot = size(boot_list[1], 2)
  length_range = length(range_arg)
  length_group = length(boot_list)

  vec_X1_list = map(boot_list) do boot
    Vector{Float64}(undef, size(boot, 1))
  end
  vec_X2_list = deepcopy(vec_X1_list)

  vec_cor = Vector{Float64}(undef, length_boot)
  vec_cor_mat = Matrix{Float64}(undef, length_boot, length_group)

  vec_calc = Vector{Float64}(undef, length_group)

  return function(sub)

    sum_orig = wbc_sum_subset(sub, X, group_list_arg, boot_list, range_arg, vec_X1_list, vec_X2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)

    below = 0
    above = 1

    for ind_perm in 1:length_perm

      for ind_total in 1:length_total
        group_vec[ind_total] = group_vec_arg[perm_mat[ind_total, ind_perm]]
      end

      val = wbc_sum_subset(sub, X, group_list, boot_list, range_arg, vec_X1_list, vec_X2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)

      below += val < sum_orig
      above += val >= sum_orig
    end

    return convert(Float64, above) / (convert(Float64, above) + convert(Float64, below))
  end
end

wbc_sum_subset = function(sub, X, group_list, boot_list, range_arg, vec_X1_list, vec_X2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)

  length_sub = length(sub)

  out = 0.

  for ind1 in 1:length_sub
    for ind2 in (ind1+1):length_sub
      val = wbc_pw(sub[ind1], sub[ind2], X, group_list, boot_list, range_arg, vec_X1_list, vec_X2_list, vec_cor, vec_cor_mat, vec_calc; wp, trans)
      out += val
    end
  end

  return out
end
