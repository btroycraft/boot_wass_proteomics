module BootPerm

import StatsBase, Random, Statistics

export split_groups, boot_reps, perm_cor, split_train, boot_reps_train

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

boot_reps = function(boot_reps, group_list)

  out = Matrix{Int}(undef, sum([length(group) for group in group_list]), boot_reps)

  vec_samp_list = [Vector{Int}(undef, length(group)) for group in group_list]

  for ind_boot in 1:boot_reps
    for ind_group in 1:length(group_list)

      @inbounds group = group_list[ind_group]
      @inbounds vec_samp = vec_samp_list[ind_group]

      @label LOOP_START

      StatsBase.sample!(group, vec_samp;
        replace = true)

      all(vec_samp[ind] == vec_samp[1] for ind in 2:length(group)) && @goto LOOP_START

      for ind_samp in 1:length(group)
        @inbounds out[group[ind_samp], ind_boot] = vec_samp[ind_samp]
      end
    end
  end

  return out
end

perm_cor = function(perm_reps, group_list, group_vec, data_mat, range)

  out = Matrix{Int}(undef, length(group_vec), perm_reps)

  vec_x_list = [Vector{Float64}(undef, length(group)) for group in group_list]
  vec_samp = collect(1:length(group_vec))

  for ind_perm in 1:perm_reps

    perm_cor_!(group_list, group_vec, data_mat, range, vec_x_list, vec_samp)

    for ind_samp in 1:length(group_vec)
      @inbounds out[ind_samp, ind_perm] = vec_samp[ind_samp]
    end
  end

  return out
end

perm_cor_! =  function(group_list, group_vec, data_mat, range, vec_x_list, vec_samp)

  @label LOOP_START
  Random.shuffle!(vec_samp)

  for ind_group in 1:length(group_list)

    @inbounds group = group_list[ind_group]
    @inbounds vec_x = vec_x_list[ind_group]

    for ind_range in 1:length(range)
      for ind_samp in 1:length(group)
        @inbounds vec_x[ind_samp] = data_mat[group_vec[vec_samp[group[ind_samp]]], range[ind_range]]
      end

      iszero(Statistics.var(vec_x)) && @goto LOOP_START
    end
  end
end

split_train = function(train_reps, train_size, test_size, group_list, group_vec, data_mat, range)

  out_mat = Matrix{Int}(undef, train_size*length(group_list) + test_size, train_reps)
  out_train = [((ind-1) * train_size+1):((ind-1)*train_size + train_size) for ind in 1:length(group_list)]
  out_test = (length(group_list)*train_size + 1):(length(group_list)*train_size + test_size)
  out_lab = Vector{Int}(undef, train_reps)

  vec_samp_list = [collect(group) for group in group_list]
  vec_samp_train_list = [Vector{Int}(undef, train_size) for _ in 1:length(group_list)]
  vec_samp_test = Vector{Int}(undef, test_size)

  vec_x_train_list = [Vector{Float64}(undef, train_size) for _ in 1:length(group_list)]
  vec_x_test = Vector{Float64}(undef, test_size)

  for ind_train in 1:train_reps

    out_lab[ind_train] = split_train_!(group_list, group_vec, train_size, test_size, data_mat, range, vec_samp_list, vec_samp_train_list, vec_samp_test, vec_x_train_list, vec_x_test)

    for ind_samp in 1:test_size
      out_mat[out_test[ind_samp], ind_train] = vec_samp_test[ind_samp]
    end

    for ind_group in 1:length(group_list)
      train = out_train[ind_group]
      vec_samp_train = vec_samp_train_list[ind_group]

      for ind_samp in 1:train_size
        out_mat[train[ind_samp], ind_train] = vec_samp_train[ind_samp]
      end
    end
  end

  return out_train, out_test, out_lab, out_mat
end

split_train_! =  function(group_list, group_vec, train_size, test_size, data_mat, range, vec_samp_list, vec_samp_train_list, vec_samp_test, vec_x_train_list, vec_x_test)

  @label LOOP_START

  for vec_samp in vec_samp_list
    Random.shuffle!(vec_samp)
  end

  lab_test = rand(1:length(group_list))
  vec_samp = vec_samp_list[lab_test]

  for ind_samp in 1:test_size
    vec_samp_test[ind_samp] = vec_samp[train_size + ind_samp]
  end

  for ind_range in 1:length(range)
    for ind_samp in 1:test_size
      vec_x_test[ind_samp] = data_mat[group_vec[vec_samp_test[ind_samp]], range[ind_range]]
    end

    iszero(Statistics.var(vec_x_test)) && @goto LOOP_START
  end

  for ind_group in 1:length(group_list)

    vec_samp = vec_samp_list[ind_group]
    vec_samp_train = vec_samp_train_list[ind_group]
    vec_x_train = vec_x_train_list[ind_group]

    for ind_samp in 1:train_size
      vec_samp_train[ind_samp] = vec_samp[ind_samp]
    end

    for ind_range in 1:length(range)
      for ind_samp in 1:train_size
        vec_x_train[ind_samp] = data_mat[group_vec[vec_samp_train[ind_samp]], range[ind_range]]
      end

      iszero(Statistics.var(vec_x_train)) && @goto LOOP_START
    end
  end

  return lab_test
end

boot_reps_train = function(boot_reps, train_list, test)

  out = Matrix{Int}(undef, sum([length(train) for train in train_list]) + length(test), boot_reps)

  vec_samp_train_list = [Vector{Int}(undef, length(train)) for train in train_list]
  vec_samp_test = Vector{Int}(undef, length(test))

  for ind_boot in 1:boot_reps

    @label LOOP_START_TEST

    StatsBase.sample!(test, vec_samp_test;
      replace = true)

    all(vec_samp_test[ind] == vec_samp_test[1] for ind in 2:length(test)) && @goto LOOP_START_TEST

    for ind_samp in 1:length(test)
      out[test[ind_samp], ind_boot] = vec_samp_test[ind_samp]
    end

    for ind_group in 1:length(train_list)

      train = train_list[ind_group]
      vec_samp_train = vec_samp_train_list[ind_group]

      @label LOOP_START_TRAIN

      StatsBase.sample!(train, vec_samp_train;
        replace = true)

      all(vec_samp_train[ind] == vec_samp_train[1] for ind in 2:length(train)) && @goto LOOP_START_TRAIN

      for ind_samp in 1:length(train)
        out[train[ind_samp], ind_boot] = vec_samp_train[ind_samp]
      end
    end
  end

  return out
end

end
