using CSV
using DataFrames
using Plots
using Distributed

include("OptimSubset.jl")
include("WassBootCor.jl")
include("SumSubset.jl")

data_gill = DataFrame!(CSV.File("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv"))
data_gill_mat = convert(Matrix, select(data_gill, Not([:sample, :location])))

data_indices = DataFrame!(CSV.File("prot_indices.csv"))

index_range = data_indices[:, :indices]
wass_boot_reps = 1000
perm_reps = 10000
keep_top = 1000
keep_bottom = keep_top
subset_size = 5
group_list = split_groups(repeat(1:4; inner = 24))
boot_list = boot_reps_cor(wass_boot_reps, group_list, data_gill_mat, index_range)
perm_list = perm_cor(perm_reps, group_list, data_gill_mat, index_range)

wbc_mat_w2 = wbc_mat(data_gill_mat, group_list, boot_list, index_range; wp = 2, trans = true)
wbc_sub = sum_subset(wbc_mat_w2)

@time max_subset(wbc_sub, 5, length(index_range);
    keep = keep_top)
best_subs = max_subset_iter(wbc_sub, 5, length(index_range);
    reps = 20,
    keep = keep_top)

(x1, y1) = dist_subset_rand(wbc_sub, 5, length(index_range);
    lower = 0,
    upper = 1,
    length_bin = 1000,
    reps = 1000000)
plot(collect(y1[1:length(x1)]), x1)

perm_func = wbc_perm(data_gill_mat, group_list, boot_list, perm_list, index_range; wp = 2, trans = true)

perm_func(best_subs[2][:, 1])

addprocs(5)

(x2, y2) = dist_subset_rand(perm_func, 5, length(index_range);
    lower = 0,
    upper = 1,
    length_bin = 100,
    reps = 200)
plot(collect(y2[1:length(x2)]), x2)
