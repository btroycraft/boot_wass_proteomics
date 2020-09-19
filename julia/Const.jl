using CSV
using DataFrames: DataFrame!, Not, select

include("BootPerm.jl")
include("WassBootCor.jl")
include("OptimSubset.jl")

data_gill = DataFrame!(CSV.File("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv"))

DATA_MAT = convert(Matrix, select(data_gill, Not([:sample, :location])))
(GROUP_LIST, GROUP_VEC) = split_groups(data_gill[:, :location])
RANGE = DataFrame!(CSV.File("prot_indices.csv"))[:, :indices]
RANGE_NAMES = names(select(data_gill, Not([:sample, :location])))[RANGE]

boot_reps = 100
perm_reps = 100

BOOT_MAT = boot_reps_cor(boot_reps, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE)
PERM_MAT = perm_cor(perm_reps, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE)

wbc_mat_w1 = wbc_mat(GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE, BOOT_MAT; wp = 1, trans = true)

wbc_sub_w1 = sum_subset_close(wbc_mat_w1)

min_subset(wbc_sub_w1, 4, length(RANGE);
    keep = 1000)

dist_subset_rand(wbc_sub_w1, 4, length(RANGE);
    lower = 0,
    upper = 1,
    num_bins = 100,
    reps = 10000)
