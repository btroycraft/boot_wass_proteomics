using CSV, DataFrames

include("CorSep.jl")
include("OptimSubset.jl")
include("BootPerm.jl")

using .CorSep, .OptimSubset, .BootPerm

data_gill = DataFrame!(CSV.File("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv"))

DATA_MAT = convert(Matrix, select(data_gill, Not([:sample, :location])))
(GROUP_LIST, GROUP_VEC) = split_groups(data_gill[:, :location])
RANGE = DataFrame!(CSV.File("prot_indices.csv"))[:, :indices]
RANGE_NAMES = names(select(data_gill, Not([:sample, :location])))[RANGE]

PERM_REPS = 1000
PERM_MAT = perm_cor(PERM_REPS, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE)

KEEP = 100

TRAIN_REPS = 1000
TRAIN_SIZE = 16
TEST_SIZE = 8
(TRAIN_LIST, TEST, TEST_LAB, TRAIN_MAT) = split_train(TRAIN_REPS, TRAIN_SIZE, TEST_SIZE, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE)

BOOT_REPS_TRAIN = 1000
BOOT_MAT_TRAIN = boot_reps_train(100*BOOT_REPS_TRAIN, TRAIN_LIST, TEST)

cor_sep_sub = cor_sep_close(100, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE)

max_subset_iter(cor_sep_sub, 50, length(RANGE); reps = 20, keep = 20, iter = 10, show = true)
