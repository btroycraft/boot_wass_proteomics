using CSV, DataFrames, Distributed

include("BootPerm.jl")
include("SumSubset.jl")
include("OptimSubset.jl")
include("WassBootCor.jl")

WORKER_IDS = addprocs(6)
using .BootPerm, .SumSubset, .OptimSubset, .WassBootCor

@everywhere WORKER_IDS begin

    include("BootPerm.jl")
    include("SumSubset.jl")
    include("OptimSubset.jl")
    include("WassBootCor.jl")

    using .BootPerm, .WassBootCor, .OptimSubset, .SumSubset
end

data_gill = DataFrame!(CSV.File("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv"))

DATA_MAT = convert(Matrix, select(data_gill, Not([:sample, :location])))
(GROUP_LIST, GROUP_VEC) = split_groups(data_gill[:, :location])
RANGE = DataFrame!(CSV.File("prot_indices.csv"))[:, :indices]
RANGE_NAMES = names(select(data_gill, Not([:sample, :location])))[RANGE]

BOOT_REPS_GEN = 10000
BOOT_MAT = boot_reps(10*BOOT_REPS_GEN, GROUP_LIST)

PERM_REPS = 1000
PERM_MAT = perm_cor(PERM_REPS, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE)
BOOT_REPS = 10

KEEP = 100

TRAIN_REPS = 1000
TRAIN_SIZE = 16
TEST_SIZE = 8
(TRAIN_LIST, TEST, TEST_LAB, TRAIN_MAT) = split_train(TRAIN_REPS, TRAIN_SIZE, TEST_SIZE, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE)

BOOT_REPS_TRAIN = 1000
BOOT_MAT_TRAIN = boot_reps_train(10*BOOT_REPS_TRAIN, TRAIN_LIST, TEST)

@everywhere WORKER_IDS begin
    DATA_MAT = $DATA_MAT
    GROUP_LIST = $GROUP_LIST
    GROUP_VEC = $GROUP_VEC
    RANGE = $RANGE
    RANGE_NAMES = $RANGE_NAMES
    BOOT_REPS_GEN = $BOOT_REPS_GEN
    BOOT_MAT = $BOOT_MAT
    PERM_REPS = $PERM_REPS
    PERM_MAT = $PERM_MAT
    BOOT_REPS = $BOOT_REPS
    KEEP = $KEEP
    TRAIN_REPS = $TRAIN_REPS
    TRAIN_SIZE = $TRAIN_SIZE
    TEST_SIZE = $TEST_SIZE
    TRAIN_LIST = $TRAIN_LIST
    TEST = $TEST
    TEST_LAB = $TEST_LAB
    TRAIN_MAT = $TRAIN_MAT
    BOOT_REPS_TRAIN = $BOOT_REPS_TRAIN
    BOOT_MAT_TRAIN = $BOOT_MAT_TRAIN
end

wbc_mat_w1 = wbc_mat(BOOT_REPS, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE, BOOT_MAT; wp = 1, trans = false)
wbc_mat_w1t = wbc_mat(BOOT_REPS, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE, BOOT_MAT; wp = 1, trans = true)
wbc_mat_w2 = wbc_mat(BOOT_REPS, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE, BOOT_MAT; wp = 2, trans = false)
wbc_mat_w2t = wbc_mat(BOOT_REPS, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE, BOOT_MAT; wp = 2, trans = true)

wbc_sub_w1 = sum_subset_close(wbc_mat_w1)
wbc_sub_w1t = sum_subset_close(wbc_mat_w1t)
wbc_sub_w2 = sum_subset_close(wbc_mat_w2)
wbc_sub_w2t = sum_subset_close(wbc_mat_w2t)

@everywhere WORKER_IDS begin
    wbc_mat_w1 = $wbc_mat_w1
    wbc_mat_w1t = $wbc_mat_w1t
    wbc_mat_w2 = $wbc_mat_w2
    wbc_mat_w2t = $wbc_mat_w2t

    wbc_sub_w1 = sum_subset_close(wbc_mat_w1)
    wbc_sub_w1t = sum_subset_close(wbc_mat_w1t)
    wbc_sub_w2 = sum_subset_close(wbc_mat_w2)
    wbc_sub_w2t = sum_subset_close(wbc_mat_w2t)
end

cor_mat_l1 = cor_mat(GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE; lp = 1, trans = false)
cor_mat_l1t = cor_mat(GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE; lp = 1, trans = true)
cor_mat_l2 = cor_mat(GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE; lp = 2, trans = false)
cor_mat_l2t = cor_mat(GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE; lp = 2, trans = true)

cor_sub_l1 = sum_subset_close(cor_mat_l1)
cor_sub_l1t = sum_subset_close(cor_mat_l1t)
cor_sub_l2 = sum_subset_close(cor_mat_l2)
cor_sub_l2t = sum_subset_close(cor_mat_l2t)

@everywhere WORKER_IDS begin
    cor_mat_l1 = $cor_mat_l1
    cor_mat_l1t = $cor_mat_l1t
    cor_mat_l2 = $cor_mat_l2
    cor_mat_l2t = $cor_mat_l2t

    cor_sub_l1 = sum_subset_close(cor_mat_l1)
    cor_sub_l1t = sum_subset_close(cor_mat_l1t)
    cor_sub_l2 = sum_subset_close(cor_mat_l2)
    cor_sub_l2t = sum_subset_close(cor_mat_l2t)
end

wbc_sub_max_w2t = max_subset_iter(wbc_sub_w2t, 5, length(RANGE); keep = KEEP, reps = 50)
cor_sub_max_l2t = max_subset_iter(cor_sub_l2t, 5, length(RANGE); keep = KEEP, reps = 50)

@everywhere WORKER_IDS begin
    wbc_sub_max_w2t = $wbc_sub_max_w2t
    cor_sub_max_l2t = $cor_sub_max_l2t
end

@everywhere run_func_w2t = function(ind)
    sub = wbc_sub_max_w2t[2][:, ind]
    return wbc_cv(sub, BOOT_REPS_TRAIN, TRAIN_LIST, TEST, TEST_LAB, TRAIN_MAT, GROUP_VEC, DATA_MAT, RANGE, BOOT_MAT_TRAIN;
        wp = 2,
        trans = true)
end

wbc_cv_w2t = pmap(run_func_w2t, WorkerPool(WORKER_IDS), 1:KEEP)
cor_cv_l2t = pmap(WorkerPool(WORKER_IDS), 1:KEEP) do ind
    sub = cor_sub_max_l2t[2][:, ind]
    return cor_cv(sub, TRAIN_LIST, TEST, TEST_LAB, TRAIN_MAT, GROUP_VEC, DATA_MAT, RANGE;
        lp = 2,
        trans = true)
end


end
