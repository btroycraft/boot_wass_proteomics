using Distributed

num_workers = 6

addprocs(num_workers)

@everywhere using CSV, DataFrames, DelimitedFiles

@everywhere include("BootPerm.jl")
@everywhere using .BootPerm

data_gill = DataFrame!(CSV.File("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv"))
ind_max = convert(Vector{Int}, select(DataFrame!(CSV.File("ind_max_total.csv")), Not(:val))[1, :])

DATA_MAT = convert(Matrix, select(data_gill, Not([:sample, :location])))
(GROUP_LIST, GROUP_VEC) = split_groups(data_gill[:, :location])
PROTEINS = names(select(data_gill, Not([:sample, :location])))
RANGE = ind_max

SUB_SIZE = 20
KEEP = 100
REPS = 200
MAX_ITER = 10^5

@everywhere include("OptimSubset.jl")
@everywhere include("TopSep.jl")

@everywhere using .OptimSubset, .TopSep

max_list_total_betti0 = let

    (DIST_MAT_LIST, DIST_SUB_MAT_LIST) = top_sep_alloc(SUB_SIZE, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE)

    @everywhere workers() TOP_SEP = TopSepCall(0, $DIST_MAT_LIST, deepcopy($DIST_SUB_MAT_LIST))

    optim_subset_iter(x -> TOP_SEP(x), SUB_SIZE, length(RANGE);
        reps = REPS,
        type = "max",
        keep = KEEP,
        max_iter = MAX_ITER,
        worker_ids = workers())
end

CSV.write("ind_top_total_betti0.csv", hcat(DataFrame(val = map(x->x[1], max_list_total_betti0)), DataFrame(hcat(map(x->x[2], max_list_total_betti0)...)')))

max_list_total_betti1 = let

    (DIST_MAT_LIST, DIST_SUB_MAT_LIST) = top_sep_alloc(SUB_SIZE, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE)

    @everywhere workers() TOP_SEP = TopSepCall(1, $DIST_MAT_LIST, deepcopy($DIST_SUB_MAT_LIST))

    optim_subset_iter(x -> TOP_SEP(x), SUB_SIZE, length(RANGE);
        reps = REPS,
        type = "max",
        keep = KEEP,
        max_iter = MAX_ITER,
        worker_ids = workers())
end

CSV.write("ind_top_total_betti1.csv", hcat(DataFrame(val = map(x->x[1], max_list_total_betti1)), DataFrame(hcat(map(x->x[2], max_list_total_betti1)...)')))

rmprocs(workers())
