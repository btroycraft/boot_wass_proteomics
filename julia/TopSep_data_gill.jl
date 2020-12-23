using Distributed

num_workers = 0

addprocs(num_workers)

@everywhere using CSV, DataFrames, DelimitedFiles

@everywhere include("BootPerm.jl")
@everywhere using .BootPerm

data_gill = DataFrame!(CSV.File("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv"))

DATA_MAT = convert(Matrix, select(data_gill, Not([:sample, :location])))
(GROUP_LIST, GROUP_VEC) = split_groups(data_gill[:, :location])
PROTEINS = names(select(data_gill, Not([:sample, :location])))
RANGE = 1:size(DATA_MAT, 2)

SUB_SIZE = 20
KEEP = 100
REPS = 3
MAX_ITER = 100

@everywhere include("OptimSubset.jl")
@everywhere include("TopSep.jl")

@everywhere using .OptimSubset, .TopSep

max_list_total = let

    (DIST_MAT_LIST, DIST_SUB_MAT_LIST) = top_sep_alloc(SUB_SIZE, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE)

    @everywhere workers() TOP_SEP = TopSepCall(0, $DIST_MAT_LIST, deepcopy($DIST_SUB_MAT_LIST))

    optim_subset_iter(x -> TOP_SEP(x), SUB_SIZE, length(RANGE);
        reps = REPS,
        type = "max",
        keep = KEEP,
        max_iter = MAX_ITER,
        worker_ids = workers())
end

CSV.write("ind_max_total.csv", hcat(DataFrame(val = map(x->x[1], max_list_total)), DataFrame(hcat(map(x->x[2], max_list_total)...)')))

max_list_rosaria = let

    (COR_MAT_LIST, VEC_COR_LIST) = cor_sep_group_alloc(SUB_SIZE, 1, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE;
        trans = true)

    @everywhere workers() COR_SEP = CorSepGroupCall($COR_MAT_LIST, deepcopy($VEC_COR_LIST))

    optim_subset_iter(x -> COR_SEP(x), SUB_SIZE, length(RANGE);
        reps = REPS,
        type = "max",
        keep = KEEP,
        max_iter = MAX_ITER,
        worker_ids = workers())
end

CSV.write("ind_max_rosaria.csv", hcat(DataFrame(val = map(x->x[1], max_list_rosaria)), DataFrame(hcat(map(x->x[2], max_list_rosaria)...)')))

max_list_bodega = let

    (COR_MAT_LIST, VEC_COR_LIST) = cor_sep_group_alloc(SUB_SIZE, 2, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE;
        trans = true)

    @everywhere workers() COR_SEP = CorSepGroupCall($COR_MAT_LIST, deepcopy($VEC_COR_LIST))

    optim_subset_iter(x -> COR_SEP(x), SUB_SIZE, length(RANGE);
        reps = REPS,
        type = "max",
        keep = KEEP,
        max_iter = MAX_ITER,
        worker_ids = workers())
end

CSV.write("ind_max_bodega.csv", hcat(DataFrame(val = map(x->x[1], max_list_bodega)), DataFrame(hcat(map(x->x[2], max_list_bodega)...)')))

max_list_solano = let

    (COR_MAT_LIST, VEC_COR_LIST) = cor_sep_group_alloc(SUB_SIZE, 3, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE;
        trans = true)

    @everywhere workers() COR_SEP = CorSepGroupCall($COR_MAT_LIST, deepcopy($VEC_COR_LIST))

    optim_subset_iter(x -> COR_SEP(x), SUB_SIZE, length(RANGE);
        reps = REPS,
        type = "max",
        keep = KEEP,
        max_iter = MAX_ITER,
        worker_ids = workers())
end

CSV.write("ind_max_solano.csv", hcat(DataFrame(val = map(x->x[1], max_list_solano)), DataFrame(hcat(map(x->x[2], max_list_solano)...)')))

max_list_westchester = let

    (COR_MAT_LIST, VEC_COR_LIST) = cor_sep_group_alloc(SUB_SIZE, 4, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE;
        trans = true)

    @everywhere workers() COR_SEP = CorSepGroupCall($COR_MAT_LIST, deepcopy($VEC_COR_LIST))

    optim_subset_iter(x -> COR_SEP(x), SUB_SIZE, length(RANGE);
        reps = REPS,
        type = "max",
        keep = KEEP,
        max_iter = MAX_ITER,
        worker_ids = workers())
end

CSV.write("ind_max_westchester.csv", hcat(DataFrame(val = map(x->x[1], max_list_westchester)), DataFrame(hcat(map(x->x[2], max_list_westchester)...)')))


rmprocs(workers())
