using Distributed

addprocs(7)

@everywhere using CSV, DataFrames, DelimitedFiles

@everywhere include("BootPerm.jl")
@everywhere using .BootPerm

data_gill = DataFrame!(CSV.File("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv"))

DATA_MAT = convert(Matrix, select(data_gill, Not([:sample, :location])))
(GROUP_LIST, GROUP_VEC) = split_groups(data_gill[:, :location])
PROTEINS = names(select(data_gill, Not([:sample, :location])))
RANGE = 1:size(DATA_MAT, 2)

SUB_SIZE = 100
KEEP = 100
REPS = 200
MAX_ITER = 10^5

@everywhere include("OptimSubset.jl")
@everywhere include("CorSep.jl")

@everywhere using .OptimSubset, .CorSep

max_list_total = let

    (COR_MAT_LIST, VEC_COR_LIST) = cor_sep_alloc(SUB_SIZE, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE;
        trans = true)

    @everywhere workers() COR_SEP = CorSepCall($COR_MAT_LIST, deepcopy($VEC_COR_LIST))

    optim_subset_iter(x -> COR_SEP(x), SUB_SIZE, length(RANGE);
        reps = REPS,
        type = "max",
        keep = KEEP,
        max_iter = MAX_ITER,
        worker_ids = workers())
end

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

CSV.write("ind_max_total.csv", hcat(DataFrame(val = map(x->x[1], max_list_total)), DataFrame(hcat(map(x->x[2], max_list_total)...)')))
CSV.write("ind_max_rosaria.csv", hcat(DataFrame(val = map(x->x[1], max_list_rosaria)), DataFrame(hcat(map(x->x[2], max_list_rosaria)...)')))
CSV.write("ind_max_bodega.csv", hcat(DataFrame(val = map(x->x[1], max_list_bodega)), DataFrame(hcat(map(x->x[2], max_list_bodega)...)')))
CSV.write("ind_max_solano.csv", hcat(DataFrame(val = map(x->x[1], max_list_solano)), DataFrame(hcat(map(x->x[2], max_list_solano)...)')))
CSV.write("ind_max_westchester.csv", hcat(DataFrame(val = map(x->x[1], max_list_westchester)), DataFrame(hcat(map(x->x[2], max_list_westchester)...)')))
