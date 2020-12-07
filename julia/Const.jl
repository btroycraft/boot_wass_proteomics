using CSV, DataFrames, Distributed, ParallelDataTransfer, DelimitedFiles

NUM_PROCS = 7
addprocs(NUM_PROCS-1)

include("BootPerm.jl")

using .BootPerm

data_gill = DataFrame!(CSV.File("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv"))

DATA_MAT = convert(Matrix, select(data_gill, Not([:sample, :location])))
(GROUP_LIST, GROUP_VEC) = split_groups(data_gill[:, :location])
PROTEINS = names(select(data_gill, Not([:sample, :location])))
RANGE = 1:size(DATA_MAT, 2)

SUB_SIZE = 100
KEEP = 100
REPS_TOTAL = 21
REPS = convert(Int, ceil(REPS_TOTAL/NUM_PROCS))
MAX_ITER = 10^5

@everywhere begin
    include("OptimSubset.jl")
    include("CorSep.jl")
end

@everywhere using .OptimSubset, .CorSep

@everywhere workers() begin

    SUB_SIZE = $SUB_SIZE
    KEEP = $KEEP
    REPS = $REPS
    MAX_ITER = $MAX_ITER

    GROUP_LIST = $GROUP_LIST
    GROUP_VEC = $GROUP_VEC
    DATA_MAT = $DATA_MAT
    RANGE = $RANGE
end

@everywhere (COR_MAT_LIST, VEC_COR_LIST) = cor_sep_alloc(SUB_SIZE, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE;
        trans = true)
@everywhere workers() COR_MAT_LIST = $COR_MAT_LIST
@everywhere max_list_t = optim_subset_iter(SUB_SIZE, length(RANGE);
    keep = KEEP,
    max_iter = MAX_ITER,
    reps = REPS) do sub
        cor_sep!(sub, GROUP_LIST, GROUP_VEC, DATA_MAT, RANGE, COR_MAT_LIST, VEC_COR_LIST; trans = true)
    end

temp = [getfrom(id, :max_list_t) for id in procs()]
rmprocs(workers())
max_list = vcat(temp...)
max_list = sort!(max_list;
    by = x -> x[1],
    rev = true)
max_list = max_list[1:KEEP]

CSV.write("upd.csv", hcat(DataFrame(val = map(x->x[1], max_list)), DataFrame(hcat(map(x->x[2], max_list)...)')))
