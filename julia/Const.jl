using CSV, DataFrames, Distributed
Distributed.worker_timeout() = 300

WORKER_IDS = addprocs(4)

@everywhere begin
    include("OptimSubset.jl")
    include("BootPerm.jl")
    include("CorSep.jl")
end

@everywhere using .OptimSubset, .BootPerm, .CorSep

data_gill = DataFrame!(CSV.File("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv"))

DATA_MAT = convert(Matrix, select(data_gill, Not([:sample, :location])))
(GROUP_LIST, GROUP_VEC) = split_groups(data_gill[:, :location])
PROTEINS = names(select(data_gill, Not([:sample, :location])))

SUB_SIZE = 100

@everywhere cor_sep_sub = cor_sep_close($SUB_SIZE, $GROUP_LIST, $GROUP_VEC, $DATA_MAT;
    trans = true)

max_list = optim_subset_iter(x -> cor_sep_sub(x), SUB_SIZE, size(DATA_MAT, 2);
    reps = 100,
    keep = 20,
    max_iter = 300,
    type = "max",
    worker_ids = workers())
max_list
