using CSV, DataFrames, Distributed, ParallelDataTransfer

addprocs(4)

include("BootPerm.jl")

using .BootPerm

data_gill = DataFrame!(CSV.File("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv"))

DATA_MAT = convert(Matrix, select(data_gill, Not([:sample, :location])))
(GROUP_LIST, GROUP_VEC) = split_groups(data_gill[:, :location])
PROTEINS = names(select(data_gill, Not([:sample, :location])))

SUB_SIZE = 100
KEEP = 20

@everywhere begin
    include("OptimSubset.jl")
    include("CorSep.jl")
end

@everywhere using .OptimSubset, .CorSep

@everywhere workers() begin
    cor_sep_sub = cor_sep_close($SUB_SIZE, $GROUP_LIST, $GROUP_VEC, $DATA_MAT;
        trans = true)

    max_list_t = optim_subset_iter(cor_sep_sub, $SUB_SIZE, size($DATA_MAT, 2);
        reps = 1,
        keep = KEEP,
        max_iter = 1000,
        type = "max",
        progress = true)
end

max_list = [getfrom(id, :max_list_t) for id in workers()]
keep_out = min(KEEP, sum(map(length, max_list))]
max_list = vcat(max_list...)
max_list = sort!(max_list;
    by = x -> x[1],
    rev = true)
max_list = max_list[1:keep_out]
