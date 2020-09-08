using CSV
using DataFrames

include("OptimSubset.jl")
using .OptimSubset

include("WassBootCor.jl")
using .WassBootCor

include("SumSubset.jl")
using .SumSubset

data_gill = DataFrame!(CSV.File("cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv"))

data_gill_mat = convert(Matrix, select(data_gill, Not([:sample, :location])))

wass_boot_w2_sep_mat = wass_boot_cor_sep(data_gill_mat; groups = repeat(1:4; inner = 24), boot_reps = 1000, wp = "w2", transformed = true)
wass_boot_w1_sep_mat = wass_boot_cor_sep(data_gill_mat; groups = repeat(1:4; inner = 24), boot_reps = 1000, wp = "w1", transformed = true)

cor_l2_sep_mat = cor_sep(data_gill_mat; groups = repeat(1:4; inner = 24), lp = "l2", transformed = true)
cor_l1_sep_mat = cor_sep(data_gill_mat; groups = repeat(1:4; inner = 24), lp = "l1", transformed = true)

wass_boot_w2_subset = sum_subset(wass_boot_w2_sep_mat)
wass_boot_w2_sep_subset = sum_subset(wass_boot_w2_sep_mat)
wass_boot_w1_sep_subset = sum_subset(wass_boot_w1_sep_mat)
cor_12_sep_subset = sum_subset(cor_l2_sep_mat)
cor_11_sep_subset = sum_subset(cor_l1_sep_mat)
