setwd("C:/Users/jamsh/OneDrive/Desktop/RTG_proteomics_Local")
source("selecting_proteins.R")
library(readxl)
library(Matrix)
library(igraph)
data <- read_excel("C:/Users/jamsh/OneDrive/Desktop/RTG_proteomics_Local/cDIA_MXLSAKBH-Exp1-2-3-4_Gill_with Volcanoplots_SelectedSubset (2).xlsx")



sim_id = c()
for(i in 1:200){
  num = as.character(i)
  a = paste("sim",num, sep = "")
  sim_id = c(sim_id, rep(a,24))
}


vec_cor <- function(df){
  temp = cor(df[1:5])
  temp[which(lower.tri(temp))]
}

wass_by_pair <- function(pair){
  corr_mid <- apply(cbind(pair$ak, pair$bh, pair$ls, pair$m), 1, median)
  wass_ak <- mean(abs(corr_mid - pair$ak))
  wass_bh <- mean(abs(corr_mid - pair$bh))
  wass_ls <- mean(abs(corr_mid - pair$ls))
  wass_m <- mean(abs(corr_mid - pair$m))
  
  return(wass_ak + wass_bh + wass_ls + wass_m)
}


set.seed(2020)

#Grab subset of 5 proteins
sub_prots = sig_prots_special[1:5]

#convert data to numeric
labeled_data = data.frame(lapply(labeled_data, function(x) as.numeric(as.character(x))))





rep = 200
sample_size = 24
total_sim = rep*sample_size

#Bootstrap 24 fish from each environment (with replacement)
alaska_boot = labeled_data[sample(1:24, replace = TRUE,size = total_sim),]
bodega_boot = labeled_data[sample(25:48, replace = TRUE,size = total_sim),]
solano_boot = labeled_data[sample(49:72, replace = TRUE, size = total_sim),]
mexico_boot = labeled_data[sample(73:96, replace = TRUE, size = total_sim),]



#Making data frames of 5 proteins, with sim_id variable
df_ak = cbind(data.frame(alaska_boot[,sub_prots]),sim_id)
df_bh = cbind(data.frame(bodega_boot[,sub_prots]),sim_id)
df_ls = cbind(data.frame(solano_boot[,sub_prots]),sim_id)
df_m = cbind(data.frame(mexico_boot[,sub_prots]),sim_id)

#Find correlations for each of the environments
ak_corrs = by(df_ak[1:5],df_ak$sim_id, FUN = vec_cor)

bh_corrs = by(df_bh[1:5],df_bh$sim_id, FUN = vec_cor)

ls_corrs = by(df_ls[1:5],df_ls$sim_id, FUN = vec_cor)

m_corrs = by(df_m[1:5],df_m$sim_id, FUN = vec_cor)


comb.pairs = combn(sub_prots,2)
name.pairs = apply(comb.pairs, 2, paste, collapse="-")

#Reshape into data frame 
ak_corrs = as.data.frame(do.call("rbind",ak_corrs))
bh_corrs = as.data.frame(do.call("rbind",bh_corrs))
ls_corrs = as.data.frame(do.call("rbind",ls_corrs))
m_corrs = as.data.frame(do.call("rbind",m_corrs))


#Make into list of pairs with 200 correlations for each env
corr.distributions.by.pair <- list(
  list(ak = ak_corrs$V1,bh = bh_corrs$V1,ls = ls_corrs$V1,m = m_corrs$V1),
  list(ak = ak_corrs$V2, bh =bh_corrs$V2,ls = ls_corrs$V2,m = m_corrs$V2),
  list(ak = ak_corrs$V3,bh =bh_corrs$V3,ls = ls_corrs$V3,m = m_corrs$V3),
  list(ak = ak_corrs$V4,bh =bh_corrs$V4,ls = ls_corrs$V4,m = m_corrs$V4),
  list(ak = ak_corrs$V5,bh =bh_corrs$V5,ls = ls_corrs$V5,m = m_corrs$V5),
  list(ak = ak_corrs$V6,bh =bh_corrs$V6,ls = ls_corrs$V6,m = m_corrs$V6),
  list(ak = ak_corrs$V7,bh =bh_corrs$V7,ls = ls_corrs$V7,m = m_corrs$V7),
  list(ak = ak_corrs$V8,bh =bh_corrs$V8,ls = ls_corrs$V8,m = m_corrs$V8),
  list(ak = ak_corrs$V9,bh =bh_corrs$V9,ls = ls_corrs$V9,m = m_corrs$V9),
  list(ak = ak_corrs$V10,bh =bh_corrs$V10,ls = ls_corrs$V10,m = m_corrs$V10)
)

#Assign the names
names(corr.distributions.by.pair) <- name.pairs

#Compute Wass Distance
wass_dist = as.vector(unlist(lapply(corr.distributions.by.pair, FUN = wass_by_pair)))


edge_vec = unlist(strsplit(name.pairs, split = "-"))

#Create network
net <- graph(edges = edge_vec, directed = F)

E(net)$width = (wass_dist)^2*10
E(net)$type = ifelse(wass_dist > cut.off, "hyperlink", "none")

#set cutoff for green edge
cut.off = 0.8

plot(net, edge.color=c("slategrey", "#008c3c")[(E(net)$type=="hyperlink")+1])

