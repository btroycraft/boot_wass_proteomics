setwd("C:/Users/jamsh/OneDrive/Desktop/RTG_proteomics_Local")
source("selecting_proteins.R")
source('wass.R')
source('gen_boot.R')
source('sum_combn_mat.R')
library(readxl)
library(Matrix)
library(igraph)
library(ggplot2)

get.labs <- function(ind, n){
  triag <- (0:(n-1))*(1:n)/2
  ind2 <- sum(triag < k)+1
  ind1 <- k - triag[ind2-1]
  
  return(c(ind1, ind2))
}

boot_reps <- 20

data <- read.table('cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv', header=TRUE, sep=',')
#data <- as.data.frame(readxl::read_excel("./cDIA_MXLSAKBH-Exp1-2-3-4_Gill_with Volcanoplots_SelectedSubset.xlsx"))
data <- data[, 1:200]

protein_names <- colnames(data)
protein_names <- protein_names[!(protein_names %in% c('location', 'sample'))]

pops <- lapply(split(data[, !(colnames(data) %in% c('sample', 'location'))], data$location), function(df){
  df <- as.matrix(df)
  return( df )
})

pop_boots <- lapply(pops, function(pop){
  
  t(replicate(boot_reps, {
    repeat({
      X_boot <- gen.boot(pop)
      cor_boot <- try(cor(X_boot), silent=TRUE)
      if(!any(is.na(cor_boot))){
        break;
      }
    })
    return( cor_boot[which(upper.tri(cor_boot))] )
  }))
})

wass_pw_w1 <- local({
  
  wass_pw <- matrix(0, ncol(pops[[1]]), ncol(pops[[1]]))
  wass_pw[upper.tri(wass_pw)] <- sapply(1:ncol(pop_boots[[1]]), function(ind){
    
    wass_pw <- wass_by_pair_w1(lapply(pop_boots, function(pop_boot) pop_boot[, ind]))
    
    return( wass_pw )
  })
  wass_pw <- wass_pw + t(wass_pw)
  
  return( wass_pw )
})

system.time({
  best_sub_pw_w1_1 <- local({
    best_sub <- apply(combn(1:ncol(pops[[1]]), 4), 2, function(ind){
      return( c(ind, wass_sum_subset_pw(wass_pw_w1, ind)) )
    })
    return( best_sub[, which.max(best_sub[5, ])] )
  })
})
best_sub_pw_w1_1

system.time({
  best_sub_pw_w1_2 <- sum.combn.mat(wass_pw_w1, 6, c(100, 100))
})
  best_sub_pw_w1_2
