source("selecting_proteins.R")
source('wass.R')
source('gen_boot.R')
source('sum_combn_mat.R')
library(readxl)
library(Matrix)
library(parallel)


get.labs <- function(ind, n){
  triag <- (0:(n-1))*(1:n)/2
  ind2 <- sum(triag < k)+1
  ind1 <- k - triag[ind2-1]
  
  return(c(ind1, ind2))
}


vec_cor <- function(df){
  temp = cor(df[1:5])
  temp[which(lower.tri(temp))]
}

compute_var_sum <- function(permute_data) {
  env_corrs = by(permute_data, permute_data$label, FUN = vec_cor)
  env_corrs = as.data.frame(do.call("rbind",env_corrs))
  sum_var = sum(unlist(lapply(env_corrs, var)))
  return(sum_var)
}


boot_reps <- 24

data <- read.table('cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv', header=TRUE, sep=',')
#data <- as.data.frame(readxl::read_excel("./cDIA_MXLSAKBH-Exp1-2-3-4_Gill_with Volcanoplots_SelectedSubset.xlsx"))
data <- cbind(data[c(sig_prots_special)],location = data[,2])

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

best_sub_pw_w1_2 <- sum.combn.mat(wass_pw_w1, 5, c(50, 50))



prot_data = labeled_data[, 1:1505]
permute_data = data.frame(lapply(prot_data, function(x) as.numeric(as.character(x))))
permute_data = cbind(permute_data, label = labeled_data$label)
observed_data = cbind(data.frame(lapply(prot_data, function(x) as.numeric(as.character(x)))), label = labeled_data$label)



perms = data.frame(replicate(1000,sample(observed_data$label)))


perm_test_var_par <- function(permute_data, best_sub, perms, num_cores = 1) {
  stats <- mclapply(data.frame(best_sub), function(ind){
    sub_prots = sig_prots_special[ind]
    summed_variances = mclapply(perms, function(perm){
      permute_data["label"] = perm
      permute_data2 = cbind(permute_data[,sub_prots],label = permute_data$label)
      sum_var = compute_var_sum(permute_data2)
      return(sum_var)
    }, mc.cores = num_cores)
    observed_data_sub = cbind(observed_data[,sub_prots], label = observed_data$label)
    test_stat = compute_var_sum(observed_data_sub)
    pval = sum(test_stat < summed_variances) / ncol(perms)
    return(c(test_stat,pval,paste(sub_prots, collapse = "-"), paste(ind, collapse = " ")))
  }, mc.cores = num_cores)
  stats <- data.frame(do.call("rbind", stats))
  colnames(stats) <- c("test.stat","p.val","sub", "sub_indices")
  return(stats)
}

stats_min = perm_test_var_par(permute_data, best_sub_pw_w1_2$min$indices, perms)
stats_max = perm_test_var_par(permute_data, best_sub_pw_w1_2$max$indices, perms)

write.csv(stats_min,file = "min_perms.csv")
write.csv(stats_max,file = "max_perms.csv")
