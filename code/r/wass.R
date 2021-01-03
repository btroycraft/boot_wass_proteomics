wass_by_pair_w1 <- function(pairs){
  
  sorted <- lapply(pairs, sort)
  sorted <- do.call(cbind, sorted)
  corr_mid <- apply(sorted, 1, median)
  wass <- sweep(sorted, 1, corr_mid)
  wass <- apply(wass, 2, function(x) mean(abs(x)) )
  
  return( sum(wass) )
}

wass_by_pair_w2 <- function(pairs){
  
  sorted <- lapply(pairs, sort)
  sorted <- do.call(cbind, sorted)
  corr_mid <- apply(sorted, 1, mean)
  wass <- sweep(sorted, 1, corr_mid)
  wass <- apply(wass, 2, function(x) mean(x^2) )
  
  return( sum(wass) )
}

wass_sum_subset_joint <- function(samp_list, p=2){
  
  wass_high_dim <- function(X, Y, q=2){
    
    D <- apply(Y, 1, function(y){
      apply(X, 1, function(x){
        sqrt(sum((x-y)^2))^q
      })
    })
    
    match <- as.numeric(clue::solve_LSAP(D))
    
    return( mean(D[cbind(1:nrow(X), match)])^(1/p) )
  }
  
  comb_ind <- combn(length(samp_list), 2)
  
  wass_pairs <- apply(comb_ind, 2, function(ind) wass_high_dim(samp_list[[ind[1]]], samp_list[[ind[2]]], p) )
  
  return( sum(wass_pairs^p) )
}

wass_sum_subset_pw <- function(mat, ind) sum(mat[t(combn(ind, 2))])