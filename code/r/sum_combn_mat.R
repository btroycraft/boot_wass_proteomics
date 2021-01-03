sum.combn.mat <- function(X, k=3, num_out=c(1L, 1L)){
  
  k <- as.integer(k)
  num_out <- as.integer(num_out)
  
  dyn.load('sum_combn_mat')
  
  list_out <- .Call('sum_combn_mat', X[upper.tri(X)], nrow(X), k, num_out, package='sum_combn_mat')
  
  names(list_out) <- c('min', 'max')
  names(list_out$min) <- c('indices', 'sum')
  names(list_out$max) <- c('indices', 'sum')
  
  dim(list_out$min$indices) <- c(k, num_out[1])
  dim(list_out$max$indices) <- c(k, num_out[2])
  
  dyn.unload('sum_combn_mat')
  
  return(list_out)
}
