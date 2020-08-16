sum.combn.mat <- function(X, k){
  dyn.load('sum_combn_mat.dll')
  list_out <- .Call('sum_combn_mat', X, as.integer(k))
  names(list_out) <- c('indices', 'maxsum')
  dyn.unload('sum_combn_mat.dll')
  return(list_out)
}
