RIPSER <- function(D, q = 1L, t = 2){
  diag <- as.data.frame(ripserr::vietoris_rips(as.dist(D), max_dim = q, threshold = t))
  names(diag) <- c('Dim', 'Birth', 'Death')
  return( diag )
}

SUB.DIAG <- function(diag, q = 1L){
  diag <- subset(diag, Dim == q)
  row.names(diag) <- seq_len(nrow(diag))
  return( diag )
}