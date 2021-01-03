temp <- function(X, reps, indices){
  dyn.load('new')
  sum <- .Call('get_cor', X, as.integer(reps), as.integer(indices))
  dyn.unload('new')
  return(sum)
}

setwd('~/GitHub/boot_wass_proteomics/r')

system('Rcmd SHLIB new.c')

n <- 100
d <- 5
X <- matrix(rnorm(n*d), n, d)
indices <- sample(nrow(X), nrow(X), replace=TRUE)
X_boot <- X[indices, ]
reps <- table(indices)
indices <- as.integer(names(reps))
reps <- as.integer(reps)
names(reps) <- NULL

cor(X_boot)
