data_gill <- read.csv('cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv')
indices <- read.csv('../julia/selected_subsets.csv')[, 1]

data_list <- split(data_gill[, !(colnames(data_gill) %in% c('sample', 'location'))], data_gill$location)
data_list <- lapply(data_list, function(data){
  data[, indices]
})

get.ripser <- function(X, t, q){
  
  diag <- as.data.frame(ripserr::vietoris_rips(as.dist(1-cor(X)), threshold = t, max_dim = q))
  colnames(diag) <- c('dim', 'birth', 'death')
  
  return(diag)
}

sub.diag <- function(diag, q){
  return( diag[diag$dim == q, c('birth', 'death')] )
}

pbetti <- function(diag, r, s, q) sum(diag$dim == q & diag$birth <= r & diag$death > s)

x <- seq(0, 1, length.out = 10^3)

b0_list <- lapply(data_list, function(data){
  diag <- get.ripser(data, 1, 1)
  sapply(x, function(x){
    pbetti(diag, 1-x, 1-x, 0)
  })
})

bound_b0 <- lapply(data_list, function(data){
  boot_reps <- replicate(200, {
    X_boot <- data[sample(nrow(data), replace = TRUE), ]
    diag <- get.ripser(X_boot, 1, 1)
    sapply(x, function(x){
      pbetti(diag, 1-x, 1-x, 0)
    })
  })
  
  apply(boot_reps, 1, function(row){
    quantile(abs(row - mean(row)), .95)
  })
})

pdf('betti0.pdf', 7, 7)
  for(ind1 in 1:3){
    for(ind2 in (ind1+1):4){
      plot(c(), xlim=c(0, 1), ylim=c(0, 100), xlab = 'Correlation', ylab = 'Betti 0',
           main = paste(names(data_list)[ind1], '
                        vs.
                        ', names(data_list)[ind2]))
      
      y_l <- b0_list[[ind1]] - bound_b0[[ind1]]
      y_u <- b0_list[[ind1]] + bound_b0[[ind1]]
      polygon(c(x,rev(x)), c(y_u, rev(y_l)), col=ind1)
      
      y_l <- b0_list[[ind2]] - bound_b0[[ind2]]
      y_u <- b0_list[[ind2]] + bound_b0[[ind2]]
      polygon(c(x,rev(x)), c(y_u, rev(y_l)), col=ind2)
    }
  }
dev.off()

b1_list <- lapply(data_list, function(data){
  diag <- get.ripser(data, 1, 1)
  sapply(x, function(x){
    pbetti(diag, 1-x, 1-x, 1)
  })
})

bound_b1 <- lapply(data_list, function(data){
  boot_reps <- replicate(200, {
    X_boot <- data[sample(nrow(data), replace = TRUE), ]
    diag <- get.ripser(X_boot, 1, 1)
    sapply(x, function(x){
      pbetti(diag, 1-x, 1-x, 1)
    })
  })
  
  apply(boot_reps, 1, function(row){
    quantile(abs(row - mean(row)), .95)
  })
})

pdf('betti1.pdf', 7, 7)
for(ind1 in 1:3){
  for(ind2 in (ind1+1):4){
    plot(c(), xlim=c(0, 1), ylim=c(0, 30), xlab = 'Correlation', ylab = 'Betti 1',
         main = paste(names(data_list)[ind1], '
                        vs.
                        ', names(data_list)[ind2]))
    
    y_l <- b1_list[[ind1]] - bound_b1[[ind1]]
    y_u <- b1_list[[ind1]] + bound_b1[[ind1]]
    polygon(c(x,rev(x)), c(y_u, rev(y_l)), col=ind1)
    
    y_l <- b1_list[[ind2]] - bound_b1[[ind2]]
    y_u <- b1_list[[ind2]] + bound_b1[[ind2]]
    polygon(c(x,rev(x)), c(y_u, rev(y_l)), col=ind2)
  }
}
dev.off()
