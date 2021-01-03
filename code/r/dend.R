data_gill <- read.csv('cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv')
data_list <- split(data_gill[, !(colnames(data_gill) %in% c('sample', 'location'))], data_gill$location)
cor_list <- lapply(data_list, cor)

indices <- read.csv('../julia/selected_subsets.csv')

get.dend <- function(cor_mat){
  
  cor_vec <- cbind(cor_mat[lower.tri(cor_mat)], t(combn(1:nrow(cor_mat), 2)))
  cor_vec <- cor_vec[order(cor_vec[, 1], decreasing = TRUE), ]
  
  group <- 1:nrow(cor_mat)
  group_list <- lapply(group, function(x) list(birth = 1, death = nrow(cor_mat)+1, clust = x))
  change <- numeric(nrow(cor_mat)+1)
  change[1] <- 1
  change[length(change)] <- -1
  
  ind_change <- 2
  for(ind in 1:nrow(cor_vec)){
    i <- cor_vec[ind, 2]
    j <- cor_vec[ind, 3]
    val <- cor_vec[ind, 1]
    if(group[i] != group[j]){
      group_pair <- sort(group[c(i, j)])
      
      group_list[[group_pair[1]]]$death <- ind_change+1
      group_list[[group_pair[2]]]$death <- ind_change+1
      group_list[[group_pair[1]]] <- list(birth = ind_change+1,
                                         death = nrow(cor_mat)+2,
                                         clust = c(group_list[[group_pair[1]]]$clust, group_list[[group_pair[2]]]$clust),
                                         group_list[[group_pair[1]]],
                                         group_list[[group_pair[2]]])
      group_list[[group_pair[2]]] <- NA
      
      group[group == group_pair[2]] <- group_pair[1]
      
      change[ind_change] <- val
      ind_change <- ind_change + 1
    }
  }
  
  group_list <- group_list[!is.na(group_list)]
  ord <- order(group_list[[1]]$clust)
  
  get.dend_ <- function(node, ord){
    node$clust <- ord[node$clust]
    if(length(node) == 5){
      node[[4]] <- get.dend_(node[[4]], ord)
      node[[5]] <- get.dend_(node[[5]], ord)
    }
    
    return(node)
  }
  
  group_list <- list(index = group_list[[1]]$clust,
                     change = change,
                     tree = get.dend_(group_list[[1]], ord))
  
  
  return( group_list )
}

dend.to.mat <- function(dend){
  
  n <- length(dend$change)-1
  
  mat <- matrix(0, n, n+1)
  
  dend.to.mat_ <- function(tree){
    val <- mean(tree$clust)
    mat[tree$clust, tree$birth:(tree$death-1)] <<- val
    
    if(length(tree) == 5){
      mat <- dend.to.mat_(tree[[4]])
      mat <- dend.to.mat_(tree[[5]])
    }
  }
  
  dend.to.mat_(dend$tree)
  
  return(mat[, -1])
}

dend.mat.x <- function(x, dend){
  
  mat <- dend.to.mat(dend)
  mat <- mat[, ncol(mat):1]
  
  change <- dend$change
  change <- change[length(change):1]
  
  sapply(x, function(x){
    
    ind_l <- max(which(change <= x))
    
    return( mat[, ind_l] )
  })
}

dend.mat.grad.x <- function(x, dend){
  
  mat <- dend.to.mat(dend)
  mat <- mat[, ncol(mat):1]
  
  change <- dend$change
  change <- change[length(change):1]
  
  sapply(x, function(x){
    
    ind_l <- max(which(change <= x))
    
    if(ind_l == 1){
      return( rep((ncol(mat)+1)/2, nrow(mat)) )
    }
    
    x_l <- change[ind_l]
    x_u <- change[ind_l+1]
    
    a <- (x-x_l)/(x_u-x_l)
    
    val_l <- mat[, ind_l-1]
    val_u <- mat[, ind_l]
    
    return( ((1-a)*val_l+a*val_u) )
  })
}

plot.dend <- function(dend){
  
  plot
  
  change <- dend$change
  
  plot.dend.rec <- function(dend, col){
    
  }
}

pdf('dend.pdf', 7, 7)
  for(ind in 1:length(cor_list)){
    
    cor_mat = cor_list[[ind]]
    
    dend <- get.dend(cor_mat[indices[, 1], indices[, 1]])
    x <- seq(.2, 1, length.out=10^3); x <- x[c(-1, -length(x))]
    dend_mat <- dend.mat.x(x, dend)
    
    image(x = x, y = 1:nrow(dend_mat), z = t(dend_mat),
          col=hcl.colors(10^3, "YlOrRd", rev = TRUE),
          xlab='cor',
          ylab='index',
          main=names(cor_list)[ind])
  }
dev.off()


