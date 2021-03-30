library(ggplot2)

cor.dend <- function(cor_mat){
  
  cor_melt <-
    t(combn(nrow(cor_mat), 2)) %>%
      data.frame(cor_mat[lower.tri(cor_mat)]) %>%
      set_names(c('i', 'j', 'val')) %>%
      arrange(desc(val))
  
  group_list <-
    seq_len(nrow(cor_mat)) %>%
      lapply(. %>%
        list(birth = 1, death = nrow(cor_mat)+1, clust = .)
      )
  
  change <- c(1, rep(0, nrow(cor_mat)-1), -1)
  
  local({
    
    group_work <- seq_len(nrow(cor_mat))
    ind_change <- 2
    
    lapply(seq_len(nrow(cor_vec)-1), function(._){
      
      i <- cor_melt$i[._]
      j <- cor_melt$j[._]
      val <- cor_melt$val[._]
      
      if(group_work[i] != group_work[j]){
        group_pair <- sort(group_work[c(i, j)])
        
        group_list[[group_pair[1]]]$death <<- ind_change+1
        group_list[[group_pair[2]]]$death <<- ind_change+1
        group_list[[group_pair[1]]] <<-
          list(
            birth = ind_change+1,
            death = nrow(cor_mat)+2,
            clust = c(
              group_list[[group_pair[1]]]$clust,
              group_list[[group_pair[2]]]$clust),
            group_list[[group_pair[1]]],
            group_list[[group_pair[2]]]
          )
        group_list[[group_pair[2]]] <<- NA
        
        group_work[group_work == group_pair[2]] <<- group_pair[1]
        
        change[ind_change] <<- val
        ind_change <<- ind_change + 1
      }
    })
  })
  
  group_list <- group_list[!is.na(group_list)]
  
  ord <- order(group_list$clust)
  rnk <- rank(group_list$clust)
  
  .cor.dend <- function(node){
    node$clust <- rnk[node$clust]
    if(length(node) == 5){
      node[[4]] <- .cor.dend(node[[4]])
      node[[4]] <- .cor.dend(node[[5]])
    }
    
    return(node)
  }
  
  return(list(
    lab = ord,
    change = change,
    tree = .cor.dend(group_list[[1]])
  ))
}

dend.to.mat <- function(dend){
  
  n <- length(dend$change)-1
  
  mat <- matrix(0, n, n+1)
  
  .dend.to.mat <- function(tree){
    val <- mean(tree$clust)
    mat[tree$clust, tree$birth:(tree$death-1)] <<- val
    
    if(length(tree) == 5){
      mat <<- .dend.to.mat(tree$child1)
      mat <<- .dend.to.mat(tree$child2)
    }
  }
  
  .dend.to.mat(dend$tree)
  
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


