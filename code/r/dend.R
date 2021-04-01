COR.DEND <- function(cor_mat){
  
  cor_melt <- local({
    combn <- rbind(
      combn(nrow(cor_mat), 2),
      cor_mat[lower.tri(cor_mat)]
    )
    combn <- combn[, order(combn[3, ], decreasing = TRUE)]
    lapply(seq_len(ncol(combn)), function(._){
      setNames(combn[, ._], c('i', 'j', 'val'))
    })
  })
  
  group_list <-
    local({
      group_list <-
        lapply(seq_len(nrow(cor_mat)), function(._){
          list(birth = 1, death = -1, clust = ._)
        })
      group_work <- seq_len(nrow(cor_mat))
      for(._ in cor_melt){
        if(group_work[._['i']] != group_work[._['j']]){
          group_pair <- sort(group_work[._[c('i', 'j')]])
          group_list[[group_pair[1]]]$death <-
            group_list[[group_pair[2]]]$death <-
              unname(._['val'])
          group_list[[group_pair[1]]] <-
            list(
              group_list[[group_pair[1]]],
              group_list[[group_pair[2]]],
              birth = unname(._['val']),
              death = -1,
              clust = c(
                group_list[[group_pair[1]]]$clust,
                group_list[[group_pair[2]]]$clust
              )
            )
          group_list[[group_pair[2]]] <- NA
          group_work[group_work == group_pair[2]] <- group_pair[1]
        }
      }
      group_list[[1]]
    })
  ord <- order(group_list$clust)
  reorder <- function(._){
    ._$clust <- ord[._$clust]
    if(length(._) == 5){
      ._[[1]] <- reorder(._[[1]])
      ._[[2]] <- reorder(._[[2]])
    }
    return(._)
  }
  return(list(
    lab = if(is.null(rownames(cor_mat))){
      group_list$clust
    } else {
      rownames(cor_mat)[group_list$clust]
    },
    tree = reorder(group_list)
  ))
}

GG.DEND <- function(dend, labels = TRUE, main = NULL, col = c('cornsilk', 'firebrick4'), xlim = c(-1, 1)){
  
  require(ggplot2)
  
  dend.unwrap <- function(._){
    out <- list(c(
      xmin = ._$death,
      xmax = ._$birth,
      ymin = min(._$clust)-.5,
      ymax = max(._$clust)+.5,
      val = mean(._$clust)
    ))
    if(length(._) == 5){
      out <- c(out, dend.unwrap(._[[1]]), dend.unwrap(._[[2]]))
    }
    return(out)
  }
  
  plot_frame <- as.data.frame(do.call(rbind, dend.unwrap(dend$tree)))
  plot_out <-
    ggplot(data = plot_frame) +
      coord_cartesian(xlim = xlim) +
      geom_rect(
        aes(
          xmin = xmin,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax,
          fill = val
        ),
        show.legend = FALSE
      ) +
      theme_classic() +
      xlab('Correlation') +
      scale_fill_gradient(low = col[1], high = col[2]) +
      theme(
        axis.line.y = element_blank(),
        plot.title = element_text(hjust = .5)
      )
    
  if(labels == TRUE){
    plot_out <- plot_out +
      scale_y_continuous(
        breaks = dend$tree$clust,
        labels = dend$lab,
        position = 'right'
      )
  } else {
    plot_out <- plot_out +
      scale_y_continuous(
        breaks = dend$tree$clust,
        position = 'right'
      ) +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }
  
  if(!is.null(main)){
    plot_out <- plot_out + labs(title = main)
  }
  
  plot_out
}
