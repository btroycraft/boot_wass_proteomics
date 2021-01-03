library(ggplot2)
library(igraph)
library(ggrepel)
library(MASS)

data_gill <- read.csv('cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv')
prot_names <- names(data_gill)[!(names(data_gill) %in% c('sample', 'location'))]

data_list <- split(data_gill[!(names(data_gill) %in% c('sample', 'location'))], data_gill$location)
cor_list <- lapply(data_list, cor)
dist_list <- lapply(cor_list, function(x) 1-x)

names_list <- {
  x <- c('total', 'bodega', 'rosaria', 'solano', 'westchester')
  names(x) <- c('Total', names(data_list))
  x
}

for(name in names_list){
  
  indices <- {
    x <- read.csv(sprintf('../data/ind_top_%s_betti0.csv', name))
    x <- x[names(x) != 'val']
    as.character(x[1, ])
  }
  mds_list <- lapply(dist_list, function(dist){
    x <- isoMDS(as.dist(dist[indices, indices]),
           maxit = 10^3,
           tol = 10^-10,
           trace = FALSE)
    x$points
  })
  plot_frame <- do.call(rbind, lapply(1:length(mds_list), function(ind){
    cbind(as.data.frame(mds_list[[ind]]),
          group = names(data_list)[ind],
          label = row.names(mds_list[[ind]]))
  }))
  plot <- ggplot(data = plot_frame, aes(x = V1, y = V2, label = label)) +
    geom_point() +
    geom_label_repel(size = 3,
                     position = 'nudge',
                     label.size = 0,
                     ) +
    facet_wrap(vars(group)) +
    xlab('') +
    ylab('')
  ggsave(sprintf('../plots/mds_%s_betti0.pdf', name), plot, height = 7, width = 7) 
}

for(name in names_list){
  
  indices <- {
    x <- read.csv(sprintf('../data/ind_top_%s_betti1.csv', name))
    x <- x[names(x) != 'val']
    as.character(x[1, ])
  }
  mds_list <- lapply(dist_list, function(dist){
    x <- isoMDS(as.dist(dist[indices, indices]),
                maxit = 10^3,
                tol = 10^-10,
                trace = FALSE)
    x$points
  })
  plot_frame <- do.call(rbind, lapply(1:length(mds_list), function(ind){
    cbind(as.data.frame(mds_list[[ind]]),
          group = names(data_list)[ind],
          label = row.names(mds_list[[ind]]))
  }))
  plot <- ggplot(data = plot_frame, aes(x = V1, y = V2, label = label)) +
    geom_point() +
    geom_label_repel(size = 3,
                     position = 'nudge',
                     label.size = 0,
    ) +
    facet_wrap(vars(group)) +
    xlab('') +
    ylab('')
  ggsave(sprintf('../plots/mds_%s_betti1.pdf', name), plot, height = 7, width = 7) 
}

for(name in names_list){
  
  indices <- {
    x <- read.csv(sprintf('../data/ind_max_%s.csv', name))
    x <- x[names(x) != 'val']
    as.character(x[1, ])
  }
  mds_list <- lapply(dist_list, function(dist){
    x <- isoMDS(as.dist(dist[indices, indices]),
                maxit = 10^3,
                tol = 10^-10,
                trace = FALSE)
    x$points
  })
  plot_frame <- do.call(rbind, lapply(1:length(mds_list), function(ind){
    cbind(as.data.frame(mds_list[[ind]]),
          group = names(data_list)[ind],
          label = row.names(mds_list[[ind]]))
  }))
  plot <- ggplot(data = plot_frame, aes(x = V1, y = V2, label = label)) +
    geom_point() +
    geom_label_repel(size = 3,
                     position = 'nudge',
                     label.size = 0,
    ) +
    facet_wrap(vars(group)) +
    xlab('') +
    ylab('')
  ggsave(sprintf('../plots/mds_%s.pdf', name), plot, height = 7, width = 7) 
}
