library(ggplot2)
library(igraph)
library(ggrepel)
library(MASS)
library(reshape2)

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
  plot_frame <- do.call(rbind, lapply(1:length(dist_list), function(ind){
    x <- dist_list[[ind]][indices, indices]
    x <- melt(x)
    x$group <- names(dist_list)[ind]; x
  }))
  plot <- ggplot(data = plot_frame, aes(x = Var1, y = Var2, fill = value)) +
    geom_raster() +
    facet_wrap(vars(group)) +
    xlab('') +
    ylab('')
  ggsave(sprintf('../plots/dist_%s_betti0.pdf', name), plot, height = 7, width = 7) 
}

for(name in names_list){
  
  indices <- {
    x <- read.csv(sprintf('../data/ind_top_%s_betti1.csv', name))
    x <- x[names(x) != 'val']
    as.character(x[1, ])
  }
  plot_frame <- do.call(rbind, lapply(1:length(dist_list), function(ind){
    x <- dist_list[[ind]][indices, indices]
    x <- melt(x)
    x$group <- names(dist_list)[ind]; x
  }))
  plot <- ggplot(data = plot_frame, aes(x = Var1, y = Var2, fill = value)) +
    geom_raster() +
    facet_wrap(vars(group)) +
    xlab('') +
    ylab('')
  ggsave(sprintf('../plots/dist_%s_betti1.pdf', name), plot, height = 7, width = 7) 
}
