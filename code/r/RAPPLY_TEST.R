library(dplyr)
library(magrittr)
library(MASS)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(viridis)

source('r/load_data.R')



plot_dend <- function(._){
  
  indices <- IND_LIST[[._]]
  data_split <- DATA_GILL %>%
    dplyr::select(all_of(indices)) %>%
    split(DATA_GILL$location)
  
  folder_name <-
    c(list(sep = '_'), as.list(._)) %>%
    do.call(what = paste)
  
  plot_list <- lapply(names(data_split), function(name){
    data <- data_split[[name]]
    cor_mat <- cor(data)
    dend <- COR.DEND(cor_mat)
    GG.DEND(dend, labels = FALSE, xlim=c(0, 1), main = name)
  })
  
  plot_out <- do.call(grid.arrange, plot_list)
  
  plot_name <- sprintf('plots/plot_dend_%s.pdf', folder_name)
  ggsave(plot_name, plot_out, width = 7, height = 7)
  
  return(NULL)
}



plot_mds_string <- function(._){
  
  indices <- IND_LIST[[._]]
  data_split <- DATA_GILL %>%
    dplyr::select(all_of(indices)) %>%
    split(DATA_GILL$location)
  
  folder_name <-
    c(list(sep = '_'), as.list(._)) %>%
    do.call(what = paste)
  
  interactions <- 
    folder_name %>%
    sprintf(fmt = 'data/stringdb/%s/string_interactions.csv') %>%
    read.csv(stringsAsFactors = FALSE)
  clusters <- 
    folder_name %>%
    sprintf(fmt = 'data/stringdb/%s/string_nodes.csv') %>%
    read.csv(stringsAsFactors = FALSE)
  clusters$mcl_cluster <- as.factor(clusters$mcl_cluster)
  plot_name <- sprintf('plots/plot_mcl_mds_%s.pdf', folder_name)
  
  plot_list <- lapply(names(data_split), function(name){
    data <- data_split[[name]]
    cor_mat <- cor(data)
    GG.MDS(cor_mat, ._, name, clusters)
  })
  plots_mds <- do.call(grid.arrange, plot_list)
  plot_mcl_cluster <- GG.MCL.CLUSTER(clusters, ._, interactions)
  plot_out <- grid.arrange(plots_mds, plot_mcl_cluster)
  ggsave(plot_name, plot_out, width = 13, height = 13)
}

GG.MDS <- function(cor_mat, ._, env_name, clusters){
  matrix_ones = matrix(rep(1, nrow(cor_mat)^2), nrow(cor_mat), nrow(cor_mat))
  dist = matrix_ones - cor_mat
  fit = isoMDS(dist, k = 2)
  if("top" %in% ._){
    max_overlap = Inf
    alpha = 0.8
  } else{
    max_overlap = 20
    alpha = 0.5
  }
  mds_data <- data.frame(fit$points)
  mds_cluster_data <- merge(mds_data, clusters, by.x = 0, by.y = "node")
  ggplot(mds_cluster_data, aes(x = X1, y = X2)) +
    geom_point(aes(fill = mcl_cluster), pch = 21, size = 3, alpha = alpha) +
    theme(plot.title = element_text(size=20, hjust = 0.5, face = "bold"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "none") +
    labs(title = env_name)  +
    xlim(-2, 2) + ylim(-2,2)+ 
    scale_color_viridis(discrete = TRUE, option = "D") +
    scale_fill_viridis(discrete = TRUE)
  
}

GG.MCL.CLUSTER <- function(clusters, ._, interactions){
  temp <- merge(interactions, clusters,  by.x = "node1" , by.y = "node")
  if("top" %in% ._){
    max_overlap = Inf
    edge_data <- merge(temp, clusters, by.x = "node2", by.y = "node")
    alpha = 1
  } else{
    max_overlap = 20
    edge_data <- merge(temp, clusters, by.x = "node2", by.y = "node") %>%
      filter(combined_score > 0.9)
    alpha = 0.6
  }
  
  
  ggplot(clusters, aes(x = x_position, y = y_position))  +
    geom_segment(data=edge_data,aes(x=x_position.x,xend = x_position.y,
                                            y=y_position.x,yend = y_position.y),
                                            colour="royalblue2",
                                            size = edge_data$combined_score*1.7,
                                            alpha = alpha) +
    geom_point(aes(fill = mcl_cluster), size = 5, pch = 21) +
    geom_text_repel(label = clusters$node,
                    max.overlaps = max_overlap,
                    fontface = "bold",
                    size = 5.5) +
    theme(plot.title = element_text(size=20, hjust = 0.5, face = "bold"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position = "none") +
    labs(title = "MCL Clustering")  +
    theme(
      axis.text.x = element_blank(),  # remove x-axis text
      axis.text.y = element_blank(), # remove y-axis text
      axis.ticks = element_blank(),  # remove axis ticks
      axis.title.x = element_blank(), # remove x-axis labels
      axis.title.y = element_blank()) + 
    scale_color_viridis(discrete = TRUE, option = "D") +
    scale_fill_viridis(discrete = TRUE)
  
  
}



rapply(REF_LIST[c('max', 'top')], how = 'replace', plot_mds_string)



