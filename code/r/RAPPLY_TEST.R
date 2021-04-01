library(dplyr)
library(magrittr)

source('r/load_data.R')
source('r/dend.R')

plot_mds_string <- function(._){
  
  indices <- IND_LIST[[._]]
  data <- DATA_GILL %>%
    select(all_of(ind))
  
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
  
  plot_name <- sprintf('plots/plot_dend_%s.pdf', folder_name)
}

plot_dend <- function(._){
  
  indices <- IND_LIST[[._]]
  data_split <- DATA_GILL %>%
    select(all_of(indices)) %>%
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

rapply(REF_LIST[c('max', 'top')], how = 'replace', plot_dend)