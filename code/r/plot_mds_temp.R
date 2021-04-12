library(gridExtra)
library(cowplot)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggrepel)

source('r/funcs.R')
source('r/load_data.R')

col_list <-
  seq(0, 1, length.out = 49) %>%
  {c(
    .[seq.int(from = 1, by = 16, length.out = 3)],
    .[seq.int(from = 9, by = 16, length.out = 3)],
    .[seq.int(from = 5, by = 8, length.out = 6)],
    .[seq.int(from = 3, by = 4, length.out = 12)],
    .[seq.int(from = 2, by = 2, length.out = 24)]
  )} %>%
  hsv(s = 1, v = 1)

rapply(REF_LIST$max, how = 'replace', function(._){
  
  data_split <-
    DATA_GILL %>%
    select(IND_LIST[[._]]) %>%
    split(DATA_GILL$location)
  
  string_nodes <-
    ._ %>%
    as.list %>%
    c(list(sep = '_' )) %>%
    do.call(what = paste) %>%
    sprintf(fmt = 'data/stringdb/%s/string_nodes.csv') %>%
    read.csv(stringsAsFactors = FALSE)
  
  col_list_node <- 
    sapply(IND_LIST[[._]], function(._){
      clust <- string_nodes$mcl_cluster[string_nodes$node == ._]
      if(length(clust) > 0){
        if(sum(string_nodes$mcl_cluster == clust) > 1){
          clust
        } else {
          0
        }
      } else {
        0
      }
    })
  
  string_interactions <-
    ._ %>%
    as.list %>%
    c(list(sep = '_' )) %>%
    do.call(what = paste) %>%
    sprintf(fmt = 'data/stringdb/%s/string_interactions.csv') %>%
    read.csv(stringsAsFactors = FALSE)
  
  mds_out <-
    lapply(names(data_split), function(._){
      
      data_split[[._]] %>%
      {1 - cor(.)} %>%
      cmdscale %>%
      as.data.frame %>%
      mutate(col = col_list_node) %>%
      ggplot() +
        geom_point(
          aes(
            x = V1 - (max(V1) + min(V1)) / 2,
            y = V2 - (max(V2) + min(V2)) / 2,
          ),
          show.legend = FALSE,
          col = 'black',
          size = 1.2
        ) +
        geom_point(
          aes(
            x = V1 - (max(V1) + min(V1)) / 2,
            y = V2 - (max(V2) + min(V2)) / 2,
            col = as.character(col)
          ),
          show.legend = FALSE,
          size = .7
        ) +
        coord_fixed(
          ratio = 1,
          xlim = c(-1, 1),
          ylim = c(-1, 1)
        ) +
        labs(title = ._) +
        scale_color_manual(values = c('gray80', unique(col_list))) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = .5),
          axis.title = element_blank(),
          axis.text =  element_blank()
        )
    }) %>%
    {plot_grid(plotlist = ., align = 'hv', nrow = 2)}
    
    edge_list <-
      string_interactions %>%
      merge(string_nodes, by.x = 'node1', by.y = 'node') %>%
      rename(x_pos1 = x_position, y_pos1 = y_position) %>%
      select(
        c('node2', 'combined_score', 'x_pos1', 'y_pos1')
      ) %>%
      merge(string_nodes, by.x = 'node2', by.y = 'node') %>%
      rename(x_pos2 = x_position, y_pos2 = y_position) %>%
      select(
        c('combined_score', 'x_pos1', 'y_pos1', 'x_pos2', 'y_pos2')
      ) %>%
      filter(combined_score >= .5)
  
    string_out <-
      string_nodes %>%
      mutate(col = col_list_node[node]) %>%
      ggplot() +
        geom_segment(
          aes(
            x = x_pos1,
            y = y_pos1,
            xend = x_pos2,
            yend = y_pos2
          ),
          data = edge_list
        ) +
        geom_point(
          aes(
            x = x_position,
            y = y_position
          ),
          col = 'black',
          size = 3,
          show.legend = FALSE
        ) +
        geom_point(
          aes(
            x = x_position,
            y = y_position,
            col = as.character(col)
          ),
          size = 2,
          show.legend = FALSE
        ) +
        coord_fixed(
          ratio = 1,
          xlim = c(0, 1),
          ylim = c(0, 1)
        ) +
        scale_color_manual(values = c('gray80', unique(col_list))) +
        geom_label_repel(
          aes(
            x = x_position,
            y = y_position,
            label = node
          ),
          size = 2,
          label.padding = unit(1.5, 'pt')
        ) +
        theme_bw() +
        theme(
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text =  element_blank()
        )
    
    plot_out <- grid.arrange(mds_out, string_out, nrow = 1)
    
  ._ %>%
    as.list %>%
    c(list(sep = '_')) %>%
    do.call(what = paste) %>%
    sprintf(fmt = 'plots/stringmds_%s.pdf') %>%
    ggsave(plot = plot_out, width = 12.5, height = 6.5)
})

rapply(REF_LIST$top, how = 'replace', function(._){
  
  data_split <-
    DATA_GILL %>%
    select(IND_LIST[[._]]) %>%
    split(DATA_GILL$location)
  
  string_nodes <-
    ._ %>%
    as.list %>%
    c(list(sep = '_' )) %>%
    do.call(what = paste) %>%
    sprintf(fmt = 'data/stringdb/%s/string_nodes.csv') %>%
    read.csv(stringsAsFactors = FALSE)
  
  col_list_node <- 
    sapply(IND_LIST[[._]], function(._){
      clust <- string_nodes$mcl_cluster[string_nodes$node == ._]
      if(length(clust) > 0){
        if(sum(string_nodes$mcl_cluster == clust) > 1){
          clust
        } else {
          0
        }
      } else {
        0
      }
    })
  
  string_interactions <-
    ._ %>%
    as.list %>%
    c(list(sep = '_' )) %>%
    do.call(what = paste) %>%
    sprintf(fmt = 'data/stringdb/%s/string_interactions.csv') %>%
    read.csv(stringsAsFactors = FALSE)
  
  mds_out <-
    lapply(names(data_split), function(._){
      
      data_split[[._]] %>%
        {1 - cor(.)} %>%
        cmdscale %>%
        as.data.frame %>%
        mutate(col = col_list_node, node = row.names(.)) %>%
        ggplot() +
        geom_point(
          aes(
            x = V1 - (max(V1) + min(V1)) / 2,
            y = V2 - (max(V2) + min(V2)) / 2,
          ),
          show.legend = FALSE,
          col = 'black',
          size = 1.2
        ) +
        geom_point(
          aes(
            x = V1 - (max(V1) + min(V1)) / 2,
            y = V2 - (max(V2) + min(V2)) / 2,
            col = as.character(col)
          ),
          show.legend = FALSE,
          size = .7
        ) +
        geom_label_repel(
          aes(
            x = V1 - (max(V1) + min(V1)) / 2,
            y = V2 - (max(V2) + min(V2)) / 2,
            label = node
          ),
          size = 2,
          label.padding = unit(1.5, 'pt')
        ) +
        coord_fixed(
          ratio = 1,
          xlim = c(-1, 1),
          ylim = c(-1, 1)
        ) +
        labs(title = ._) +
        scale_color_manual(values = c('gray80', unique(col_list))) +
        theme_bw() +
        theme(
          plot.title = element_text(hjust = .5),
          axis.title = element_blank(),
          axis.text =  element_blank()
        )
    }) %>%
    {plot_grid(plotlist = ., align = 'hv', nrow = 2)}
  
  edge_list <-
    string_interactions %>%
    merge(string_nodes, by.x = 'node1', by.y = 'node') %>%
    rename(x_pos1 = x_position, y_pos1 = y_position) %>%
    select(
      c('node2', 'combined_score', 'x_pos1', 'y_pos1')
    ) %>%
    merge(string_nodes, by.x = 'node2', by.y = 'node') %>%
    rename(x_pos2 = x_position, y_pos2 = y_position) %>%
    select(
      c('combined_score', 'x_pos1', 'y_pos1', 'x_pos2', 'y_pos2')
    ) %>%
    filter(combined_score >= .5)
  
  string_out <-
    string_nodes %>%
    mutate(col = col_list_node[node]) %>%
    ggplot() +
    geom_segment(
      aes(
        x = x_pos1,
        y = y_pos1,
        xend = x_pos2,
        yend = y_pos2
      ),
      data = edge_list
    ) +
    geom_point(
      aes(
        x = x_position,
        y = y_position
      ),
      col = 'black',
      size = 3,
      show.legend = FALSE
    ) +
    geom_point(
      aes(
        x = x_position,
        y = y_position,
        col = as.character(col)
      ),
      size = 2,
      show.legend = FALSE
    ) +
    coord_fixed(,
      xlim = c(0, 1),
      ylim = c(0, 1),
      ratio = 1
    ) +
    scale_color_manual(values = c('gray80', unique(col_list))) +
    geom_label_repel(
      aes(
        x = x_position,
        y = y_position,
        label = node
      ),
      size = 2,
      label.padding = unit(1.5, 'pt')
    ) +
    theme_bw() +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text =  element_blank()
    )
  
  plot_out <- grid.arrange(mds_out, string_out, nrow = 1)
  
  ._ %>%
    as.list %>%
    c(list(sep = '_')) %>%
    do.call(what = paste) %>%
    sprintf(fmt = 'plots/stringmds_%s.pdf') %>%
    ggsave(plot = plot_out, width = 12.5, height = 6.5)
})
