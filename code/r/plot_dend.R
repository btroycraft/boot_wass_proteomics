library(gridExtra)
library(cowplot)
library(dplyr)
library(magrittr)
library(ggplot2)

source('r/funcs.R')
source('r/load_data.R')

rapply(REF_LIST$max, how = 'replace', function(._){
   
  data_split <-
    DATA_GILL %>%
    select(IND_LIST[[._]]) %>%
    split(DATA_GILL$location)
  
  pbetti_frame <-
    lapply(names(data_split), function(._){
      data_split[[._]] %>%
        {1-cor(.)} %>%
        RIPSER(q = 0, t = 2) %>%
        {sapply(seq(0, 2, length.out = 10^3), function(._){
          PBETTI(., r = ._, q = 0)
        })} %>%
        {data.frame(
          cor = seq(1, -1, length.out = 10^3),
          pbetti = pmax(., 1),
          location = ._
        )}
    }) %>%
    do.call(what = rbind)
  
  plot_out <-
    lapply(names(data_split), function(._){
      
      plot1 <- 
        data_split[[._]] %>%
        cor %>%
        COR.DEND %>%
        GG.DEND(labels = FALSE, outline = TRUE, xlim = c(0, 1))
      
      col_list <-
        rep('gray70', length(data_split)) %>%
        set_names(names(data_split))
      col_list[._] <- 'black'
      
      plot2 <- 
        pbetti_frame %>%
        ggplot() +
        labs(title = ._) +
        geom_line(
          aes(
            x = cor,
            y = pbetti,
            group = location
          ),
          show.legend = FALSE,
          color = 'gray70'
        ) +
        geom_line(
          aes(
            x = cor,
            y = pbetti,
            group = location
          ),
          data =
            pbetti_frame %>%
            filter(location == ._),
          show.legend = FALSE,
          color = 'black',
          size = 1
        ) +
        theme_bw() +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = .5),
          panel.grid.major.x = element_line(color = 'gray91', size = .4),
          panel.grid.major.y = element_line(color = 'gray91', size = .4),
          panel.grid.minor.x = element_line(color = 'gray91', size = .2),
          panel.grid.minor.y = element_line(color = 'gray91', size = .2)
        ) +
        ylab('Betti-0') +
        coord_cartesian(xlim = c(0, 1))
      
      align_plots(plot2, plot1, align = 'v') %>%
        {arrangeGrob(grobs = ., nrow = 2, heights = c(1, 1.5))}
    }) %>%
    arrangeGrob(grobs = ., nrow = 1)
  
  ._ %>%
    as.list %>%
    c(list(sep = '_')) %>%
    do.call(what = paste) %>%
    sprintf(fmt = 'plots/dend_%s.pdf') %>%
    ggsave(plot = plot_out, width = 6.5, scale = 2)
})


rapply(REF_LIST$top, how = 'replace', function(._){
  
  data_split <-
    DATA_GILL %>%
    select(IND_LIST[[._]]) %>%
    split(DATA_GILL$location)
  
  pbetti_frame <-
    lapply(names(data_split), function(._){
      data_split[[._]] %>%
        {1-cor(.)} %>%
        RIPSER(q = 0, t = 2) %>%
        {sapply(seq(0, 2, length.out = 10^3), function(._){
          PBETTI(., r = ._, q = 0)
        })} %>%
        data.frame(
          cor = seq(1, -1, length.out = 10^3),
          pbetti = pmax(., 1),
          location = ._
        )
    }) %>%
    do.call(what = rbind)
  
  plot_out <-
    lapply(names(data_split), function(._){
      
      plot1 <- 
        data_split[[._]] %>%
        cor %>%
        COR.DEND %>%
        GG.DEND(labels = TRUE, outline = TRUE, xlim = c(-.6, 1))
      
      col_list <-
        rep('gray70', length(data_split)) %>%
        set_names(names(data_split))
      col_list[._] <- 'black'
      
      plot2 <- 
        pbetti_frame %>%
        ggplot() +
        labs(title = ._) +
        geom_line(
          aes(
            x = cor,
            y = pbetti,
            group = location
          ),
          show.legend = FALSE,
          color = 'gray70'
        ) +
        geom_line(
          aes(
            x = cor,
            y = pbetti,
            group = location
          ),
          data =
            pbetti_frame %>%
            filter(location == ._),
          show.legend = FALSE,
          color = 'black',
          size = 1
        ) +
        theme_bw() +
        theme(
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = .5),
          panel.grid.major.x = element_line(color = 'gray91', size = .4),
          panel.grid.major.y = element_line(color = 'gray91', size = .4),
          panel.grid.minor.x = element_line(color = 'gray91', size = .2),
          panel.grid.minor.y = element_line(color = 'gray91', size = .2)
        ) +
        ylab('Betti-0') +
        coord_cartesian(xlim = c(-.6, 1))
      
      align_plots(plot2, plot1, align = 'v') %>%
        {arrangeGrob(grobs = ., nrow = 2, heights = c(1, 1.5))}
    }) %>%
    arrangeGrob(grobs = ., nrow = 1)
  
  ._ %>%
    as.list %>%
    c(list(sep = '_')) %>%
    do.call(what = paste) %>%
    sprintf(fmt = 'plots/dend_%s.pdf') %>%
    ggsave(plot = plot_out, width = 6.5, scale = 2)
})
