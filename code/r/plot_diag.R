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
        RIPSER(q = 1, t = 2) %>%
        {sapply(seq(0, 2, length.out = 10^3), function(._){
          PBETTI(., r = ._, q = 1)
        })} %>%
        {data.frame(
          cor = seq(0, 2, length.out = 10^3),
          pbetti = .,
          location = ._
        )}
    }) %>%
    do.call(what = rbind)
  
  plot_aligned <-
    lapply(names(data_split), function(._){
      
      plot1 <- 
        data_split[[._]] %>%
        {1-cor(.)} %>%
        RIPSER(q = 2, t = 2) %>%
        filter(Dim >= 1) %>%
        mutate(Dim = as.character(Dim)) %>%
        ggplot() +
          geom_point(
            aes(
              x = Birth,
              y = Death,
              group = Dim,
              color = Dim,
              shape = Dim
            ),
            show.legend = FALSE
          ) +
          geom_abline(slope = 1, intercept= 0) +
          coord_fixed(ratio = 1, xlim = c(0, 1.7), ylim = c(0, 1.7)) +
          theme_classic() +
          theme(
            plot.title = element_text(hjust = .5),
            panel.grid.major.x = element_line(color = 'gray91', size = .4),
            panel.grid.major.y = element_line(color = 'gray91', size = .4),
            panel.grid.minor.x = element_line(color = 'gray91', size = .2),
            panel.grid.minor.y = element_line(color = 'gray91', size = .2)
          ) +
          scale_x_continuous(
            breaks = c(0, .5, 1, 1.5),
            labels = c('1.0', '0.5', '0', '-0.5'),
            position = 'bottom'
          ) +
          scale_y_continuous(
            breaks = c(0, .5, 1, 1.5),
            labels = c('1.0', '0.5', '0', '-0.5'),
            position = 'left'
          ) +
          scale_color_manual(values = c('red', 'blue'))
      
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
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = .5),
          panel.grid.major.x = element_line(color = 'gray91', size = .4),
          panel.grid.major.y = element_line(color = 'gray91', size = .4),
          panel.grid.minor.x = element_line(color = 'gray91', size = .2),
          panel.grid.minor.y = element_line(color = 'gray91', size = .2)
        ) +
        ylab('Betti-1') +
        xlab('Correlation') +
        coord_cartesian(xlim = c(0, 1.7), ylim = c(0, 25))
      
      list(plot2, plot1)
    }) %>%
    do.call(what = c) %>%
    plot_grid(
      plotlist = .,
      align = 'v',
      nrow = 2,
      byrow = FALSE,
      rel_heights = c(1, 2))
  
  plot_out <-
    {ggplot(data = data.frame(Dim = c('1', '2'))) +
      geom_point(
        aes(
          x = c(1, 2),
          y = c(1, 2),
          group = Dim,
          color = Dim,
          shape = Dim,
        ),
        show.legend = TRUE
      ) +
      scale_color_manual(
        values = c('red', 'blue')
      ) +
      theme(
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"),
        legend.position = c(.5, .7)
      )
    } %>%
    get_legend %>%
    list(plot_aligned, ggplot(), .) %>%
    arrangeGrob(grobs = ., nrow = 2, widths = c(1, 1, 1, 1, .3),
      layout_matrix = matrix(nrow = 2, ncol = 5, byrow = TRUE,
        c(
          1, 1, 1, 1, 2,
          1, 1, 1, 1, 3
        )                     
      )
    )
  
  ._ %>%
    as.list %>%
    c(list(sep = '_')) %>%
    do.call(what = paste) %>%
    sprintf(fmt = 'plots/diag_%s.pdf') %>%
    ggsave(plot = plot_out, width = 6.5, height = 2.3, scale = 2.1)
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
        RIPSER(q = 1, t = 2) %>%
        {sapply(seq(0, 2, length.out = 10^3), function(._){
          PBETTI(., r = ._, q = 1)
        })} %>%
        {data.frame(
          cor = seq(0, 2, length.out = 10^3),
          pbetti = .,
          location = ._
        )}
    }) %>%
    do.call(what = rbind)
  
  plot_aligned <-
    lapply(names(data_split), function(._){
      
      plot1 <- 
        data_split[[._]] %>%
        {1-cor(.)} %>%
        RIPSER(q = 2, t = 2) %>%
        filter(Dim >= 1) %>%
        mutate(Dim = as.character(Dim)) %>%
        ggplot() +
        geom_point(
          aes(
            x = Birth,
            y = Death,
            group = Dim,
            color = Dim,
            shape = Dim
          ),
          show.legend = FALSE
        ) +
        geom_abline(slope = 1, intercept= 0) +
        coord_fixed(ratio = 1, xlim = c(0, 1.7), ylim = c(0, 1.7)) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust = .5),
          panel.grid.major.x = element_line(color = 'gray91', size = .4),
          panel.grid.major.y = element_line(color = 'gray91', size = .4),
          panel.grid.minor.x = element_line(color = 'gray91', size = .2),
          panel.grid.minor.y = element_line(color = 'gray91', size = .2)
        ) +
        scale_x_continuous(
          breaks = c(0, .5, 1, 1.5),
          labels = c('1.0', '0.5', '0', '-0.5'),
          position = 'bottom'
        ) +
        scale_y_continuous(
          breaks = c(0, .5, 1, 1.5),
          labels = c('1.0', '0.5', '0', '-0.5'),
          position = 'left'
        ) +
        scale_color_manual(values = c('red', 'blue'))
      
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
          axis.text.x = element_blank(),
          plot.title = element_text(hjust = .5),
          panel.grid.major.x = element_line(color = 'gray91', size = .4),
          panel.grid.major.y = element_line(color = 'gray91', size = .4),
          panel.grid.minor.x = element_line(color = 'gray91', size = .2),
          panel.grid.minor.y = element_line(color = 'gray91', size = .2)
        ) +
        ylab('Betti-1') +
        xlab('Correlation') +
        coord_cartesian(xlim = c(0, 1.7), ylim = c(0, 15))
      
      list(plot2, plot1)
    }) %>%
    do.call(what = c) %>%
    plot_grid(
      plotlist = .,
      align = 'v',
      nrow = 2,
      byrow = FALSE,
      rel_heights = c(1, 2))
  
  plot_out <-
    {ggplot(data = data.frame(Dim = c('1', '2'))) +
        geom_point(
          aes(
            x = c(1, 2),
            y = c(1, 2),
            group = Dim,
            color = Dim,
            shape = Dim,
          ),
          show.legend = TRUE
        ) +
        scale_color_manual(
          values = c('red', 'blue')
        ) +
        theme(
          legend.background = element_blank(),
          legend.box.background = element_rect(colour = "black"),
          legend.position = c(.5, .7)
        )
    } %>%
    get_legend %>%
    list(plot_aligned, ggplot(), .) %>%
    arrangeGrob(grobs = ., nrow = 2, widths = c(1, 1, 1, 1, .3),
                layout_matrix = matrix(nrow = 2, ncol = 5, byrow = TRUE,
                                       c(
                                         1, 1, 1, 1, 2,
                                         1, 1, 1, 1, 3
                                       )                     
                )
    )
  
  ._ %>%
    as.list %>%
    c(list(sep = '_')) %>%
    do.call(what = paste) %>%
    sprintf(fmt = 'plots/diag_%s.pdf') %>%
    ggsave(plot = plot_out, width = 6.5, height = 2.3, scale = 2.1)
})
