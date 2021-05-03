library(gridExtra)
library(cowplot)
library(dplyr)
library(magrittr)
library(ggplot2)
library(gtable)
library(grid)

source('r/funcs.R')
source('r/load_data.R')


box_list <- seq.int(1, 4) %>%
  set_names(
    c(
      'rosaria',
      'bodega',
      'solano',
      'westchester'
    )
  )

plot_out <-
  lapply(
    REF_LIST$top$betti0[
      c('rosaria',
        'bodega',
        'solano',
        'westchester')
    ],
    function(._){
    
      data_split <-
        DATA_GILL %>%
        select(IND_LIST[[._]]) %>%
        split(DATA_GILL$location)
      data_split <- data_split[
        NAME_LIST[
          c('rosaria',
            'bodega',
            'solano',
            'westchester')
        ]
      ]
      
      plot_list <- lapply(names(data_split), function(._){
        data_split[[._]] %>%
          cor %>%
          COR.DEND %>%
          GG.DEND(labels = FALSE, outline = TRUE, xlim = c(-.6, 1), xlab = NULL)
      })
      
      box <- box_list[last(._)]
      title <- textGrob(NAME_LIST[last(._)], gp = gpar(fontsize = 10))
      padding <- unit(0.5,"line")
      
      arrangeGrob(grobs = plot_list, nrow = 1) %>%
        gtable_add_rows(
          heights = grobHeight(title) + padding, pos = 0
        ) %>%
        gtable_add_grob(
          list(title),
          t = 1, l = 1, r = 4
        ) %>%
        gtable_add_grob(
          rectGrob(
            width = .98,
            height = .98,
            gp = gpar(lwd = 5, fill = NA, col = 'red')),
          name = 'box',
          2, box, 2, box
        ) %>% 
        gtable_add_grob(
          rectGrob(
            gp = gpar(lwd = 3, fill = NA)),
          name = 'box_outer',
          1, 1, 2, 4
        )
    }
  ) %>%
  arrangeGrob(
    grobs = .,
    nrow = 2,
    layout_matrix = matrix(seq.int(1, 4), nrow = 2)
  )
    
'plots/dend_4x4.pdf' %>%
ggsave(plot = plot_out, width = 6.5, height = 1.7, scale = 1.4)

