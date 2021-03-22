library(tidyverse)

source('funcs.R')

data_gill <- read.csv('../data/cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv')

name_list <-
  c('total', 'bodega', 'rosaria', 'solano', 'westchester') %>%
  set_names(c('Total', levels(data_gill$location)))

data_list <- local({
  data <-
    data_gill %>%
    select(!any_of(c('sample', 'location'))) %>%
    split(data_gill$location)
  
  list(
    all = data,
    max = {
      {
        name_list %>%
        lapply(
          . %>%
          sprintf(fmt = '../data/ind_max_%s.csv') %>%
          read.csv %>%
          .[1, -1] %>%
          unlist %>%
          as.character
        )
      } %>%
      lapply(
        (. %>% lapply(data, . %>% select)
      )
      
      lapply(. %>%
        sprintf('../data/ind_max_%s.csv', .) %>%
        read.csv %>%
        extract(1, -1) %>%
        unlist %>%
        as.character
      ) %>%
      lapply(. %>% {lapply(data, subset, select = .)})
      )
        
        
    }
  name_list 
    lapply(subset, select =
      name_list %>%
      lapply(. %>%
        sprintf('../data/ind_top_%s_betti0.csv', .) %>%
        read.csv %>%
        extract(1, -1) %>%
        unlist %>%
        as.character
      )
    split(f = data_gill$location)
  
  data_list %>%
    lapply(. %>%
      lapply(subset, select =
        name_list %>%
        lapply(. %>%
          sprintf('../data/ind_top_%s_betti0.csv', .) %>%
          read.csv %>%
          extract(1, -1) %>%
          unlist %>%
          as.character %.%
        )
      )
  lapply()
  name_list[1]%>%
    sprintf(fmt = '../data/ind_top_%s_betti0.csv') %>%
    read.csv %>%
    extract(1, -1) %>% 
    unlist %>%
    as.character %>%
    ?extract
    lapply(data_list, . %>% subset(select = .))
    
  lapply(data_list, subset, select = 
      name_list[1]%>%
      sprintf(fmt = '../data/ind_top_%s_betti0.csv') %>%
      read.csv %>%
      extract(1, -1) %>% 
      unlist %>%
      as.character
  max =
    name_list %>%
    lapply(. %>%
      sprintf(fmt = '../data/ind_top_%s_betti0.csv') %>%
      read.csv %>%
      extract(1, -1) %>% 
      unlist %>%
      as.character #%>%
      #lapply(data_list, subset, select = .)
    )
    local({
    file_list <- lapply(name_list, sprintf, fmt = '../data/ind_top_%s_betti0.csv')
    lapply(file_list, function(file){
      ind <- as.character(unlist(read.csv(sprintf(file))[1, -1]))
      lapply(data_list, subset, select = ind)
    })
  })
  {list(
    all = .,
    max = local({
      file_list <- lapply(name_list, sprintf, fmt = '../data/ind_top_%s_betti0.csv')
      lapply(file_list, function(file){
        ind <- as.character(unlist(read.csv(sprintf(file))[1, -1]))
        lapply(data_list, subset, select = ind)
      })
    }),
    top = list(
      betti0 = local({
        file_list <- lapply(name_list, sprintf, fmt = '../data/ind_top_%s_betti0.csv')
        lapply(file_list, function(file){
          ind <- as.character(unlist(read.csv(sprintf(file))[1, -1]))
          lapply(data_list, subset, select = ind)
        })
      }),
      betti1 = local({
        file_list <- lapply(name_list, function(name) sprintf('../data/ind_top_%s_betti1.csv', name) )
        lapply(file_list, function(file){
          ind <- as.character(unlist(read.csv(sprintf(file))[1, -1]))
          lapply(data_list, subset, select = ind)
        })
      })
    )
)

ind_list <- list(
  max = local({
    file_list <- lapply(name_list, function(name) sprintf('../data/ind_max_%s.csv', name) )
    return( lapply(file_list, function(file) as.character(unlist(read.csv(sprintf(file))[1, -1])) ) )
  }),
  top = list(
    betti0 = local({
      file_list <- lapply(name_list, function(name) sprintf('../data/ind_top_%s_betti0.csv', name) )
      return( lapply(file_list, function(file) as.character(unlist(read.csv(sprintf(file))[1, -1])) ) )
    }),
    betti1 = local({
      file_list <- lapply(name_list, function(name) sprintf('../data/ind_top_%s_betti1.csv', name) )
      return( lapply(file_list, function(file) as.character(unlist(read.csv(sprintf(file))[1, -1])) ) )
    })
  )
)

dist_list <- local({
  max <- lapply(names(data_list), function(name){
    data <- subset(data_list[[name]], select = ind_list$max[[name]])
    1 - cor(data)
  })
  top <- list(
    betti0 = lapply(names(data_list), function(name){
      data <- subset(data_list[[name]], select = ind_list$top$betti0[[name]])
      1 - cor(data)
    }),
    betti1 = lapply(names(data_list), function(name){
      data <- subset(data_list[[name]], select = ind_list$top$betti1[[name]])
      1 - cor(data)
    })
  )
  names(max) = names(top$betti0) = names(top$betti1) = names(data_list)
  return( list( max = max, top = top) )    
})

get.ripser <- function(D, q = 1L, t = 2){
  diag <- as.data.frame(ripserr::vietoris_rips(as.dist(D), max_dim = q, threshold = t))
  names(diag) <- c('Dim', 'Birth', 'Death')
  return( diag )
}

sub.diag <- function(diag, q){
  diag <- subset(diag, Dim == q)
  row.names(diag) <- seq_len(nrow(diag))
  return( diag )
}