library(tidyverse)

DATA_GILL <- read.csv('data/cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv', stringsAsFactors = FALSE)

NAME_LIST <- c(
  'Total' = 'total',
  'Laguna de la Bocana del Rosaria' = 'rosaria',
  'Bodega Harbor' = 'bodega',
  'Lake Solano' = 'solano',
  'Westchester Lagoon' = 'westchester'
)

IND_LIST <- list(
  all = names(DATA_GILL)[!(names(DATA_GILL) %in% c('sample', 'location'))],
  max = lapply(NAME_LIST, . %>%
    sprintf(fmt = 'data/ind_max_%s.csv') %>%
    read.csv(stringsAsFactors = FALSE) %>%
    .[1, ] %>%
    select(!val) %>%
    as.character
  ),
  top = list(
    betti0 = lapply(NAME_LIST, . %>%
      sprintf(fmt = 'data/ind_top_%s_betti0.csv') %>%
      read.csv(stringsAsFactors = FALSE) %>%
      .[1, ] %>%
      select(!val) %>%
      as.character
    ),
    betti1 = lapply(NAME_LIST, . %>%
      sprintf(fmt = 'data/ind_top_%s_betti1.csv') %>%
      read.csv(stringsAsFactors = FALSE) %>%
      .[1, ] %>%
      select(!val) %>%
      as.character
    )
  )
)

REF_LIST <- list(
  all = 'all',
  max =
    names(IND_LIST$max) %>%
    set_names(., .) %>%
    as.list %>%
    lapply(. %>% c('max', .) ),
  top = list(
    betti0 =
      names(IND_LIST$top$betti0) %>%
      set_names(., .) %>%
      lapply(. %>% c('top', 'betti0', .) ),
    betti1 =
      names(IND_LIST$top$betti1) %>%
      set_names(., .) %>%
      lapply(. %>% c('top', 'betti1', .) )
  )
)

DATA_LIST <- rapply(IND_LIST, how = 'replace', . %>% select(DATA_GILL, .) )