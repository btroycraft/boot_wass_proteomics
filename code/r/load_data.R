library(tidyverse)

DATA_GILL <- read.csv('data/cDIA_MXLSAKBH-Exp1-2-3-4_Gill_r_format.csv', stringsAsFactors = FALSE)

NAME_LIST <- c(
  total = 'Total',
  rosaria = 'Laguna de la Bocana del Rosaria',
  bodega = 'Bodega Harbor',
  solano = 'Lake Solano',
  westchester = 'Westchester Lagoon'
)

# List of references for cases where the name of the branch is needed
# rapply(REF_LIST, how = 'replace', . %>% IND_LIST[[.]] ) should
# reproduce IND_LIST

REF_LIST <- c(
  list(all = 'all'),
  list(
    max = 'max',
    top = list(
      betti0 = c('top', 'betti0'),
      betti1 = c('top', 'betti1')
    )
  ) %>%
    rapply(how = 'replace', function(._)
      lapply(names(NAME_LIST), . %>% c(._, .)) %>%
        set_names(names(NAME_LIST))
    )
)

# Loads a nested list of indices, corresponding to all selection steps
# Only the first most separating subset is considered for each

IND_LIST <- c(
  list(all = names(DATA_GILL)[!(names(DATA_GILL) %in% c('sample', 'location'))]),
  REF_LIST[c('max', 'top')] %>%
    rapply(how = 'replace', . %>%
      c(list(sep = '_')) %>%
      do.call(paste, .) %>%
      paste0('data/protein_selection/ind_', ., '.csv') %>%
      read.csv(stringsAsFactors = FALSE) %>%
      .[1, ] %>%
      select(!val) %>%
      as.character
    )
)
