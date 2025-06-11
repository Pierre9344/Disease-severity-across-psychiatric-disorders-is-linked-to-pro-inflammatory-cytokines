library(targets, quietly = T)
library(tarchetypes, quietly = T)
library(crew, quietly = T) # parallel computing
library(autometric)
# library(proffer)
# to figure the steps that are slow
# proffer::pprof(tar_make(callr_function = NULL))
library(qs2) # more efficient save format
library(quarto)

library(tidyverse, quietly = T)
library(magrittr, quietly = T)
library(here, quietly = T)
# library(janitor)
library(stringr, quietly = T)
library(R6, quietly = T)
library(utils, quietly = T) # useful for getting combinaison when comparing groups:
# utils::combn(LETTERS[1:3], 2, simplify = F) %>% sapply(function(x) paste(x, collapse = "_")

library(ggplot2, quietly = T)
library(ggiraph, quietly = T)
library(ggh4x, quietly = T)
library(ggrepel, quietly = T)

library(openxlsx, quietly = T)
library(lme4)
library(org.Hs.eg.db)

# needed here for rendering the Proteome.qmd document
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(DOSE)
