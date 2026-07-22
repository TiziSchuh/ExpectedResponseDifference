# Load packages and functions
library(tidyverse)
library(lavaan)
library(osfr)

# Download data set from Partsch et al. (2024) from their OSF repository
osf_retrieve_node("https://osf.io/pdntc/") %>%
  osf_ls_files(., path = "1_Data+analyses/2_Scale-properties+scale-validation_Tables4-6", 
               pattern = "VIA_merged_unfiltered.Rda", n_max = Inf) %>%
  osf_download(., conflicts = "overwrite")

# Load dataset
load("VIA_merged_unfiltered.Rda")

# Filter valid cases (i.e., analysis sample) from Partsch et al.'s "validation sample" and sort them by country
via <- filter(via, FILTER == 1, STUDY == 3) %>%
  arrange(COUNTRY)

# Specify model object
fa.model <- 'Pos =~ cap2 + cap4 + for4 + gra4 + hop3 + zes3 
\n Dep =~ equ2 + equ3 + int3 + kin1 + mod4 + pru4 
\n Mas =~ jud4 + lea1 + ori3 + per2 + soc2 + val2 
\n aqui =~ cap2 + 1*cap4 + 1*for4 + 1*gra4 + 1*hop3 + 1*zes3 
+ 1*equ2 + 1*equ3 + 1*int3 + 1*kin1 + 1*mod4 + 1*pru4 + 1*jud4 
+ 1*lea1 + 1*ori3 + 1*per2 + 1*soc2 + 1*val2 
\n aqui ~~ 0*Pos + 0*Dep + 0*Mas'

# Fit multi-group model (i.e., CONFIGURAL invariance measurement model) to yield group-specific parameters
fa.model.mg <-lavaan::cfa(
  model = fa.model,
  data = via,
  meanstructure = T,
  estimator = "MLR",
  missing = "fiml",
  group = "COUNTRY") # Germany: COUNTRY == 1; UK: COUNTRY == 2 

# Inspect fitted multi-group model
lavaan::summary(fa.model.mg, fit.measures = T, standardized = T, remove.unused = F)

# Fit single-group ("global") model to yield joint parameters
fa.model.sg  <-lavaan::cfa(
  model = fa.model, 
  data = via,
  meanstructure = T,
  estimator = "MLR",
  missing = "fiml")

# Inspect fitted single-group model
lavaan::summary(fa.model.sg, fit.measures = T, standardized = T, remove.unused = F)
