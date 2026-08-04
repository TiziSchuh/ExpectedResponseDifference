# Load packages and functions
library(tidyverse)
library(lavaan)
library(mirt)

# Load dataset
load("VIA_merged_unfiltered.Rda")


# Filter valid cases (i.e., analysis sample)
via <- filter(via, FILTER == 1)


# Sort data by study (2,3,4) and country (Germany, U.K.)
via <- dplyr::arrange(via, STUDY, COUNTRY)


# Recode negatively-keyed items
via <- via %>% mutate_at(
  vars(cap4, for4, gra4, equ2, int3, kin1, lea1, soc2, val2),
  .funs = list(
    R = ~recode(., .default = NA_real_,
                `1` = 5, `2` = 4, `3` = 3, `4` = 2, `5` = 1)) )


# Item names
item_names <- c(
  "cap2", "cap4_R", "for4_R", "gra4_R", "hop3", "zes3",
  "equ2_R", "equ3", "int3_R", "kin1_R", "mod4", "pru4",
  "jud4", "lea1_R", "ori3", "per2", "soc2_R", "val2_R"
)

# Use validation sample (no planned missingness)... can be changed later of course
data_val <- via[
  via$STUDY == 3,
  c("COUNTRY", item_names)
]

# Separate items and grouping variable
items_val <- data_val[, item_names]

group_val <- factor(data_val$COUNTRY)


#### Fit IRT-model with multipleGroup (required for empirical_ES) ----

# Specify the model
mod_syntax <- '
  Pos = cap2, cap4_R, for4_R, gra4_R, hop3, zes3
  Dep = equ2_R, equ3, int3_R, kin1_R, mod4, pru4
  Mas = jud4, lea1_R, ori3, per2, soc2_R, val2_R

  COV = Pos*Dep*Mas
'

mod <- mirt.model(
  mod_syntax,
  itemnames = items_val
)

# Estimate the model
mirt_mod <- multipleGroup(
  data     = items_val,
  model    = mod,
  group    = group_val,
  itemtype = "graded"
)

coef_mod <- coef(
  mirt_mod, 
  simplify = TRUE, 
  IRTpars = TRUE
  )

# WARNING: is unstable without acquiescence factor (non-pos definite latent cov-matrix; Correlation between Dep and Mas = 1)
# Reestimate with acquiescence factor

mod_acq_syntax <- '
  Pos = 1-6
  Dep = 7-12
  Mas = 13-18
  Acq = 1-18

  COV = Pos*Dep*Mas, Acq*Acq

  START = (1,5-6,8,11-13,15-16,a4,1),
          (2-4,7,9-10,14,17-18,a4,-1)

  FIXED = (1-18,a4)
'

mod_acq <- mirt.model(mod_acq_syntax)

# Estimation (with Metropolis-Hastings Robbins-Monro (MHRM) algorithm; better suited for models with >3 factors according to documentation)
mirt_mod_acq <- multipleGroup(
  data     = items_val,
  model    = mod_acq,
  group    = group_val,
  itemtype = "graded",
  method   = "MHRM"
)
# --> looks better

cf_acq <- coef(
  mirt_mod_acq,
  simplify = TRUE,
  IRTpars = FALSE # I set this as false here because it would calculate thresholds but these only make sense in unidimensional models (here with have factor + acq per item)
)


#### Calculate DIF-effect sizes ----

# Problem: only unidimensional models are supported...
# So below I estimate the three subscales separately and get effect sizes for them. First item of the scale is used as an anchor


#### For Factor Pos ----

items_Pos <- c(
  "cap2",
  "cap4_R",
  "for4_R",
  "gra4_R",
  "hop3",
  "zes3"
)


data_Pos <- items_val[, items_Pos]

group_Pos <- droplevels(factor(group_val))

# Specify a unidimensional model for the Pos scale
mod_Pos <- mirt.model("
  Pos = 1-6
")


# Define the single anchor item
anchor_Pos <- "cap2"


# Estimate the multiple-group model
#
# The anchor item defines the common metric.
# All other items are freely estimated across groups.
# The focal-group latent mean and variance are freely estimated.
mirt_mod_Pos <- multipleGroup(
  data       = data_Pos,
  model      = mod_Pos,
  group      = group_Pos,
  itemtype   = "graded",
  invariance = c(
    anchor_Pos,
    "free_means",
    "free_var"
  ),
  method  = "EM",
  verbose = TRUE
)

# Inspect the estimated item parameters
coef_Pos <- coef(
  mirt_mod_Pos,
  simplify = TRUE,
  IRTpars  = TRUE
)

coef_Pos


# Calculate item-level empirical DIF effect sizes
es_Pos_all <- as.data.frame(
  empirical_ES(
    mirt_mod_Pos,
    DIF = TRUE
  )
)

# Calculate test-level empirical DIF effect sizes
es_Pos_all_test <- as.data.frame(
  empirical_ES(
    mirt_mod_Pos,
    DIF = FALSE
  )
)


# Move the item names from the row names into a variable
es_Pos_all$item <- rownames(es_Pos_all)
rownames(es_Pos_all) <- NULL

es_Pos_all <- es_Pos_all[
  ,
  c("item", setdiff(names(es_Pos_all), "item"))
]


# Retain the five freely estimated target items
# The anchor has SIDS = UIDS = 0 by construction
es_Pos_report <- es_Pos_all[
  es_Pos_all$item != anchor_Pos,
  c("item", "SIDS", "UIDS", "SIDN", "UIDN")
]


#### Caculate SD_pool per non-anchor item and also calcualte SIDS/UIDS ---

reference_group <- levels(group_Pos)[1] # GERMAN
focal_group     <- levels(group_Pos)[2] # UK


# Function to calculate the pooled observed SD for one item
pooled_item_sd <- function(item, data, group,
                           reference_group, focal_group) {
  
  x_ref <- data[[item]][group == reference_group]
  x_foc <- data[[item]][group == focal_group]
  
  # Remove missing responses separately by group
  x_ref <- x_ref[!is.na(x_ref)]
  x_foc <- x_foc[!is.na(x_foc)]
  
  n_ref <- length(x_ref)
  n_foc <- length(x_foc)
  
  sd_ref <- sd(x_ref)
  sd_foc <- sd(x_foc)
  
  sd_pool <- sqrt(
    ((n_ref - 1) * sd_ref^2 +
       (n_foc - 1) * sd_foc^2) /
      (n_ref + n_foc - 2)
  )
  
  data.frame(
    item         = item,
    n_reference  = n_ref,
    n_focal      = n_foc,
    sd_reference = sd_ref,
    sd_focal     = sd_foc,
    sd_pooled    = sd_pool
  )
}


# Calculate the observed SDs for the five non-anchor items
observed_sd_Pos <- do.call(
  rbind,
  lapply(
    es_Pos_report$item,
    pooled_item_sd,
    data = data_Pos,
    group = group_Pos,
    reference_group = reference_group,
    focal_group = focal_group
  )
)


# Merge the observed SDs with the SIDS/UIDS results
es_Pos_final <- merge(
  es_Pos_report,
  observed_sd_Pos,
  by = "item",
  all.x = TRUE,
  sort = FALSE
)

# Standardized unsigned effect size:
# UIDS divided by pooled SD of observed 1-5 responses
es_Pos_final$UIDS_std <- with(
  es_Pos_final,
  UIDS / sd_pooled
)

# Standardized unsigned effect size (inferential):
# UIDS divided by pooled SD of observed 1-5 responses
es_Pos_final$UIDN_std <- with(
  es_Pos_final,
  UIDN / sd_pooled
)


# SIDS/UIDS
es_Pos_final$SIDS_UIDS_ratio <- with(
  es_Pos_final,
  ifelse(
    UIDS == 0,
    NA_real_,
    SIDS / UIDS
  )
)


# SIDN/UIDN
es_Pos_final$SIDN_UIDN_ratio <- with(
  es_Pos_final,
  ifelse(
    UIDN == 0,
    NA_real_,
    SIDN / UIDN
  )
)


es_Pos_final <- es_Pos_final[, c(
  "item",
  "SIDS",
  "UIDS",
  "SIDN",
  "UIDN",
  "sd_pooled",
  "UIDS_std",
  "SIDS_UIDS_ratio",
  "UIDN_std",
  "SIDN_UIDN_ratio",
  "sd_reference",
  "sd_focal"
)]


### Final item-level effect sizes
es_Pos_final



#### From test-level effect sizes above (es_Pos_all_test) get UTDS/SD_test_pool and STDS/UTDS ---

STDS_Pos <- es_Pos_all_test[1,2]
UTDS_Pos <- es_Pos_all_test[2,2]


### STDS/UTDS
STDS_UTDS_ratio_Pos <- STDS_Pos / UTDS_Pos


### UTDS/SD_test_pool
test_items_Pos <- items_Pos # all six Pos items, including the anchor

complete_Pos <- complete.cases(data_Pos[, test_items_Pos])

# Calculate observed sum scores
sumscore_Pos <- rowSums(
  data_Pos[complete_Pos, test_items_Pos, drop = FALSE]
)

group_sumscore_Pos <- droplevels(group_Pos[complete_Pos])

# Identify reference and focal groups
reference_group <- levels(group_sumscore_Pos)[1] # GERMAN
focal_group     <- levels(group_sumscore_Pos)[2] # UK

# Group-specific sum-score vectors
sumscore_ref <- sumscore_Pos[
  group_sumscore_Pos == reference_group
]

sumscore_foc <- sumscore_Pos[
  group_sumscore_Pos == focal_group
]

# Sample sizes and SDs
n_ref <- length(sumscore_ref)
n_foc <- length(sumscore_foc)

sd_test_ref <- sd(sumscore_ref)
sd_test_foc <- sd(sumscore_foc)

sd_test_pool <- sqrt(
  ((n_ref - 1) * sd_test_ref^2 +
     (n_foc - 1) * sd_test_foc^2) /
    (n_ref + n_foc - 2)
)

sd_test_pool

# Standardized unsigned test-level DIF effect
UTDS_std_Pos <- UTDS_Pos / sd_test_pool



#### For Factor Dep ----

items_Dep <- c(
  "equ2_R", 
  "equ3", 
  "int3_R", 
  "kin1_R", 
  "mod4", 
  "pru4"
)


data_Dep <- items_val[, items_Dep]

group_Dep <- droplevels(factor(group_val))

# Specify a unidimensional model for the Dep scale
mod_Dep <- mirt.model("
  Dep = 1-6
")


# Define the single anchor item
anchor_Dep <- "equ2_R"


# Estimate the multiple-group model
#
# The anchor item defines the common metric.
# All other items are freely estimated across groups.
# The focal-group latent mean and variance are freely estimated.
mirt_mod_Dep <- multipleGroup(
  data       = data_Dep,
  model      = mod_Dep,
  group      = group_Dep,
  itemtype   = "graded",
  invariance = c(
    anchor_Dep,
    "free_means",
    "free_var"
  ),
  method  = "EM",
  verbose = TRUE
)

# Inspect the estimated item parameters
coef_Dep <- coef(
  mirt_mod_Dep,
  simplify = TRUE,
  IRTpars  = TRUE
)

coef_Dep


# Calculate item-level empirical DIF effect sizes
es_Dep_all <- as.data.frame(
  empirical_ES(
    mirt_mod_Dep,
    DIF = TRUE
  )
)

# Calculate test-level empirical DIF effect sizes
es_Dep_all_test <- as.data.frame(
  empirical_ES(
    mirt_mod_Dep,
    DIF = FALSE
  )
)


# Move the item names from the row names into a variable
es_Dep_all$item <- rownames(es_Dep_all)
rownames(es_Dep_all) <- NULL

es_Dep_all <- es_Dep_all[
  ,
  c("item", setdiff(names(es_Dep_all), "item"))
]


# Retain the five freely estimated target items
# The anchor has SIDS = UIDS = 0 by construction
es_Dep_report <- es_Dep_all[
  es_Dep_all$item != anchor_Dep,
  c("item", "SIDS", "UIDS", "SIDN", "UIDN")
]


#### Caculate SD_pool per non-anchor item and also calcualte SIDS/UIDS ---

reference_group <- levels(group_Dep)[1] # GERMAN
focal_group     <- levels(group_Dep)[2] # UK


# Function to calculate the pooled observed SD for one item
pooled_item_sd <- function(item, data, group,
                           reference_group, focal_group) {
  
  x_ref <- data[[item]][group == reference_group]
  x_foc <- data[[item]][group == focal_group]
  
  # Remove missing responses separately by group
  x_ref <- x_ref[!is.na(x_ref)]
  x_foc <- x_foc[!is.na(x_foc)]
  
  n_ref <- length(x_ref)
  n_foc <- length(x_foc)
  
  sd_ref <- sd(x_ref)
  sd_foc <- sd(x_foc)
  
  sd_pool <- sqrt(
    ((n_ref - 1) * sd_ref^2 +
       (n_foc - 1) * sd_foc^2) /
      (n_ref + n_foc - 2)
  )
  
  data.frame(
    item         = item,
    n_reference  = n_ref,
    n_focal      = n_foc,
    sd_reference = sd_ref,
    sd_focal     = sd_foc,
    sd_pooled    = sd_pool
  )
}


# Calculate the observed SDs for the five non-anchor items
observed_sd_Dep <- do.call(
  rbind,
  lapply(
    es_Dep_report$item,
    pooled_item_sd,
    data = data_Dep,
    group = group_Dep,
    reference_group = reference_group,
    focal_group = focal_group
  )
)


# Merge the observed SDs with the SIDS/UIDS results
es_Dep_final <- merge(
  es_Dep_report,
  observed_sd_Dep,
  by = "item",
  all.x = TRUE,
  sort = FALSE
)

# Standardized unsigned effect size:
# UIDS divided by pooled SD of observed 1-5 responses
es_Dep_final$UIDS_std <- with(
  es_Dep_final,
  UIDS / sd_pooled
)

# Standardized unsigned effect size (inferential):
# UIDS divided by pooled SD of observed 1-5 responses
es_Dep_final$UIDN_std <- with(
  es_Dep_final,
  UIDN / sd_pooled
)


# SIDS/UIDS
es_Dep_final$SIDS_UIDS_ratio <- with(
  es_Dep_final,
  ifelse(
    UIDS == 0,
    NA_real_,
    SIDS / UIDS
  )
)


# SIDN/UIDN
es_Dep_final$SIDN_UIDN_ratio <- with(
  es_Dep_final,
  ifelse(
    UIDN == 0,
    NA_real_,
    SIDN / UIDN
  )
)


es_Dep_final <- es_Dep_final[, c(
  "item",
  "SIDS",
  "UIDS",
  "SIDN",
  "UIDN",
  "sd_pooled",
  "UIDS_std",
  "SIDS_UIDS_ratio",
  "UIDN_std",
  "SIDN_UIDN_ratio",
  "sd_reference",
  "sd_focal"
)]


### Final item-level effect sizes
es_Dep_final



#### From test-level effect sizes above (es_Dep_all_test) get UTDS/SD_test_pool and STDS/UTDS ---

STDS_Dep <- es_Dep_all_test[1,2]
UTDS_Dep <- es_Dep_all_test[2,2]


### STDS/UTDS
STDS_UTDS_ratio_Dep <- STDS_Dep / UTDS_Dep


### UTDS/SD_test_pool
test_items_Dep <- items_Dep # all six Dep items, including the anchor

complete_Dep <- complete.cases(data_Dep[, test_items_Dep])

# Calculate observed sum scores
sumscore_Dep <- rowSums(
  data_Dep[complete_Dep, test_items_Dep, drop = FALSE]
)

group_sumscore_Dep <- droplevels(group_Dep[complete_Dep])

# Identify reference and focal groups
reference_group <- levels(group_sumscore_Dep)[1] # GERMAN
focal_group     <- levels(group_sumscore_Dep)[2] # UK

# Group-specific sum-score vectors
sumscore_ref <- sumscore_Dep[
  group_sumscore_Dep == reference_group
]

sumscore_foc <- sumscore_Dep[
  group_sumscore_Dep == focal_group
]

# Sample sizes and SDs
n_ref <- length(sumscore_ref)
n_foc <- length(sumscore_foc)

sd_test_ref <- sd(sumscore_ref)
sd_test_foc <- sd(sumscore_foc)

sd_test_pool <- sqrt(
  ((n_ref - 1) * sd_test_ref^2 +
     (n_foc - 1) * sd_test_foc^2) /
    (n_ref + n_foc - 2)
)

sd_test_pool

# Standardized unsigned test-level DIF effect
UTDS_std_Dep <- UTDS_Dep / sd_test_pool


#### For Factor Mas ----

items_Mas <- c(
  "jud4", 
  "lea1_R", 
  "ori3",
  "per2", 
  "soc2_R", 
  "val2_R"
)


data_Mas <- items_val[, items_Mas]

group_Mas <- droplevels(factor(group_val))

# Specify a unidimensional model for the Mas scale
mod_Mas <- mirt.model("
  Mas = 1-6
")


# Define the single anchor item
anchor_Mas <- "jud4"


# Estimate the multiple-group model
#
# The anchor item defines the common metric.
# All other items are freely estimated across groups.
# The focal-group latent mean and variance are freely estimated.
mirt_mod_Mas <- multipleGroup(
  data       = data_Mas,
  model      = mod_Mas,
  group      = group_Mas,
  itemtype   = "graded",
  invariance = c(
    anchor_Mas,
    "free_means",
    "free_var"
  ),
  method  = "EM",
  verbose = TRUE
)

# Inspect the estimated item parameters
coef_Mas <- coef(
  mirt_mod_Mas,
  simplify = TRUE,
  IRTpars  = TRUE
)

coef_Mas


# Calculate item-level empirical DIF effect sizes
es_Mas_all <- as.data.frame(
  empirical_ES(
    mirt_mod_Mas,
    DIF = TRUE
  )
)

# Calculate test-level empirical DIF effect sizes
es_Mas_all_test <- as.data.frame(
  empirical_ES(
    mirt_mod_Mas,
    DIF = FALSE
  )
)


# Move the item names from the row names into a variable
es_Mas_all$item <- rownames(es_Mas_all)
rownames(es_Mas_all) <- NULL

es_Mas_all <- es_Mas_all[
  ,
  c("item", setdiff(names(es_Mas_all), "item"))
]


# Retain the five freely estimated target items
# The anchor has SIDS = UIDS = 0 by construction
es_Mas_report <- es_Mas_all[
  es_Mas_all$item != anchor_Mas,
  c("item", "SIDS", "UIDS", "SIDN", "UIDN")
]


#### Caculate SD_pool per non-anchor item and also calcualte SIDS/UIDS ---

reference_group <- levels(group_Mas)[1] # GERMAN
focal_group     <- levels(group_Mas)[2] # UK


# Function to calculate the pooled observed SD for one item
pooled_item_sd <- function(item, data, group,
                           reference_group, focal_group) {
  
  x_ref <- data[[item]][group == reference_group]
  x_foc <- data[[item]][group == focal_group]
  
  # Remove missing responses separately by group
  x_ref <- x_ref[!is.na(x_ref)]
  x_foc <- x_foc[!is.na(x_foc)]
  
  n_ref <- length(x_ref)
  n_foc <- length(x_foc)
  
  sd_ref <- sd(x_ref)
  sd_foc <- sd(x_foc)
  
  sd_pool <- sqrt(
    ((n_ref - 1) * sd_ref^2 +
       (n_foc - 1) * sd_foc^2) /
      (n_ref + n_foc - 2)
  )
  
  data.frame(
    item         = item,
    n_reference  = n_ref,
    n_focal      = n_foc,
    sd_reference = sd_ref,
    sd_focal     = sd_foc,
    sd_pooled    = sd_pool
  )
}


# Calculate the observed SDs for the five non-anchor items
observed_sd_Mas <- do.call(
  rbind,
  lapply(
    es_Mas_report$item,
    pooled_item_sd,
    data = data_Mas,
    group = group_Mas,
    reference_group = reference_group,
    focal_group = focal_group
  )
)


# Merge the observed SDs with the SIDS/UIDS results
es_Mas_final <- merge(
  es_Mas_report,
  observed_sd_Mas,
  by = "item",
  all.x = TRUE,
  sort = FALSE
)

# Standardized unsigned effect size:
# UIDS divided by pooled SD of observed 1-5 responses
es_Mas_final$UIDS_std <- with(
  es_Mas_final,
  UIDS / sd_pooled
)

# Standardized unsigned effect size (inferential):
# UIDS divided by pooled SD of observed 1-5 responses
es_Mas_final$UIDN_std <- with(
  es_Mas_final,
  UIDN / sd_pooled
)


# SIDS/UIDS
es_Mas_final$SIDS_UIDS_ratio <- with(
  es_Mas_final,
  ifelse(
    UIDS == 0,
    NA_real_,
    SIDS / UIDS
  )
)


# SIDN/UIDN
es_Mas_final$SIDN_UIDN_ratio <- with(
  es_Mas_final,
  ifelse(
    UIDN == 0,
    NA_real_,
    SIDN / UIDN
  )
)


es_Mas_final <- es_Mas_final[, c(
  "item",
  "SIDS",
  "UIDS",
  "SIDN",
  "UIDN",
  "sd_pooled",
  "UIDS_std",
  "SIDS_UIDS_ratio",
  "UIDN_std",
  "SIDN_UIDN_ratio",
  "sd_reference",
  "sd_focal"
)]


### Final item-level effect sizes
es_Mas_final



#### From test-level effect sizes above (es_Mas_all_test) get UTDS/SD_test_pool and STDS/UTDS ---

STDS_Mas <- es_Mas_all_test[1,2]
UTDS_Mas <- es_Mas_all_test[2,2]


### STDS/UTDS
STDS_UTDS_ratio_Mas <- STDS_Mas / UTDS_Mas


### UTDS/SD_test_pool
test_items_Mas <- items_Mas # all six Mas items, including the anchor

complete_Mas <- complete.cases(data_Mas[, test_items_Mas])

# Calculate observed sum scores
sumscore_Mas <- rowSums(
  data_Mas[complete_Mas, test_items_Mas, drop = FALSE]
)

group_sumscore_Mas <- droplevels(group_Mas[complete_Mas])

# Identify reference and focal groups
reference_group <- levels(group_sumscore_Mas)[1] # GERMAN
focal_group     <- levels(group_sumscore_Mas)[2] # UK

# Group-specific sum-score vectors
sumscore_ref <- sumscore_Mas[
  group_sumscore_Mas == reference_group
]

sumscore_foc <- sumscore_Mas[
  group_sumscore_Mas == focal_group
]

# Sample sizes and SDs
n_ref <- length(sumscore_ref)
n_foc <- length(sumscore_foc)

sd_test_ref <- sd(sumscore_ref)
sd_test_foc <- sd(sumscore_foc)

sd_test_pool <- sqrt(
  ((n_ref - 1) * sd_test_ref^2 +
     (n_foc - 1) * sd_test_foc^2) /
    (n_ref + n_foc - 2)
)

sd_test_pool

# Standardized unsigned test-level DIF effect
UTDS_std_Mas <- UTDS_Mas / sd_test_pool
