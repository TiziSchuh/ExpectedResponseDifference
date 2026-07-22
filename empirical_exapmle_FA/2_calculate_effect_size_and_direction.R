## =====================================================================
## Step 2: Calculate per-item effect size and direction
## =====================================================================

# source("0_fit_models.R") # run if step0-script has not been executed yet
# source("1_fetch_parameters.R") # run if step1-script has not been executed yet
download.file(
  "https://raw.githubusercontent.com/TiziSchuh/ExpectedDifferenceMeasures/main/lavaan_edm.R",
  destfile = "lavaan_edm.R", mode = "wb")
source("lavaan_edm.R")

## ---------------------------------------------------------------------
## 0) Align items and sanity checks
## ---------------------------------------------------------------------
items <- rownames(step1$single$tau)

stopifnot(
  identical(rownames(step1$ger$Lambda), items),
  identical(rownames(step1$uk$Lambda),  items),
  identical(colnames(step1$ger$Lambda), colnames(step1$single$Lambda)),
  identical(colnames(step1$uk$Lambda),  colnames(step1$single$Lambda)),
  all(items %in% names(step1$s_pool))
)

s_pool <- step1$s_pool[items]
p_ger  <- step1$p["ger"]
p_uk   <- step1$p["uk"]

## ---------------------------------------------------------------------
## 1) Prepare wrapper:
## Per-item expected difference between the group-specific model and 
## the single-group model, integrating over the *group's* latent
## distribution (kappa_g, Sigma_g):
##   exp_response_diff(Lambda_g,i, Lambda_i, tau_g,i, tau_i,
##                     kappa_g, Sigma_g, diffType)
## ---------------------------------------------------------------------
ed_per_item <- function(grp, diffType) {
  vapply(items, function(i) {
    exp_response_diff(
      loading1   = step1[[grp]]$Lambda[i, ],
      loading2   = step1$single$Lambda[i, ],
      intercept1 = step1[[grp]]$tau[i, 1],
      intercept2 = step1$single$tau[i, 1],
      latentMean = step1[[grp]]$kappa,
      latentCov  = step1[[grp]]$Sigma,
      diffType   = diffType
    )
  }, numeric(1))
}

## ---------------------------------------------------------------------
## 2) Calculate effect sizes
## ---------------------------------------------------------------------

# MNI in original metric

EAD_ger <- ed_per_item("ger", diffType = "absolute") 
EAD_uk  <- ed_per_item("uk",  diffType = "absolute") 

TEAD <- (p_ger * EAD_ger + p_uk * EAD_uk)

# MNI effect

SEAD_ger <- ed_per_item("ger", diffType = "absolute") / s_pool 
SEAD_uk  <- ed_per_item("uk",  diffType = "absolute") / s_pool

STEAD <- (p_ger * EAD_ger + p_uk * EAD_uk) / s_pool


## ---------------------------------------------------------------------
## 3) Calculate effect directions
## ---------------------------------------------------------------------

ED_ger <- ed_per_item("ger", diffType = "signed")
ED_uk  <- ed_per_item("uk",  diffType = "signed")

MNIdir_ger   <- ED_ger / EAD_ger
MNIdir_uk    <- ED_uk  / EAD_uk
MNIdir_total <- (p_ger * ED_ger + p_uk * ED_uk) /
                (p_ger * EAD_ger + p_uk * EAD_uk)

## ---------------------------------------------------------------------
## 4) Prepare results table
## ---------------------------------------------------------------------
edm_results <- data.frame(
  SEAD_ger      = SEAD_ger,        
  SEAD_uk       = SEAD_uk,         
  STEAD         = STEAD,           
  MNIdir_ger    = MNIdir_ger,     
  MNIdir_uk     = MNIdir_uk,      
  MNIdir_total  = MNIdir_total,   
  row.names     = items
)

round(edm_results, 3)
# round(edm_results, 10)
