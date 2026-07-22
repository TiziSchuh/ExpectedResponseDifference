## =====================================================================
## Step 1: Fetching parameters
## =====================================================================

# source("0_fit_models.R") # run if step0-script has not been executed yet

## ---------------------------------------------------------------------
## 1) Group-specific parameters from the configural model
##    tau    = item intercepts        (lavaan matrix `nu`)
##    Lambda = factor loadings        (lavaan matrix `lambda`)
##    kappa  = latent means           (lavaan matrix `alpha`)
##    Sigma  = latent (co)variances   (lavaan matrix `psi`)
## ---------------------------------------------------------------------
est_mg    <- lavInspect(fa.model.mg, what = "est")
grp_labels <- as.character(lavInspect(fa.model.mg, what = "group.label"))

i_ger <- which(grp_labels == "1")   # ger: COUNTRY == 1
i_uk  <- which(grp_labels == "2")   # uk:  COUNTRY == 2

tau_ger    <- est_mg[[i_ger]]$nu
Lambda_ger <- est_mg[[i_ger]]$lambda
kappa_ger  <- est_mg[[i_ger]]$alpha
Sigma_ger  <- est_mg[[i_ger]]$psi

tau_uk     <- est_mg[[i_uk]]$nu
Lambda_uk  <- est_mg[[i_uk]]$lambda
kappa_uk   <- est_mg[[i_uk]]$alpha
Sigma_uk   <- est_mg[[i_uk]]$psi

## NOTE: under lavaan's default configural identification the latent
## means are fixed to 0 in BOTH groups, so kappa_ger and kappa_uk will
## be zero vectors here. Sigma_ger / Sigma_uk are freely estimated
## (marker-variable identification, first loading per factor = 1;
## for `aqui` all loadings are fixed to +/-1, so its variance is free).

## ---------------------------------------------------------------------
## 2) Single-group (ger + uk pooled) parameters
## ---------------------------------------------------------------------
est_sg    <- lavInspect(fa.model.sg, what = "est")
tau_sg    <- est_sg$nu
Lambda_sg <- est_sg$lambda

## ---------------------------------------------------------------------
## 3) Pooled SD per item (observed item scores, pooled across countries)
##    Classic pooled-within-group formula:
##    s_pool = sqrt( ((n1-1)*s1^2 + (n2-1)*s2^2) / (n1 + n2 - 2) )
##    with item-specific complete-case n's per group (data contain NAs).
## ---------------------------------------------------------------------
pooled_sd <- function(x, g) {
  x1 <- x[g == 1]
  x2 <- x[g == 2]
  n1 <- sum(!is.na(x1))
  n2 <- sum(!is.na(x2))
  sqrt(((n1 - 1) * var(x1, na.rm = TRUE) +
        (n2 - 1) * var(x2, na.rm = TRUE)) / (n1 + n2 - 2))
}

items <- lavNames(fa.model.mg, type = "ov") # Observed items in the model

s_pool <- vapply(via[items], pooled_sd, numeric(1), g = via$COUNTRY)
s_pool

## ---------------------------------------------------------------------
## 4) Group proportions (sample-based)
##    Taken from the n's lavaan actually used in the configural fit
##    (rows where all model variables are missing are dropped by FIML).
## ---------------------------------------------------------------------
n_grp <- lavInspect(fa.model.mg, what = "nobs")   # ordered like grp_labels
n_ger <- n_grp[i_ger]
n_uk  <- n_grp[i_uk]

p_ger <- n_ger / (n_ger + n_uk)
p_uk  <- n_uk  / (n_ger + n_uk)

## ---------------------------------------------------------------------
## Collect everything in one object
## ---------------------------------------------------------------------
step1 <- list(
  ger    = list(tau = tau_ger, Lambda = Lambda_ger,
                kappa = kappa_ger, Sigma = Sigma_ger),
  uk     = list(tau = tau_uk, Lambda = Lambda_uk,
                kappa = kappa_uk, Sigma = Sigma_uk),
  single = list(tau = tau_sg, Lambda = Lambda_sg),
  s_pool = s_pool,
  p      = c(ger = unname(p_ger), uk = unname(p_uk)),
  n      = c(ger = unname(n_ger), uk = unname(n_uk))
)

print(step1)
