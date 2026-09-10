## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4.5,
  warning  = FALSE,
  message  = FALSE
)


## ----load---------------------------------------------------------------------
library(terradish)
library(terra)

data(melip)
melip.altitude    <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords      <- terra::unwrap(melip.coords)

# Scale continuous covariates before building the conductance surface.
# This keeps the conductance parameters on a comparable scale and avoids
# numerical overflow in exp(theta * x).
covariates <- c(melip.altitude, melip.forestcover)
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(covariates)

surface <- conductance_surface(covariates, melip.coords, directions = 8)

## ----ibe-covariates-----------------------------------------------------------
ibe_altitude <- melip.altitude
names(ibe_altitude) <- "altitude_site"

z_ibe <- pairwise_endpoint_covariates(
  ibe_altitude,            # site-level environmental data (raster or matrix)
  coords    = melip.coords,# focal-point locations (used to extract from raster)
  transform = "absdiff",  # how to convert site values into pairwise values
  scale     = TRUE         # standardize site values first
)

# Each row is one pair of sites; each column is one pairwise covariate.
dim(z_ibe)
head(z_ibe)

## ----ibe-covariates-multi-----------------------------------------------------
# Example: separate absolute differences for altitude and forest cover
z_multi <- pairwise_endpoint_covariates(
  covariates,            # the two-layer scaled raster stack
  coords    = melip.coords,
  transform = "absdiff",
  scale     = FALSE      # already scaled in the raster
)
dim(z_multi)   # one column per layer
colnames(z_multi)

## ----build-model--------------------------------------------------------------
g_ibe <- mlpe_covariates(z_ibe)

## ----fit-models---------------------------------------------------------------
# NewtonRaphsonControl limits the number of optimizer iterations.
# The vignette caps maxit at 20. In an applied analysis, allow convergence and
# verify optimizer diagnostics rather than relying on this illustrative cap.
ctrl <- NewtonRaphsonControl(maxit = 20, verbose = FALSE)

# IBD: null model, no landscape or environment predictors
fit_ibd <- terradish(
  melip.Fst ~ 1,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe,
  leverage          = FALSE,
  control           = ctrl
)

# IBE only: environmental mismatch, but no landscape resistance
fit_ibe <- terradish(
  melip.Fst ~ 1,           # no raster covariates on the conductance side
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = g_ibe, # IBE is on the measurement side
  leverage          = FALSE,
  control           = ctrl
)

# IBR only: landscape resistance, no environmental mismatch
fit_ibr <- terradish(
  melip.Fst ~ altitude + forestcover,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe,
  leverage          = FALSE,
  control           = ctrl
)

# Joint IBE + IBR: both processes in a single model
fit_joint <- terradish(
  melip.Fst ~ altitude + forestcover,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = g_ibe,
  leverage          = FALSE,
  control           = ctrl
)

## ----coef-ibr-----------------------------------------------------------------
coef(fit_ibr)   # IBR-only model
coef(fit_joint) # joint model (theta side only)

## ----summary-joint------------------------------------------------------------
summary(fit_joint)

## ----gamma-ci-----------------------------------------------------------------
summary(fit_joint)$phi_table["absdiff_altitude_site", ]

## ----response-scale-change----------------------------------------------------
mlpe_response_change(fit_joint)

## ----response-scale-change-custom---------------------------------------------
mlpe_response_change(fit_joint,
                     covariate = "absdiff_altitude_site",
                     values = c(0, 2))

## ----aic-table----------------------------------------------------------------
aic_table(
  list(fit_ibd, fit_ibe, fit_ibr, fit_joint),
  mod_names = c("IBD", "IBE only", "IBR only", "IBE + IBR (joint)")
)

## ----aicc-table---------------------------------------------------------------
aic_table(
  list(fit_ibd, fit_ibe, fit_ibr, fit_joint),
  mod_names = c("IBD", "IBE only", "IBR only", "IBE + IBR (joint)"),
  AICc = TRUE
)

## ----lrt----------------------------------------------------------------------
# Does the joint model improve on IBR-only?
anova(fit_joint, fit_ibr)

# Does IBR-only improve on IBD?
anova(fit_ibr, fit_ibd)

## ----cv-model-selection-------------------------------------------------------
cv_ctrl <- NewtonRaphsonControl(maxit = 15, verbose = FALSE)

# IBR-only: has conductance parameters, so held-out evaluation is supported
cv_ibr <- terradish_cv(
  pts = melip.coords, covariates = covariates,
  fmla = melip.Fst ~ altitude + forestcover, model = mlpe,
  prop_train = 2 / 3, seed = 1, fit_full = TRUE,
  control = cv_ctrl
)

# Joint IBE + IBR: uses the same held-out split
cv_joint <- terradish_cv(
  pts = melip.coords, covariates = covariates,
  fmla = melip.Fst ~ altitude + forestcover, model = g_ibe,
  prop_train = 2 / 3, seed = 1, fit_full = TRUE,
  control = cv_ctrl
)

# Compare only models evaluated on held-out data from the same split.
comparison <- cv_model_selection(
  list(
    list(train_mod = cv_ibr$train_mod,
         cv_loglik  = cv_ibr$cv_loglik,
         full_mod   = cv_ibr$full_mod),
    list(train_mod = cv_joint$train_mod,
         cv_loglik  = cv_joint$cv_loglik,
         full_mod   = cv_joint$full_mod)
  ),
  cv_names = c("IBR only", "IBE + IBR (joint)"),
  aic = TRUE
)

# Held-out log-likelihood (higher is better predictive performance)
comparison$loglik_tab

# AIC table from the full-data fits
comparison$AIC_tab

## ----cv-replicates------------------------------------------------------------
cv_joint_reps <- terradish_cv_replicates(
  pts        = melip.coords,
  covariates = covariates,
  fmla       = melip.Fst ~ altitude + forestcover,
  model      = g_ibe,
  seeds      = 1:5,
  fit_full   = FALSE,
  keep_fits  = FALSE,
  control    = cv_ctrl
)

summary(cv_joint_reps)

## ----obs-vs-fitted, fig.cap = "***Observed vs. fitted pairwise genetic distances from the joint model.***"----
plot(fit_joint, type = "fit")

## ----conductance-surface, fig.cap = "***Estimated conductance surface from the joint IBE + IBR model.***"----
plot(fit_joint, type = "surface", data = surface)

## ----marginal-joint, fig.cap = "***Marginal conductance associations for the joint IBE + IBR model.***"----
plot(fit_joint, type = "marginal", data = surface)

## ----marginal-ibr, fig.cap = "***Marginal conductance associations for the IBR-only model.***"----
plot(fit_ibr, type = "marginal", data = surface)

## ----ibe-ibr-support-clamp, eval = FALSE--------------------------------------
# plot(
#   fit_joint,
#   type = "marginal",
#   data = surface,
#   support = "focal",
#   support_probs = c(0.01, 0.99),
#   clamp_covariates = c("altitude", "forestcover")
# )
# 
# cond_joint_focal <- conductance(
#   surface,
#   fit_joint,
#   support = "focal",
#   support_probs = c(0.01, 0.99),
#   clamp_covariates = c("altitude", "forestcover")
# )

## ----saved-results------------------------------------------------------------
# Create a temporary directory to illustrate the pattern.
tmp_dir <- tempfile("terradish-results-")
dir.create(tmp_dir)

# Save the joint fit and a minimal AllResults object.
saveRDS(fit_joint, file.path(tmp_dir, "fit--mlpe.rds"))
saveRDS(
  list(
    effect_size = coef(fit_joint),
    pts         = melip.coords,
    covariates  = covariates
  ),
  file.path(tmp_dir, "AllResults_list.rds")
)

results_info <- terradish_results(tmp_dir)
names(results_info)

## ----saved-parameters---------------------------------------------------------
parameter_table <- terradish_parameters(
  tmp_dir,
  model      = "mlpe",
  save_table = FALSE
)

parameter_table

## ----wishart-ibe-setup--------------------------------------------------------
nu_demo <- 50
sim_wishart <- simulate_covariance_response(
  theta   = c(altitude = 0.35, forestcover = -0.25),
  formula = ~ altitude + forestcover,
  data    = surface,
  tau     = 0.8,
  sigma   = 0.2,
  nu      = nu_demo,
  seed    = 2026
)
S_wishart <- dist_from_cov(sim_wishart$covariance)
check_distance_response(S_wishart)

# Altitude kernel for the generalized Wishart likelihood
g_wishart_ibe <- wishart_covariates(
  melip.altitude,
  coords = melip.coords,
  model  = "generalized_wishart",
  scale  = TRUE        # standardize site values before forming the kernel
)

# phi starting values reveal the parameter names that will be estimated
g_wishart_ibe(diag(nrow(melip.coords)),
              matrix(0, nrow(melip.coords), nrow(melip.coords)))$phi

## ----wishart-fit-models-------------------------------------------------------
# Wishart IBD: null conductance, no environmental kernels
fit_w_ibd <- terradish(
  S_wishart ~ 1,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = generalized_wishart,
  nu                = nu_demo,
  leverage          = FALSE,
  control           = ctrl
)

# Wishart IBE only: altitude kernel, uniform conductance
fit_w_ibe <- terradish(
  S_wishart ~ 1,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = g_wishart_ibe,
  nu                = nu_demo,
  leverage          = FALSE,
  control           = ctrl
)

# Wishart IBR only: landscape covariates, no environmental kernel
fit_w_ibr <- terradish(
  S_wishart ~ altitude + forestcover,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = generalized_wishart,
  nu                = nu_demo,
  leverage          = FALSE,
  control           = ctrl
)

# Wishart joint: landscape covariates + altitude IBE kernel
fit_w_joint <- terradish(
  S_wishart ~ altitude + forestcover,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = g_wishart_ibe,
  nu                = nu_demo,
  leverage          = FALSE,
  control           = ctrl
)

## ----wishart-summary----------------------------------------------------------
summary(fit_w_joint)

## ----wishart-aic--------------------------------------------------------------
aic_table(
  list(fit_w_ibd, fit_w_ibe, fit_w_ibr, fit_w_joint),
  mod_names = c("Wishart IBD", "Wishart IBE", "Wishart IBR",
                "Wishart IBE + IBR"),
  AICc = TRUE
)

## ----wishart-lrt--------------------------------------------------------------
# Does adding the altitude IBE kernel improve on IBR alone?
anova(fit_w_joint, fit_w_ibr)

## ----quick-reference, eval = FALSE--------------------------------------------
# library(terradish)
# library(terra)
# 
# # 1. Load and scale data
# data(melip)
# melip.altitude    <- terra::unwrap(melip.altitude)
# melip.forestcover <- terra::unwrap(melip.forestcover)
# melip.coords      <- terra::unwrap(melip.coords)
# 
# covariates <- c(melip.altitude, melip.forestcover)
# names(covariates) <- c("altitude", "forestcover")
# covariates <- scale_covariates(covariates)
# surface <- conductance_surface(covariates, melip.coords, directions = 8)
# 
# # 2a. MLPE approach: build IBE covariates and measurement model
# z_ibe <- pairwise_endpoint_covariates(
#   melip.altitude, coords = melip.coords,
#   transform = "absdiff", scale = TRUE
# )
# g_ibe <- mlpe_covariates(z_ibe)
# 
# # 3. Fit IBD, IBE, IBR, and joint models (MLPE)
# fit_ibd   <- terradish(melip.Fst ~ 1,
#                        surface, loglinear_conductance, mlpe)
# fit_ibe   <- terradish(melip.Fst ~ 1,
#                        surface, loglinear_conductance, g_ibe)
# fit_ibr   <- terradish(melip.Fst ~ altitude + forestcover,
#                        surface, loglinear_conductance, mlpe)
# fit_joint <- terradish(melip.Fst ~ altitude + forestcover,
#                        surface, loglinear_conductance, g_ibe)
# 
# # 4. Compare
# aic_table(list(fit_ibd, fit_ibe, fit_ibr, fit_joint),
#           mod_names = c("IBD", "IBE", "IBR", "IBE + IBR"))
# anova(fit_joint, fit_ibr)   # does IBE add to IBR?
# anova(fit_ibr,   fit_ibd)   # does IBR improve on IBD?
# 
# # 5. Interpret theta (IBR) and gamma (IBE inside phi)
# coef(fit_joint)
# summary(fit_joint)$phi
# mlpe_response_change(fit_joint)
# 
# # 6. Visualize
# plot(fit_joint, type = "fit")                       # observed vs. fitted
# plot(fit_joint, type = "surface", data = surface)   # conductance map + CI panels
# plot(fit_joint, type = "marginal", data = surface)  # marginal covariate associations
# 
# # 2b. Wishart alternative: IBE as a covariance kernel (requires effective nu)
# # First construct an admissible squared-distance response from a positive semidefinite
# # covariance matrix; do not substitute an unchecked dissimilarity matrix.
# # S_cov <- cov_from_genetic_data(genotype_matrix, groups = population)
# # S_dist <- dist_from_cov(S_cov)
# g_wishart_ibe <- wishart_covariates(
#   melip.altitude, coords = melip.coords,
#   model = "generalized_wishart", scale = TRUE
# )
# effective_marker_df <- 1000  # replace with a justified value and sensitivity range
# fit_w_joint <- terradish(
#   S_dist ~ altitude + forestcover,
#   data = surface, conductance_model = loglinear_conductance,
#   measurement_model = g_wishart_ibe, nu = effective_marker_df
# )
# summary(fit_w_joint)  # phi: tau (IBR), lambda_altitude (IBE), sigma (nugget)

