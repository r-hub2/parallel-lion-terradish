## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4.5,
  warning  = FALSE,
  message  = FALSE
)


## ----eval = FALSE-------------------------------------------------------------
# # Released version
# install.packages("terradish")
# 
# # Development version
# remotes::install_github("wpeterman/terradish")

## ----eval = FALSE-------------------------------------------------------------
# install.packages("adegenet")

## ----load---------------------------------------------------------------------
library(terradish)
library(terra)

## ----data---------------------------------------------------------------------
data(melip)

# The rasters are stored as "packed" objects to comply with CRAN policy;
# unwrap() converts them back to usable SpatRaster / SpatVector objects.
melip.altitude    <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords      <- terra::unwrap(melip.coords)

## ----inspect-data-------------------------------------------------------------
# How many sampling sites?
nrow(melip.coords)

# Quick look at the genetic distance matrix
# (only the first 5 x 5 block for readability)
round(melip.Fst[1:5, 1:5], 3)

## ----plot-covariates, fig.cap = "***Altitude raster with sampling sites.***"----
plot(melip.altitude, main = "Altitude (rescaled)")
points(melip.coords, pch = 19, cex = 0.8, col = 'red')

## ----plot-forestcover, fig.cap = "***Forest cover raster (0 = no forest, 1 = full cover).***"----
plot(melip.forestcover, main = "Forest cover")
points(melip.coords, pch = 19, cex = 0.8, col = 'red')

## ----scale--------------------------------------------------------------------
covariates <- c(melip.altitude, melip.forestcover)
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(covariates)

## ----surface------------------------------------------------------------------
surface <- conductance_surface(covariates, melip.coords, directions = 8)

## ----surface-check------------------------------------------------------------
# How many grid cells are in the graph?
nrow(surface$x)

# The covariate data for the first few cells
head(surface$x)

## ----crop-buffer-example, eval=FALSE------------------------------------------
# surface_cropped <- conductance_surface(
#   covariates,
#   melip.coords,
#   directions = 8,
#   saveStack = TRUE,
#   crop_buffer = 0.05
# )

## ----coarse-raster-example, eval=FALSE----------------------------------------
# surface_for_coarse <- conductance_surface(
#   covariates,
#   melip.coords,
#   directions = 8,
#   saveStack = TRUE
# )
# 
# fit_coarse <- terradish(
#   melip.Fst ~ forestcover + altitude,
#   data = surface_for_coarse,
#   conductance_model = loglinear_conductance,
#   measurement_model = mlpe,
#   approximation = "coarse_raster",
#   approximation_control = list(
#     factor = c(4, 2),
#     exact_refine = TRUE
#   )
# )

## ----assess-settings-demo, eval = FALSE---------------------------------------
# # saveStack = TRUE is required for the coarse-raster probe
# surface_for_assess <- conductance_surface(
#   covariates, melip.coords,
#   directions = 8, saveStack = TRUE
# )
# 
# assessment <- terradish_assess_settings(
#   melip.Fst ~ forestcover + altitude,
#   data              = surface_for_assess,
#   conductance_model = loglinear_conductance,
#   measurement_model = mlpe,
#   probe_maxit       = 2,      # keep small: probes compare early-iteration cost
#   coarse_probe      = TRUE,   # also time the coarse-raster warm-start path
#   verbose           = TRUE    # print progress as probes run
# )
# 
# # Print a full summary: graph profile, probe results, recommended settings
# assessment

## ----assess-use, eval = FALSE-------------------------------------------------
# fit_tuned <- terradish(
#   melip.Fst ~ forestcover + altitude,
#   data              = surface_for_assess,
#   conductance_model = loglinear_conductance,
#   measurement_model = mlpe,
#   optimizer         = assessment$recommended$optimizer,
#   control           = assessment$recommended$control,
#   solver            = assessment$recommended$solver,
#   solver_control    = assessment$recommended$solver_control,
#   approximation     = assessment$recommended$approximation,
#   approximation_control = assessment$recommended$approximation_control
# )

## ----assess-edit, eval = FALSE------------------------------------------------
# my_control <- assessment$recommended$control
# my_control$maxit <- 50          # allow more iterations than the probe default
# 
# fit_tuned2 <- terradish(
#   melip.Fst ~ forestcover + altitude,
#   data    = surface_for_assess,
#   conductance_model = loglinear_conductance,
#   measurement_model = mlpe,
#   optimizer      = assessment$recommended$optimizer,
#   control        = my_control,
#   solver         = assessment$recommended$solver,
#   solver_control = assessment$recommended$solver_control
# )

## ----fit-leastsquares---------------------------------------------------------
fit_ls <- terradish(
  melip.Fst ~ forestcover + altitude,  # formula: response ~ covariates
  data              = surface,          # the conductance surface
  conductance_model = loglinear_conductance,
  measurement_model = leastsquares      # ordinary regression
)

## ----summary-ls---------------------------------------------------------------
summary(fit_ls)

## ----fit-mlpe-----------------------------------------------------------------
fit_mlpe <- terradish(
  melip.Fst ~ forestcover + altitude,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe
)

summary(fit_mlpe)

## ----fit-wishart--------------------------------------------------------------
# Construct an admissible squared-distance response from a simulated positive
# semidefinite covariance matrix. The same nu is used for simulation and fit.
nu_gw <- 50
sim_gw <- simulate_covariance_response(
  theta   = c(forestcover = 0.4, altitude = -0.3),
  formula = ~ forestcover + altitude,
  data    = surface,
  tau     = 0.8,
  sigma   = 0.2,
  nu      = nu_gw,
  seed    = 1
)
gw_distance <- dist_from_cov(sim_gw$covariance)
gw_check <- check_distance_response(gw_distance)
gw_check[c("admissible", "relative_min_eigenvalue", "n_negative")]

fit_gw <- terradish(
  gw_distance ~ forestcover + altitude,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = generalized_wishart,
  nu                = nu_gw
)

summary(fit_gw)

## ----compare-measurement-models-----------------------------------------------
# Valid: both use a Gaussian likelihood on pairwise distances.
aic_table(
  list(fit_ls, fit_mlpe),
  mod_names = c("leastsquares", "mlpe")
)

## ----coef---------------------------------------------------------------------
coef(fit_mlpe)

## ----interpret----------------------------------------------------------------
theta <- coef(fit_mlpe)
cat("One SD increase in forest cover multiplies conductance by",
    round(exp(theta["forestcover"]), 2), "\n")
cat("One SD increase in altitude multiplies conductance by",
    round(exp(theta["altitude"]), 2), "\n")

## ----plot-fit, fig.cap = "***Observed vs. resistance-distance-fitted genetic distances.***"----
plot(fit_mlpe, type = "fit")

## ----conductance-surface, fig.cap = "***Estimated conductance surface with 95% confidence interval bounds.***"----
plot(fit_mlpe, type = "surface", data = surface)

## ----marginal-effects, fig.cap = "***Marginal association of each covariate with predicted genetic distance, with predictive bands.***"----
plot(fit_mlpe, data = surface)

## ----marginal-conductance, fig.cap = "***Marginal association of each covariate with conductance, with a 95% CI.***"----
plot(fit_mlpe, type = "marginal", data = surface)

## ----support-clamp, eval = FALSE----------------------------------------------
# # Clamp selected covariates to focal-site quantile support before plotting.
# plot(
#   fit_mlpe,
#   type = "marginal",
#   data = surface,
#   support = "focal",
#   support_probs = c(0.01, 0.99),
#   clamp_covariates = c("altitude", "forestcover")
# )
# 
# # Apply the same support rule to raster conductance prediction.
# cond_focal <- conductance(
#   surface,
#   fit_mlpe,
#   support = "focal",
#   support_probs = c(0.01, 0.99),
#   clamp_covariates = c("altitude", "forestcover")
# )

## ----full-workflow, eval = FALSE----------------------------------------------
# library(terradish)
# library(terra)
# 
# # 1. Load data
# data(melip)
# melip.altitude    <- terra::unwrap(melip.altitude)
# melip.forestcover <- terra::unwrap(melip.forestcover)
# melip.coords      <- terra::unwrap(melip.coords)
# 
# # 2. Scale covariates and build surface
# covariates <- c(melip.forestcover, melip.altitude)
# names(covariates) <- c("forestcover", "altitude")
# covariates <- scale_covariates(covariates)
# 
# surface <- conductance_surface(covariates, melip.coords,
#                                directions = 8, saveStack = TRUE)
# 
# # 2b. (Optional) Profile graph and identify fastest settings
# #     Most useful for large graphs (> 10 000 vertices) or spline models.
# assessment <- terradish_assess_settings(
#   melip.Fst ~ forestcover + altitude,
#   data              = surface,
#   conductance_model = loglinear_conductance,
#   measurement_model = mlpe,
#   probe_maxit       = 2
# )
# cat(assessment$comparison$summary, "\n")  # one-line verdict
# rec <- assessment$recommended             # reuse for all fits below
# 
# # 3. Fit model - choose a measurement model:
# #   leastsquares      : fast, independent errors (distance matrix)
# #   mlpe              : shared-site correlation correction (distance matrix)
# #   generalized_wishart: Wishart model for an admissible squared-distance response;
# #                        requires nu
# #   wishart_covariance : Wishart model, requires nu (covariance matrix)
# #   For Wishart IBE or site-level environmental covariance kernels, use
# #   wishart_covariates(); pair_subset_measurement_model() is for pairwise
# #   regression rows and does not apply to the full-matrix Wishart likelihood.
# fit <- terradish(melip.Fst ~ forestcover + altitude,
#                  data              = surface,
#                  conductance_model = loglinear_conductance,
#                  measurement_model = mlpe,
#                  optimizer         = rec$optimizer,
#                  control           = rec$control,
#                  solver            = rec$solver,
#                  solver_control    = rec$solver_control,
#                  approximation     = rec$approximation,
#                  approximation_control = rec$approximation_control)
# 
# # For generalized_wishart, derive and check a squared-distance response:
# # S_cov <- cov_from_genetic_data(genotype_matrix, groups = population)
# # S_dist <- dist_from_cov(S_cov)
# # check_distance_response(S_dist)
# # fit_gw <- terradish(S_dist ~ forestcover + altitude, data = surface,
# #                     conductance_model = loglinear_conductance,
# #                     measurement_model = generalized_wishart,
# #                     nu = effective_marker_df)
# 
# # 4. Summarize and visualize
# summary(fit)
# 
# plot(fit, type = "fit")                                # observed vs. fitted scatter
# plot(fit, type = "surface", data = surface)            # conductance map + 95% CI panels
# plot(fit, data = surface)                              # marginal associations on genetic distance (default)
# plot(fit, type = "marginal", data = surface)           # marginal associations on conductance scale
# 
# # 5. Optional: focal-support clamping for robust tails
# plot(fit, type = "marginal", data = surface,
#      support = "focal", support_probs = c(0.01, 0.99),
#      clamp_covariates = c("forestcover", "altitude"))
# 
# cond <- conductance(surface, fit,
#                     support = "focal", support_probs = c(0.01, 0.99),
#                     clamp_covariates = c("forestcover", "altitude"))

