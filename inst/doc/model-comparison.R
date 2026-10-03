## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4.5,
  warning  = FALSE,
  message  = FALSE,
  # Evaluate fits during ordinary builds so the source tarball includes
  # the complete rendered guide and its figures.
  eval = TRUE
)


## ----load-and-fit-------------------------------------------------------------
library(terradish)
library(terra)

data(melip)
melip.altitude    <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords      <- terra::unwrap(melip.coords)

covariates <- c(melip.altitude, melip.forestcover)
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(covariates)

surface <- conductance_surface(
  covariates,
  melip.coords,
  directions = 8,
  saveStack = TRUE
)

# Lightweight one-time assessment: reuse recommendations for all fits below.
assessment <- terradish_assess_settings(
  melip.Fst ~ forestcover + altitude,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe,
  probe_maxit       = 1,
  verbose           = FALSE
)
rec <- assessment$recommended

# Full model: both covariates
fit_full <- terradish(
  melip.Fst ~ forestcover + altitude,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

## ----assess-for-comparison, eval = FALSE--------------------------------------
# assessment <- terradish_assess_settings(
#   melip.Fst ~ forestcover + altitude,
#   data              = surface,
#   conductance_model = loglinear_conductance,
#   measurement_model = mlpe,
#   probe_maxit       = 2,
#   verbose           = FALSE
# )
# 
# # One-line summary: did the assessment differ from defaults?
# cat(assessment$comparison$summary, "\n")
# 
# # Extract recommended settings and reuse for every model in this vignette
# rec <- assessment$recommended

## ----lrt----------------------------------------------------------------------
# Reduced model: only forest cover, no altitude
fit_fc_only <- terradish(
  melip.Fst ~ forestcover,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

anova(fit_full, fit_fc_only)

## ----ibd-test-----------------------------------------------------------------
fit_ibd <- terradish(
  melip.Fst ~ 1,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

# Is the full model significantly better than pure IBD?
anova(fit_full, fit_ibd)

## ----interaction-test---------------------------------------------------------
fit_interaction <- terradish(
  melip.Fst ~ forestcover * altitude,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

# Is the interaction term worth keeping?
anova(fit_full, fit_interaction)

## ----fit-all-candidates-------------------------------------------------------
fit_altitude_only <- terradish(
  melip.Fst ~ altitude,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

## ----aic-table----------------------------------------------------------------
aic_table(
  list(fit_ibd, fit_altitude_only, fit_fc_only, fit_full, fit_interaction),
  mod_names = c("IBD", "Altitude only", "Forestcover only",
                "Full (A + FC)", "Interaction (A * FC)")
)

## ----nu-sensitivity, eval = FALSE---------------------------------------------
# # fit_wishart is a converged Wishart fit at the declared primary nu.
# sensitivity <- lapply(c(10, 30, 100), function(value) {
#   terradish_rescale_nu(fit_wishart, nu = value)
# })
# lapply(sensitivity, confint)

## ----aicc-table---------------------------------------------------------------
aic_table(
  list(fit_ibd, fit_altitude_only, fit_fc_only, fit_full, fit_interaction),
  mod_names = c("IBD", "Altitude only", "Forestcover only",
                "Full", "Interaction"),
  AICc = TRUE
)

## ----bic-table----------------------------------------------------------------
aic_table(
  list(fit_ibd, fit_altitude_only, fit_fc_only, fit_full, fit_interaction),
  mod_names = c("IBD", "Altitude only", "Forestcover only",
                "Full", "Interaction"),
  BIC = TRUE
)

## ----likelihood-grid----------------------------------------------------------
# Define a 21 x 21 grid of (forestcover, altitude) combinations
theta_grid <- as.matrix(
  expand.grid(
    forestcover = seq(-1, 1, length.out = 21),
    altitude    = seq(-1, 1, length.out = 21)
  )
)

grid_result <- terradish_grid(
  theta     = theta_grid,
  formula   = melip.Fst ~ forestcover + altitude,
  data      = surface,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe
)

## ----plot-likelihood-surface, fig.cap = "***Log-likelihood surface. Red dot marks the MLE.***"----
# Build a data frame for plotting
grid_df <- data.frame(
  forestcover = theta_grid[, "forestcover"],
  altitude    = theta_grid[, "altitude"],
  loglik      = grid_result$loglik
)

# Base-R image plot (no ggplot2 dependency)
loglik_mat <- matrix(
  grid_df$loglik,
  nrow = length(unique(grid_df$forestcover)),
  ncol = length(unique(grid_df$altitude))
)

image(
  x    = sort(unique(grid_df$forestcover)),
  y    = sort(unique(grid_df$altitude)),
  z    = loglik_mat,
  xlab = expression(theta[forestcover]),
  ylab = expression(theta[altitude]),
  main = "Log-likelihood surface",
  col  = hcl.colors(64, "Blues", rev = TRUE)
)
contour(
  x    = sort(unique(grid_df$forestcover)),
  y    = sort(unique(grid_df$altitude)),
  z    = loglik_mat,
  add  = TRUE
)
points(coef(fit_full)["forestcover"],
       coef(fit_full)["altitude"],
       pch = 19, col = "red", cex = 1.5)
legend("topright", legend = "MLE", pch = 19, col = "red", bty = "n")

## ----likelihood-grid-coarse, eval=FALSE---------------------------------------
# grid_coarse <- terradish_grid(
#   theta = theta_grid,
#   formula = melip.Fst ~ forestcover + altitude,
#   data = surface,
#   conductance_model = loglinear_conductance,
#   measurement_model = mlpe,
#   approximation = "coarse_raster",
#   approximation_control = list(factor = 2L)
# )

## ----ci-conductance, fig.cap = "***Fitted conductance with 95% confidence bounds.***"----
plot(fit_full, type = "surface", data = surface)

## ----marginal-conductance, fig.cap = "***Marginal association of each covariate with conductance.***"----
plot(fit_full, type = "marginal", data = surface)

## ----support-clamp-marginal, eval = FALSE-------------------------------------
# plot(
#   fit_full,
#   type = "marginal",
#   data = surface,
#   support = "focal",
#   support_probs = c(0.01, 0.99),
#   clamp_covariates = c("altitude", "forestcover")
# )

## ----support-clamp-conductance, eval = FALSE----------------------------------
# cond_full_focal <- conductance(
#   surface,
#   fit_full,
#   support = "focal",
#   support_probs = c(0.01, 0.99),
#   clamp_covariates = c("altitude", "forestcover")
# )

## ----cv-spatial-folds, eval = FALSE-------------------------------------------
# # melip.coords is geographic, so project it before measuring planar separation.
# # Center a local distance projection on the sampled Brazilian sites.
# lonlat <- terra::crds(terra::project(melip.coords, "EPSG:4326"))
# local_crs <- sprintf("+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m",
#                      mean(lonlat[, 2]), mean(lonlat[, 1]))
# coords_projected <- terra::project(melip.coords, local_crs)
# folds_spatial <- terradish_folds(
#   coords_projected,
#   k = 5,
#   method = "spatial_kmeans",
#   seed = 42
# )
# 
# cv_spatial <- terradish_cv_folds(
#   data = surface,
#   formulas = list(
#     "Forest cover" = melip.Fst ~ forestcover,
#     "Forest cover + altitude" = melip.Fst ~ forestcover + altitude
#   ),
#   folds = folds_spatial,
#   model = mlpe,
#   nuisance = "fixed",
#   keep_fits = "slim",
#   checkpoint = file.path(tempdir(), "terradish-spatial-cv.rds")
# )
# 
# cv_spatial$summary
# cv_spatial$results

## ----cv-repeats, eval = FALSE-------------------------------------------------
# folds_repeated <- lapply(c(42, 91, 137), function(seed) {
#   terradish_folds(coords_projected, k = 5, seed = seed)
# })
# cv_shape <- terradish_cv_folds(
#   data = surface,
#   formulas = list(linear = melip.Fst ~ forestcover,
#                   spline = melip.Fst ~ s(forestcover, df = 3)),
#   folds = folds_repeated,
#   model = mlpe,
#   conductance_model = list(linear = loglinear_conductance,
#                            spline = smooth_loglinear_conductance),
#   nuisance = "fixed"
# )
# cv_shape$summary
# cv_shape$results

## ----qref-mc, eval = FALSE----------------------------------------------------
# library(terradish)
# 
# # 1. Build the surface once and reuse it for every candidate model
# surface <- conductance_surface(covariates, coords, directions = 8, saveStack = TRUE)
# 
# # 2. Profile the graph once; reuse the recommended settings for every fit so
# #    models differ only in their formula, never in their optimizer
# rec <- terradish_assess_settings(melip.Fst ~ forestcover + altitude,
#                                  data = surface)$recommended
# 
# fit_full <- terradish(melip.Fst ~ forestcover + altitude, data = surface,
#                       conductance_model = loglinear_conductance,
#                       measurement_model = mlpe,
#                       optimizer = rec$optimizer, control = rec$control)
# fit_alt  <- terradish(melip.Fst ~ altitude, data = surface,
#                       conductance_model = loglinear_conductance,
#                       measurement_model = mlpe,
#                       optimizer = rec$optimizer, control = rec$control)
# 
# # 3. Nested hypothesis test: does forest cover add likelihood support?
# anova(fit_alt, fit_full)          # interpret with nesting and boundary checks
# 
# # 4. Rank the whole set; AICc when sites are few relative to parameters
# aic_table(
#   list(fit_full, fit_alt),
#   mod_names = c("Full", "Altitude"),
#   AICc = TRUE
# )                                 # small delta is a heuristic within this set
# 
# # 5. Fixed-domain transfer across spatial blocks. Prefer design-based blocks;
# #    otherwise cluster projected coordinates without standardizing them.
# # Center a local distance projection on the sampled Brazilian sites.
# lonlat <- terra::crds(terra::project(coords, "EPSG:4326"))
# local_crs <- sprintf("+proj=aeqd +lat_0=%f +lon_0=%f +datum=WGS84 +units=m",
#                      mean(lonlat[, 2]), mean(lonlat[, 1]))
# coords_projected <- terra::project(coords, local_crs)
# folds <- terradish_folds(coords_projected, k = 5, seed = 42)
# cv <- terradish_cv_folds(
#   data = surface,
#   formulas = list(full = melip.Fst ~ forestcover + altitude,
#                   altitude = melip.Fst ~ altitude),
#   folds = folds,
#   model = mlpe, nuisance = "fixed"
# )
# cv$summary                        # paired differences on common successful folds
# 
# # 6. Map the surface with uncertainty
# plot(fit_full, type = "surface", data = surface)

