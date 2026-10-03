## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4.25,
  warning  = FALSE,
  message  = FALSE
)

## ----load---------------------------------------------------------------------
library(terradish)

data(melip)
melip.altitude    <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords      <- terra::unwrap(melip.coords)

covariates <- c(melip.altitude, melip.forestcover)
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(covariates)
surface    <- conductance_surface(covariates, melip.coords, directions = 8)

## ----solver-calls, eval = FALSE-----------------------------------------------
# # direct (the default): exact sparse Cholesky factorization
# fit <- terradish(melip.Fst ~ forestcover + altitude, surface,
#                  conductance_model = loglinear_conductance,
#                  solver = "direct")
# 
# # amg: algebraic multigrid, the large-N path
# fit <- terradish(melip.Fst ~ forestcover + altitude, surface,
#                  conductance_model = loglinear_conductance,
#                  solver = "amg")

## ----measurement-control, eval = FALSE----------------------------------------
# fit_smooth_gw <- terradish(
#   admissible_squared_distance ~ s(forestcover, df = 4) +
#     s(altitude, df = 4),
#   surface,
#   conductance_model = smooth_loglinear_conductance,
#   measurement_model = generalized_wishart,
#   nu = effective_marker_df,
#   optimizer = "bfgs",
#   control = NewtonRaphsonControl(maxit = 80, verbose = TRUE),
#   measurement_control = NewtonRaphsonControl(
#     maxit = 40,
#     ctol = 1e-8,
#     ftol = 1e-10,
#     verbose = FALSE
#   )
# )
# 
# fit_smooth_gw$fit$subproblem

## ----curvature-calls, eval = FALSE--------------------------------------------
# # Generalized Wishart requires an admissible squared-distance response that is
# # coherently related to a centered positive-semidefinite covariance matrix.
# # exact (default): full observed Hessian
# fit <- terradish(admissible_squared_distance ~ forestcover + altitude, surface,
#                  measurement_model = generalized_wishart,
#                  nu = effective_marker_df,
#                  curvature = "exact")
# 
# # gauss_newton: information-based curvature
# fit <- terradish(admissible_squared_distance ~ forestcover + altitude, surface,
#                  measurement_model = generalized_wishart,
#                  nu = effective_marker_df,
#                  curvature = "gauss_newton")

## ----curvature-benchmark-run, eval = FALSE------------------------------------
# source(system.file("benchmarks", "benchmark-curvature-large-landscape.R",
#                    package = "terradish"))
# 
# bench <- run_curvature_large_landscape_benchmark(CONFIG)
# bench$comparison

## ----crop-sensitivity, eval = FALSE-------------------------------------------
# # Scale on the original domain once, before cropping.
# buffers <- c(2, 5, 10) * max(terra::res(covariates))
# fits_buffer <- lapply(buffers, function(buffer) {
#   cropped_surface <- conductance_surface(
#     covariates, melip.coords, directions = 8, crop_buffer = buffer
#   )
#   terradish(melip.Fst ~ forestcover + altitude, cropped_surface,
#             measurement_model = mlpe)
# })
# lapply(fits_buffer, coef)
# lapply(fits_buffer, confint)

## ----coarse-warm-start, eval = FALSE------------------------------------------
# fit <- terradish(
#   melip.Fst ~ forestcover + altitude, surface,
#   measurement_model = mlpe,
#   approximation = "coarse_raster",
#   approximation_control = list(factor = c(4, 2), exact_refine = TRUE)
# )
# summary(fit)

## ----slim-fit, eval = FALSE---------------------------------------------------
# fit_slim <- slim_terradish(fit)
# fit_slim$storage
# 
# saveRDS(fit_slim, "terradish-fit-slim.rds")
# 
# # Equivalent shorthand at fit time
# fit_slim <- terradish(
#   admissible_squared_distance ~ cov1 + cov2,
#   surface,
#   measurement_model = generalized_wishart,
#   nu = effective_marker_df,
#   slim = TRUE
# )

## ----quick-ref, eval = FALSE--------------------------------------------------
# # 1. Build the surface as usual
# surface <- conductance_surface(covariates, coords, directions = 8)
# 
# # 2. Fit at scale: use amg once the raster is too large for the direct solver
# # Generalized Wishart requires an admissible squared-distance response.
# fit <- terradish(admissible_squared_distance ~ cov1 + cov2, surface,
#                  conductance_model = loglinear_conductance,
#                  measurement_model = generalized_wishart,
#                  nu = effective_marker_df,
#                  solver    = "amg",            # large-N linear solver
#                  curvature = "gauss_newton")   # information-based standard errors
# 
# summary(fit)                                   # conditional estimates and SEs
# 
# # Optional archival copy when prediction from the saved object is unnecessary
# fit_archive <- slim_terradish(fit)
# saveRDS(fit_archive, "fit-archive.rds")
# 

