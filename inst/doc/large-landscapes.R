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

## ----reduce-------------------------------------------------------------------
# a conductance vector to reduce (evaluate the model at some theta)
model <- loglinear_conductance(~ forestcover + altitude, surface$x)
cond  <- model(c(0.3, -0.2))$conductance

red <- terradish_kron_reduce_tiled(surface, cond)
red$method                 # which strategy "auto" chose
red$peak$interior          # largest single factorization (vertices)
red$n_interior             # interior vertices in total

# forcing nested dissection bounds the largest factorization far below the
# whole interior (this is the path that keeps memory in check at scale)
nested <- terradish_kron_reduce_tiled(surface, cond, method = "nested")
nested$peak$interior

# both are identical to the single-shot reduction
single <- terradish_kron_reduce(surface, cond)
P <- match(single$boundary, red$boundary)
max(abs(as.matrix(red$laplacian[P, P]) -
        as.matrix(single$laplacian))) / max(abs(as.matrix(single$laplacian)))

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
# # 3. (Advanced) exact memory-bounded reduction onto the focal sites
# cond_map <- conductance(surface, fit)[["est"]]
# cond <- terra::values(cond_map, mat = FALSE)
# cond <- cond[is.finite(cond)]                   # active-cell conductance vector
# red  <- terradish_kron_reduce_tiled(surface, cond)  # auto: fast unless it would not fit
# red$method                                     # "direct" or "nested"

