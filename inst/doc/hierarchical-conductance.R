## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4.5,
  warning  = FALSE,
  message  = FALSE
)
# Results are pre-computed (model fitting is too slow for a CRAN vignette build).
# The chunk below loads a saved object; the fitting code is shown with eval=FALSE.
res <- tryCatch(
  readRDS("vignette-hierarchical.rds"),
  error = function(e) NULL
)
have_res <- !is.null(res)

## ----build-surface, eval = FALSE----------------------------------------------
# library(terradish)
# library(terra)
# 
# data(melip)
# covs <- c(terra::scale(terra::unwrap(melip.altitude)),
#           terra::scale(terra::unwrap(melip.forestcover)))
# names(covs) <- c("altitude", "forestcover")
# surface <- conductance_surface(covs, terra::unwrap(melip.coords), directions = 8)

## ----simulate, eval = FALSE---------------------------------------------------
# set.seed(20)
# 
# n  <- length(surface$demes)                       # 37 focal sites
# Vc <- surface$vertex_coordinates                  # cell centers of the graph
# X  <- as.matrix(surface$x[, c("altitude", "forestcover")])
# 
# # The UNMAPPED feature: a smooth Gaussian bump in log-conductance that no
# # raster in the model encodes. Centering it keeps it out of the intercept.
# xr <- range(Vc[, 1]); yr <- range(Vc[, 2])
# blob <- 1.6 * exp(-(((Vc[, 1] - (xr[1] + 0.35 * diff(xr)))^2 +
#                      (Vc[, 2] - (yr[1] + 0.6 * diff(yr)))^2) /
#                     (2 * (0.18 * sqrt(diff(xr)^2 + diff(yr)^2))^2)))
# blob <- blob - mean(blob)
# 
# # True log-conductance = 0.5*altitude - 0.6*forestcover + 1.0*blob.
# # .design_loglinear_model() builds a conductance model straight from a design
# # matrix, which lets us put the blob in the truth without giving it to the fit.
# E <- as.matrix(terradish_algorithm(
#   terradish:::.design_loglinear_model(cbind(X, blob)),
#   leastsquares, surface, S = diag(n), theta = c(0.5, -0.6, 1.0),
#   objective = FALSE, gradient = FALSE, hessian = FALSE,
#   partial = FALSE)$covariance)
# 
# # Draw a genetic covariance matrix from that truth. With real data this is
# # where cov_from_genetic_data() or cov_from_biallelic() would supply S, and nu
# # would be your effective marker count.
# nu <- 1500
# S  <- rWishart(1, df = nu, Sigma = (E + 0.25 * diag(n)) / nu)[, , 1]

## ----fit, eval = FALSE--------------------------------------------------------
# fit_h <- terradish_hierarchical(
#   S ~ altitude + forestcover,
#   data              = surface,
#   measurement_model = wishart_covariance,
#   nu                = 1500,
#   field_resolution  = 6,
#   tau2              = "reml"
# )
# 
# fit_plain <- terradish(
#   S ~ altitude + forestcover, data = surface,
#   conductance_model = loglinear_conductance,
#   measurement_model = wishart_covariance, nu = 1500
# )

## ----show-coef, echo = FALSE, eval = have_res---------------------------------
knitr::kable(res$coef_table, digits = 3,
             caption = "Conductance coefficients (true: altitude 0.50, forestcover -0.60)")

## ----coef-text, eval = FALSE--------------------------------------------------
# # With the field: altitude ~0.54, forestcover ~-0.63 (close to truth 0.50, -0.60)
# coef(fit_h)
# # Without the field: forestcover ~-0.43 in this particular simulation.
# coef(fit_plain)

## ----recovery-text, echo = FALSE, results = "asis", eval = have_res-----------
cat(sprintf(
  "In this controlled simulation, the fitted field correlates **%.2f** with the generating blob, and the selected field variance is **tau^2 = %.3g**.\n",
  res$cor_field_blob, res$tau2))

## ----field-map, eval = FALSE--------------------------------------------------
# field_raster <- conductance_field(fit_h, surface, type = "field")
# plot(field_raster)            # the smooth residual conductance surface u
# 
# # Or the full fitted surface, exp(X theta + u):
# plot(conductance_field(fit_h, surface, type = "conductance"))

## ----field-plot, echo = FALSE, eval = have_res, fig.cap = "Fitted smooth residual field u in the controlled simulation. The peak overlaps the omitted generating feature."----
if (!is.null(res$field_raster)) {
  terra::plot(terra::unwrap(res$field_raster), main = "Fitted field (u)")
}

## ----tau2-table, echo = FALSE, eval = have_res--------------------------------
knitr::kable(res$tau2_selection, digits = 3,
             caption = "Laplace marginal log-likelihood across the tau^2 grid (interior maximum).")

## ----quickref, eval = FALSE---------------------------------------------------
# library(terradish)
# 
# # 1. Build the conductance surface (as for any terradish workflow)
# surface <- conductance_surface(covariates, coords, directions = 8)
# 
# # 2. Fit the hierarchical model (covariates + smooth field), REML tau^2
# fit <- terradish_hierarchical(
#   S ~ forestcover + altitude,
#   data              = surface,
#   measurement_model = wishart_covariance,   # or generalized_wishart / mlpe
#   nu                = nu,
#   field_resolution  = 6,
#   tau2              = "reml"
# )
# 
# # 3. Inspect conditional coefficients and the residual-field diagnostic
# summary(fit)                                  # theta + tau^2 profile
# plot(conductance_field(fit, surface))         # the field u

