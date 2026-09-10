## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width  = 6,
  fig.height = 4.5,
  warning   = FALSE,
  message   = FALSE
)


## ----load---------------------------------------------------------------------
library(terradish)
library(terra)

data(melip)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.altitude <- terra::unwrap(melip.altitude)
melip.coords <- terra::unwrap(melip.coords)

keep <- 1:8
melip.Fst <- melip.Fst[keep, keep, drop = FALSE]
coords <- melip.coords[keep]

## ----surface------------------------------------------------------------------
covariates <- melip.forestcover
names(covariates) <- "forestcover"

surface <- conductance_surface(
  covariates,
  coords,
  directions = 8,
  saveStack = TRUE
)

## ----gaussian-model-----------------------------------------------------------
gaussian_model <- gaussian_smoothed_loglinear_conductance(surface)

## ----fit----------------------------------------------------------------------
fit_gaussian <- terradish(
  melip.Fst ~ forestcover,
  data = surface,
  conductance_model = gaussian_model,
  measurement_model = mlpe,
  optimizer = "auto",
  leverage = FALSE,
  control = NewtonRaphsonControl(maxit = 12, verbose = FALSE)
)

fit_gaussian

## ----gaussian-coarse-start, eval=FALSE----------------------------------------
# fit_gaussian_coarse <- terradish(
#   melip.Fst ~ forestcover,
#   data = surface,
#   conductance_model = gaussian_model,
#   measurement_model = mlpe,
#   approximation = "coarse_raster",
#   approximation_control = list(
#     factor = c(4, 2),
#     exact_refine = TRUE
#   ),
#   optimizer = "auto",
#   leverage = FALSE,
#   control = NewtonRaphsonControl(maxit = 12, verbose = FALSE)
# )

## ----coef---------------------------------------------------------------------
coef(fit_gaussian)

## ----sigma--------------------------------------------------------------------
coef(fit_gaussian)[["sigma.forestcover"]]

## ----sigma-summary------------------------------------------------------------
gaussian_scale_summary(fit_gaussian)

## ----sigma-plot, fig.cap = "***Fitted Gaussian kernel with uncertainty band. The dashed line marks the distance containing 90% of the one-dimensional kernel mass.***"----
plot(fit_gaussian, type = "sigma")

## ----sigma-plot-km, fig.cap = "***Fitted Gaussian kernel with x-axis relabeled in kilometers.***"----
plot(
  fit_gaussian,
  type = "sigma",
  distance_per_map_unit = 0.001,
  distance_unit = "km"
)

## ----fixed-fit----------------------------------------------------------------
fit_fixed <- terradish(
  melip.Fst ~ forestcover,
  data = surface,
  conductance_model = loglinear_conductance,
  measurement_model = mlpe,
  optimizer = "auto",
  leverage = FALSE,
  control = NewtonRaphsonControl(maxit = 12, verbose = FALSE)
)

fit_fixed

## ----compare-table------------------------------------------------------------
aic_table(
  list(fit_gaussian, fit_fixed),
  mod_names = c("Gaussian scale-aware", "Fixed raster")
)

## ----summary------------------------------------------------------------------
summary(fit_gaussian)

## ----conductance--------------------------------------------------------------
fitted_conductance <- conductance(surface, fit_gaussian)
fitted_conductance

## ----conductance-plot, fig.cap = "***Fitted conductance surface from the Gaussian scale-aware model, with sampling sites shown in red.***"----
plot(fitted_conductance[[1]], main = "Fitted conductance")
points(coords, pch = 19, col = "red")

## ----gaussian-support-clamp, eval = FALSE-------------------------------------
# fitted_conductance_focal <- conductance(
#   surface,
#   fit_gaussian,
#   support = "focal",
#   support_probs = c(0.01, 0.99),
#   clamp_covariates = "forestcover"
# )

## ----surface-two--------------------------------------------------------------
covariates_two <- c(melip.altitude, melip.forestcover)
names(covariates_two) <- c("altitude", "forestcover")

surface_two <- conductance_surface(
  covariates_two,
  coords,
  directions = 8,
  saveStack = TRUE
)

## ----fit-two------------------------------------------------------------------
fit_gaussian_two <- terradish(
  melip.Fst ~ forestcover + altitude,
  data = surface_two,
  conductance_model = gaussian_smoothed_loglinear_conductance(surface_two),
  measurement_model = mlpe,
  optimizer = "auto",
  leverage = FALSE,
  control = NewtonRaphsonControl(maxit = 12, verbose = FALSE)
)

fit_gaussian_two

## ----sigma-two----------------------------------------------------------------
gaussian_scale_summary(fit_gaussian_two)

## ----qref-gs, eval = FALSE----------------------------------------------------
# library(terradish)
# 
# # 1. Build the surface from the ORIGINAL raster; do not pre-smooth it
# surface <- conductance_surface(covariates, coords, directions = 8,
#                                saveStack = TRUE)
# 
# # 2. Scale-aware conductance model: sigma is estimated, not fixed. The model
# #    is built from the surface itself; sigma_upper bounds the search in map
# #    units, and defaults to the retained raster extent's diagonal when left NULL.
# cm <- gaussian_smoothed_loglinear_conductance(surface)
# 
# # 3. Fit; theta and sigma.<layer> are estimated jointly
# fit <- terradish(melip.Fst ~ altitude, data = surface,
#                  conductance_model = cm, measurement_model = mlpe)
# 
# # 4. Read the fitted scale. `sigma` is in map units, so reproject a
# #    longitude/latitude raster first or it will be reported in degrees.
# gaussian_scale_summary(fit)   # sigma + kernel-mass distances in map units
# 
# # 5. Sanity check against a fixed-raster fit of the same covariate
# fit_fixed <- terradish(melip.Fst ~ altitude, data = surface,
#                        conductance_model = loglinear_conductance,
#                        measurement_model = mlpe)
# aic_table(list(fit, fit_fixed),
#           mod_names = c("Scale-aware", "Fixed raster"))
# 
# # 6. Recover the conductance surface at the estimated scale
# plot(conductance(surface, fit))

