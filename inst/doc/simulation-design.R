## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4.5,
  warning  = FALSE,
  message  = FALSE
)

## ----template-----------------------------------------------------------------
library(terradish)
library(terra)

data(melip)
melip.altitude    <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords      <- terra::unwrap(melip.coords)

covariates <- c(melip.altitude, melip.forestcover)
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(covariates)

surface <- conductance_surface(covariates, melip.coords,
                               directions = 8, saveStack = TRUE)

cat("Candidate focal sites:", length(surface$demes), "\n")

## ----truth--------------------------------------------------------------------
theta_true <- c(forestcover = 0.6, altitude = -0.4)
tau_true   <- 2   # strength of the landscape signal
sigma_true <- 1   # nugget (non-landscape variance)

## ----single-sim---------------------------------------------------------------
sim <- simulate_covariance_response(
  theta             = theta_true,
  formula           = ~ forestcover + altitude,
  data              = surface,
  conductance_model = loglinear_conductance,
  tau               = tau_true,
  sigma             = sigma_true,
  nu                = 100,     # 100 effective independent marker contributions
  seed              = 1
)

dim(sim$covariance)        # one row/column per focal site
round(sim$covariance[1:4, 1:4], 3)

## ----single-fit---------------------------------------------------------------
fit <- terradish(
  sim$covariance ~ forestcover + altitude,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = wishart_covariance,
  nu                = 100
)

rbind(true      = theta_true,
      estimated = round(coef(fit), 3))

## ----nu-fit-sensitivity, eval = FALSE-----------------------------------------
# # These runs use the same seeds and therefore the same simulated responses.
# # Increase nsim after checking runtime and convergence on a small pilot.
# design <- list(
#   theta = theta_true, formula = ~ forestcover + altitude, data = surface,
#   sample_sizes = 20, strategies = "spacefill",
#   tau = tau_true, sigma = sigma_true, nu = 100,
#   nsim = 100, seed = 905,
#   control = NewtonRaphsonControl(maxit = 100)
# )
# matched <- do.call(covariance_response_power, design)
# overstated <- do.call(covariance_response_power,
#                      c(design, list(nu_fit = 2000)))
# rbind(matched = matched$parameter_summary,
#       overstated = overstated$parameter_summary)

## ----power-run, eval = FALSE--------------------------------------------------
# power <- covariance_response_power(
#   theta             = theta_true,
#   formula           = ~ forestcover + altitude,
#   data              = surface,
#   sample_sizes      = c(10, 20, 34),
#   strategies        = c("spacefill", "random"),
#   conductance_model = loglinear_conductance,
#   fit_models = list(
#     full          = list(formula = ~ forestcover + altitude),
#     altitude_only = list(formula = ~ altitude)
#   ),
#   tau        = tau_true,
#   sigma      = sigma_true,
#   nu         = 100,
#   nu_fit     = 100,   # change separately to assess overstated precision
#   nsim       = 20,    # increase for smoother estimates in a real study
#   n_designs  = 2,     # random designs evaluated per sample size
#   seed       = 1,
#   control    = NewtonRaphsonControl(maxit = 100, verbose = FALSE)
# )
# 
# power

## ----power-summary, eval = FALSE----------------------------------------------
# power$summary[, c("strategy", "sample_size", "model",
#                   "fit_rate", "conductance_power",
#                   "mean_conductance_cor", "selected_AIC_rate")]

## ----power-params, eval = FALSE-----------------------------------------------
# power$parameter_summary[
#   power$parameter_summary$model == "full",
#   c("strategy", "sample_size", "parameter",
#     "power", "bias", "rmse", "coverage")
# ]

## ----quick-reference, eval = FALSE--------------------------------------------
# library(terradish)
# library(terra)
# 
# # 1. Build a landscape template (saveStack = TRUE enables site subsetting)
# surface <- conductance_surface(covariates, coords,
#                                directions = 8, saveStack = TRUE)
# 
# # 2. One simulated covariance matrix from specified generating parameters
# sim <- simulate_covariance_response(
#   theta = c(x1 = 0.6, x2 = -0.4), formula = ~ x1 + x2,
#   data = surface, tau = 2, sigma = 1, nu = 100, seed = 1
# )
# 
# # 3. Refit to check recovery
# fit <- terradish(sim$covariance ~ x1 + x2, data = surface,
#                  conductance_model = loglinear_conductance,
#                  measurement_model = wishart_covariance, nu = 100)
# coef(fit)
# 
# # 4. Full design assessment across sample sizes / strategies / models
# power <- covariance_response_power(
#   theta = c(x1 = 0.6, x2 = -0.4), formula = ~ x1 + x2,
#   data = surface, sample_sizes = c(10, 20, 30),
#   strategies = c("spacefill", "random"),
#   tau = 2, sigma = 1, nu = 100, nu_fit = 100, nsim = 200, seed = 1
# )
# power$summary            # fit rate, conductance recovery, model selection
# power$parameter_summary  # power, bias, RMSE, coverage per coefficient

