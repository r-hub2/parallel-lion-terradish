## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
                      fig.width = 6, fig.height = 4.5, warning = FALSE, message = FALSE)
suppressPackageStartupMessages(library(terradish))
res <- tryCatch(readRDS("vignette-directional.rds"), error = function(e) NULL)
have_res <- !is.null(res)
if (have_res) {
  DIM <- res$example$DIM
  r <- terra::rast(nrows = DIM, ncols = DIM, xmin = 0, xmax = DIM,
                   ymin = 0, ymax = DIM)
  covs <- c(
    terra::setValues(r, res$example$v1_values),
    terra::setValues(r, res$example$elev_values)
  )
  names(covs) <- c("v1", "elev")
  coords <- terra::xyFromCell(r, res$example$focal_cells)
  surface <- conductance_surface(covs, coords, directions = 8)
  dir_cov <- edge_gradient(covs[["elev"]], surface)
  dir_cov$d <- cbind(elev = dir_cov$d)
  edge_rates <- directed_rates(res$fit, data = surface, directional = dir_cov)
}

## ----edges, eval = FALSE------------------------------------------------------
# library(terradish)
# library(terra)
# 
# # A coarse/deme-scale conductance surface is the recommended starting point.
# surface <- conductance_surface(covariates, coords, directions = 8)
# 
# # Directional covariate from an elevation layer (downhill = positive drop)
# dir_cov <- edge_gradient(elevation, surface)
# 
# # Optional, but helpful: name the directional covariate so the coefficient
# # appears as gamma_elev rather than gamma1 in printed output.
# dir_cov$d <- cbind(elev = dir_cov$d)

## ----fit, eval = FALSE--------------------------------------------------------
# fit <- terradish_directed(
#   S ~ habitat,                    # symmetric conductance covariate(s)
#   data              = surface,
#   directional       = dir_cov,    # directional covariate(s)
#   measurement_model = wishart_covariance,   # match this to the response type
#   nu                = nu,
#   solver            = "sparse_lu_cpp"
# )
# summary(fit)        # theta (symmetric) + gamma (directional), z-tests, AIC

## ----cached-summary, eval = have_res------------------------------------------
summary(res$fit)

## ----lrt, eval = FALSE--------------------------------------------------------
# # Refit the same model with the directional coefficient fixed at zero.
# fit_reversible <- terradish_directed(
#   S ~ habitat,
#   data              = surface,
#   directional       = dir_cov,
#   measurement_model = wishart_covariance,
#   nu                = nu,
#   gamma_bound       = 0,
#   solver            = "sparse_lu_cpp"
# )
# 
# lrt_stat <- 2 * (as.numeric(logLik(fit)) -
#                    as.numeric(logLik(fit_reversible)))
# 
# # The test df is the number of constrained directional coefficients.
# # Set this explicitly because the reduced object still carries a gamma slot
# # fixed at zero.
# lrt_df <- length(fit$gamma)
# lrt_p  <- pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)
# 
# data.frame(statistic = lrt_stat, df = lrt_df, p_value = lrt_p)

## ----lrt-table, echo = FALSE, eval = have_res---------------------------------
knitr::kable(
  res$lrt_table,
  digits = 3,
  caption = "Likelihood-ratio comparison of the directed model against the reversible special case."
)

## ----recovery-table, echo = FALSE, eval = have_res----------------------------
knitr::kable(
  res$coef_table,
  digits = 3,
  caption = "Directed-model recovery on simulated data."
)

## ----directional-visuals, eval = FALSE----------------------------------------
# # One row per undirected graph edge, with fitted rates in both directions.
# edge_rates <- directed_rates(fit, data = surface, directional = dir_cov)
# 
# # Arrows point from the lower-rate endpoint toward the higher-rate endpoint.
# # Width and color increase with |log(rate_ab / rate_ba)|.
# plot(fit, type = "directional", data = surface, directional = dir_cov)
# 
# # Overlay directional arrows on the fitted reversible conductance surface.
# plot(fit, type = "combined", data = surface, directional = dir_cov)

## ----edge-summary, echo = FALSE, eval = have_res------------------------------
knitr::kable(
  res$edge_summary,
  digits = 3,
  caption = "Summary of fitted directional edge-rate bias in the cached example."
)

## ----directional-plot, echo = FALSE, fig.width = 6, fig.height = 5, eval = have_res----
plot(res$fit, type = "directional", data = surface, directional = dir_cov,
     min_abs_log_ratio = 0.05)

## ----combined-plot, echo = FALSE, fig.width = 6, fig.height = 5, eval = have_res----
plot(res$fit, type = "combined", data = surface, directional = dir_cov,
     min_abs_log_ratio = 0.05)

## ----directional-filter, eval = FALSE-----------------------------------------
# plot(fit, type = "combined", data = surface, directional = dir_cov,
#      min_abs_log_ratio = 0.6)

## ----recovery-show, eval = FALSE----------------------------------------------
# fit <- terradish_directed(S ~ v1, data = surface, directional = dir_cov,
#                           measurement_model = wishart_covariance, nu = nu,
#                           solver = "sparse_lu_cpp")
# coef(fit)           # theta (~truth) and gamma (~truth, positive)

## ----grad-check, echo = FALSE, results = "asis", eval = have_res--------------
cat(sprintf("The reverse-mode (transpose-solve) gradient matches numerical differentiation to **%.0e**, and gamma is recovered with the correct (positive) sign.\n", res$grad_err))

## ----quickref, eval = FALSE---------------------------------------------------
# library(terradish)
# surface <- conductance_surface(covariates, coords, directions = 8)  # coarse/deme-scale
# dir_cov <- edge_gradient(elevation, surface)                        # directional covariate
# fit <- terradish_directed(S ~ habitat, data = surface, directional = dir_cov,
#                           measurement_model = wishart_covariance, nu = nu,
#                           solver = "sparse_lu_cpp")
# summary(fit)                       # theta (symmetric) + gamma (directional)
# plot(fit, data = surface)          # symmetric conductance surface
# edge_rates <- directed_rates(fit, data = surface, directional = dir_cov)
# head(edge_rates)                   # fitted edge rates and log-rate ratios
# plot(fit, type = "combined", data = surface, directional = dir_cov)

