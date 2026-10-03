## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>",
  fig.width = 6, fig.height = 4.5, warning = FALSE, message = FALSE)

## ----landscape----------------------------------------------------------------
library(terradish)
library(terra)

# A small projected grid with one continuous landscape covariate.
r <- rast(nrows = 8, ncols = 8, xmin = 0, xmax = 8,
          ymin = 0, ymax = 8, crs = "EPSG:3857")
values(r) <- sin(seq_len(64) / 5) + seq_len(64) / 64
names(r) <- "x"
coords <- xyFromCell(r, c(1, 5, 8, 12, 20, 27, 33, 37, 45, 53, 57, 64))
surface <- conductance_surface(r, coords, directions = 8, saveStack = TRUE)

# Environmental values are site attributes. They need not be graph covariates.
environment <- data.frame(climate = sin(seq_len(nrow(coords))))
environment_pairs <- pairwise_endpoint_covariates(
  environment, transform = "absdiff", scale = TRUE
)
z <- pairwise_covariates(environment_pairs, geographic = dist(coords))
head(z)

## ----measurement-models-------------------------------------------------------
g_mlpe <- mlpe_covariates(z)
g_wishart <- wishart_covariates(z, model = "wishart_covariance")

## ----response-----------------------------------------------------------------
nu_demo <- 1000
base <- simulate_covariance_response(c(x = 0.8), ~x, surface,
  tau = 1, sigma = 0.1, nu = nu_demo, seed = 301)
kernels <- attr(g_wishart, "kernel_covariates")
Sigma <- base$Sigma + 0.1 * kernels[, , 1] + 0.05 * kernels[, , 2]
set.seed(301)
S_cov <- stats::rWishart(1, df = nu_demo, Sigma = Sigma)[, , 1] / nu_demo
S_dist <- dist_from_cov(S_cov)
check_distance_response(S_dist)$admissible

## ----fit-models---------------------------------------------------------------
ctrl <- NewtonRaphsonControl(maxit = 100, verbose = FALSE)
fit_ibd <- terradish(S_dist ~ 1, surface,
  measurement_model = mlpe, control = ctrl)
fit_ibe <- terradish(S_dist ~ 1, surface,
  measurement_model = g_mlpe, control = ctrl)
fit_ibr <- terradish(S_dist ~ x, surface,
  measurement_model = mlpe, control = ctrl)
fit_joint <- terradish(S_dist ~ x, surface,
  measurement_model = g_mlpe, control = ctrl)

fit_w_ibd <- terradish(S_cov ~ 1, surface,
  measurement_model = wishart_covariance, nu = nu_demo, control = ctrl)
fit_w_ibe <- terradish(S_cov ~ 1, surface,
  measurement_model = g_wishart, nu = nu_demo, control = ctrl)
fit_w_ibr <- terradish(S_cov ~ x, surface,
  measurement_model = wishart_covariance, nu = nu_demo, control = ctrl)
fit_w_joint <- terradish(S_cov ~ x, surface,
  measurement_model = g_wishart, nu = nu_demo, control = ctrl)

fits <- list(mlpe_ibd = fit_ibd, mlpe_uniform_ibe = fit_ibe,
  mlpe_ibr = fit_ibr, mlpe_joint = fit_joint,
  wishart_ibd = fit_w_ibd, wishart_uniform_ibe = fit_w_ibe,
  wishart_ibr = fit_w_ibr, wishart_joint = fit_w_joint)
vapply(fits, function(fit) fit$convergence$code, numeric(1))

## ----conditional-coefficients-------------------------------------------------
rbind(MLPE_without_pairs = coef(fit_ibr), MLPE_with_pairs = coef(fit_joint),
      Wishart_without_pairs = coef(fit_w_ibr), Wishart_with_pairs = coef(fit_w_joint))
summary(fit_joint)
summary(fit_w_joint)

## ----response-change----------------------------------------------------------
mlpe_response_change(fit_joint)

## ----ratios-------------------------------------------------------------------
ratios <- rbind(terradish_ibe_ratio(fit_joint),
                terradish_ibe_ratio(fit_w_joint))
ratios

## ----information-criteria-----------------------------------------------------
labels <- c("IBD", "IBD + IBE (uniform surface)", "IBR", "IBR + pairwise terms")
aic_table(list(fit_ibd, fit_ibe, fit_ibr, fit_joint), mod_names = labels)
aic_table(list(fit_w_ibd, fit_w_ibe, fit_w_ibr, fit_w_joint), mod_names = labels)

## ----nested-tests-------------------------------------------------------------
anova(fit_ibr, fit_joint)
anova(fit_w_ibr, fit_w_joint)

## ----nu-sensitivity-----------------------------------------------------------
fit_w_lower_nu <- terradish_rescale_nu(fit_w_joint, nu = 50)
rbind(primary = coef(fit_w_joint), sensitivity = coef(fit_w_lower_nu))
terradish_ibe_ratio(fit_w_lower_nu)

## ----fixed-cv, eval = FALSE---------------------------------------------------
# # Use the identical site folds for the plain and extended measurement models.
# folds <- terradish_folds(coords, k = 3, method = "random", seed = 42)
# cv_mlpe <- terradish_cv_folds(surface,
#   formulas = list(IBR = S_dist ~ x, joint = S_dist ~ x), folds = folds,
#   model = list(IBR = mlpe, joint = g_mlpe), nuisance = "fixed", control = ctrl)
# cv_wishart <- terradish_cv_folds(surface,
#   formulas = list(IBR = S_cov ~ x, joint = S_cov ~ x), folds = folds,
#   model = list(IBR = wishart_covariance, joint = g_wishart),
#   nu = nu_demo, nuisance = "fixed", control = ctrl)
# cv_mlpe$summary
# cv_wishart$summary

## ----fit-plot, fig.cap = "Observed and fitted pairwise distances for the joint MLPE model."----
plot(fit_joint, type = "fit")

## ----surface-plot, fig.cap = "Relative conductance conditional on the pairwise climate and geographic terms."----
plot(fit_joint, type = "surface", data = surface)

## ----marginal-plot, fig.cap = "Conditional conductance association with landscape x."----
plot(fit_joint, type = "marginal", data = surface)

## ----quick-reference, eval = FALSE--------------------------------------------
# # 1. Preserve identical transforms and site ordering across families.
# z <- pairwise_covariates(
#   pairwise_endpoint_covariates(environment, transform = "absdiff", scale = TRUE),
#   geographic = dist(coords)
# )
# g_mlpe <- mlpe_covariates(z)
# g_wishart <- wishart_covariates(z, model = "wishart_covariance")
# 
# # 2. Use covariance-derived distances when comparing ratio interpretations.
# S_dist <- dist_from_cov(S_cov)
# fit_m <- terradish(S_dist ~ x, surface, measurement_model = g_mlpe)
# fit_w <- terradish(S_cov ~ x, surface, measurement_model = g_wishart, nu = nu_demo)
# 
# # 3. Inspect convergence, conditional surface changes, and ratio uncertainty.
# summary(fit_m)
# summary(fit_w)
# rbind(terradish_ibe_ratio(fit_m), terradish_ibe_ratio(fit_w))
# 
# # 4. Assess predictive support and sensitivity to effective information.
# terradish_ibe_ratio(terradish_rescale_nu(fit_w, nu = 50))

