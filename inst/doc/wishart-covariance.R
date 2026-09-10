## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment  = "#>",
  fig.width  = 6,
  fig.height = 4.5,
  warning  = FALSE,
  message  = FALSE
)


## ----biallelic-example--------------------------------------------------------
library(terradish)

# Simulated counts: 4 populations, 5 biallelic loci
# Y[i, l] = number of derived alleles in population i at locus l
# N[i, l] = total haploid chromosomes sampled (here, 20 diploids → N = 40)
set.seed(42)
n_pop <- 4
n_loci <- 5
Y <- matrix(sample(0:40, n_pop * n_loci, replace = TRUE),
            nrow = n_pop, ncol = n_loci,
            dimnames = list(paste0("pop", seq_len(n_pop)),
                            paste0("snp", seq_len(n_loci))))
N <- matrix(40L, nrow = n_pop, ncol = n_loci)

S_biallelic <- cov_from_biallelic(Y, N)
round(S_biallelic, 3)

## ----feature-example----------------------------------------------------------
# Individual-level PCA scores: 12 individuals across 4 populations, 3 PCs
set.seed(7)
pca_scores <- matrix(rnorm(12 * 3), nrow = 12, ncol = 3,
                     dimnames = list(NULL, c("PC1", "PC2", "PC3")))
groups <- rep(c("pop1", "pop2", "pop3", "pop4"), each = 3)

S_pca <- cov_from_genetic_data(pca_scores, groups = groups)
round(S_pca, 3)

## ----microsatellite-example---------------------------------------------------
# 9 individuals, 2 microsatellite loci, each stored in two allele-call columns
alleles <- data.frame(
  loc1_a = c(100, 100, 102, 102, 104, 104, 100, 102, 104),
  loc1_b = c(100, 102, 102, 104, 104, 100, 104, 100, 102),
  loc2_a = c(200, 202, 200, 202, 204, 204, 200, 204, 202),
  loc2_b = c(202, 202, 204, 204, 204, 200, 202, 200, 204)
)
groups_ms <- rep(c("popA", "popB", "popC"), each = 3)
loci_ms <- c("loc1", "loc1", "loc2", "loc2")

S_ms <- cov_from_genetic_data(
  alleles,
  groups  = groups_ms,
  input   = "allele_calls",
  loci    = loci_ms
)
round(S_ms, 3)

## ----microsatellite-nu--------------------------------------------------------
# Larger sensitivity value: sum(K_l - 1) across loci
msat_nu_upper <- sum(vapply(split(seq_along(loci_ms), loci_ms), function(j) {
  alleles_l <- unlist(alleles[, j, drop = FALSE], use.names = FALSE)
  length(unique(alleles_l[!is.na(alleles_l)])) - 1L
}, integer(1)))

n_loci <- length(unique(loci_ms))   # conservative primary value

cat("Primary value (number of loci):             nu =", n_loci, "\n")
cat("Larger sensitivity value (sum(K_l - 1)):   nu =", msat_nu_upper, "\n")

## ----check-eigenvalues--------------------------------------------------------
ev <- eigen(S_ms, symmetric = TRUE, only.values = TRUE)$values
all(ev > 0)  # TRUE → positive definite; FALSE → semi-definite or indefinite

## ----landscape-setup----------------------------------------------------------
library(terra)

data(melip)
melip.altitude    <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords      <- terra::unwrap(melip.coords)

covariates <- c(melip.altitude, melip.forestcover)
names(covariates) <- c("altitude", "forestcover")
covariates <- scale_covariates(covariates)

surface <- conductance_surface(covariates, melip.coords, directions = 8)

n_sites <- nrow(melip.coords)
cat("Number of sampling sites:", n_sites, "\n")

## ----simulate-covariance------------------------------------------------------
# Use conductance parameters roughly consistent with the melip fit in the
# getting-started vignette (forestcover positive, altitude negative).
# nu = 1000 is an illustrative effective marker count for this simulation.
sim <- simulate_covariance_response(
  theta             = c(forestcover = 1.0, altitude = -0.5),
  formula           = ~ forestcover + altitude,
  data              = surface,
  conductance_model = loglinear_conductance,
  tau               = 1.0,
  sigma             = 0.5,
  nu                = 1000,
  seed              = 123
)

# This is what would come from cov_from_genetic_data() applied to real genotypes.
S_cov <- sim$covariance
dim(S_cov)
round(S_cov[1:5, 1:5], 4)

## ----real-data-sketch, eval = FALSE-------------------------------------------
# # With biallelic SNP count matrices (Y = derived allele counts, N = haploid sample sizes)
# S_cov <- cov_from_biallelic(Y, N)
# 
# # With individual dosage / PCA scores grouped by population
# S_cov <- cov_from_genetic_data(dosage_matrix, groups = pop_vector,
#                                 normalize = "features")
# 
# # With microsatellite allele calls
# S_cov <- cov_from_genetic_data(allele_df, groups = pop_vector,
#                                 input = "allele_calls",
#                                 loci  = locus_id_vector)

## ----fit-wishart-covariance---------------------------------------------------
fit_wc <- terradish(
  S_cov ~ forestcover + altitude,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = wishart_covariance,
  nu                = 1000
)

summary(fit_wc)

## ----nu-invariance------------------------------------------------------------
fit_nu <- function(nu) {
  f <- terradish(
    S_cov ~ forestcover + altitude,
    data              = surface,
    conductance_model = loglinear_conductance,
    measurement_model = wishart_covariance,
    nu                = nu
  )
  list(theta = coef(f),
       se    = sqrt(diag(solve(f$fit$hessian))))
}

f100   <- fit_nu(100)
f1000  <- fit_nu(1000)
f10000 <- fit_nu(10000)

# Point estimates agree across three orders of magnitude
rbind(`nu=100`   = f100$theta,
      `nu=1000`  = f1000$theta,
      `nu=10000` = f10000$theta)

## ----nu-se--------------------------------------------------------------------
rbind(`nu=100`   = f100$se,
      `nu=1000`  = f1000$se,
      `nu=10000` = f10000$se)

# Ratio confirms the 1/sqrt(nu) scaling: a 10x increase in nu
# shrinks SEs by sqrt(10) ~ 3.16
f100$se / f1000$se

## ----recovery-----------------------------------------------------------------
cat("True conductance parameters:\n")
cat("  forestcover:", 1.0, "  altitude:", -0.5, "\n\n")

cat("Estimated conductance parameters:\n")
print(round(coef(fit_wc), 3))

cat("\nTrue nuisance parameters (on the variance scale):\n")
cat("  tau   =", sim$tau, "\n")
cat("  sigma =", sim$sigma, " (nugget variance; fitted model stores log(sigma))\n\n")

cat("Estimated nuisance parameters (phi vector, log scale for sigma):\n")
phi_est <- round(fit_wc$fit$phi, 3)
print(phi_est)

cat("\nFor comparison, exp(phi['sigma']) =",
    round(exp(fit_wc$fit$phi["sigma"]), 3),
    " vs. true sigma =", sim$sigma, "\n")

## ----plot-fit, fig.cap = "***Observed vs. fitted genetic covariance.***"------
plot(fit_wc, type = "fit")

## ----plot-surface, fig.cap = "***Estimated conductance surface with 95% CI bounds.***"----
plot(fit_wc, type = "surface", data = surface)

## ----plot-marginal, fig.cap = "***Marginal associations on the fitted conductance scale.***"----
plot(fit_wc, type = "marginal", data = surface)

## ----wishart-support-clamp, eval = FALSE--------------------------------------
# plot(
#   fit_wc,
#   type = "marginal",
#   data = surface,
#   support = "focal",
#   support_probs = c(0.01, 0.99),
#   clamp_covariates = c("forestcover", "altitude")
# )
# 
# cond_wc_focal <- conductance(
#   surface,
#   fit_wc,
#   support = "focal",
#   support_probs = c(0.01, 0.99),
#   clamp_covariates = c("forestcover", "altitude")
# )

## ----drift-simulate-----------------------------------------------------------
# An illustrative site-level covariate for diagonal variance.
n_site <- length(surface$demes)
site_variance <- scale(seq_len(n_site))[, 1]

# True conductance-implied covariance E at the conductance truth used above.
E_true <- as.matrix(terradish_algorithm(
  loglinear_conductance(~ forestcover + altitude, surface$x),
  leastsquares, surface, S = diag(n_site),
  theta = c(forestcover = 1.0, altitude = -0.5),
  objective = FALSE, gradient = FALSE, hessian = FALSE, partial = FALSE
)$covariance)

# Per-site nugget: baseline exp(log(0.2)), slope -0.6 on the (centered) covariate.
Z <- cbind(1, scale(site_variance, scale = FALSE))
gamma_true <- c(log(0.2), -0.6)
nugget     <- as.vector(exp(Z %*% gamma_true))
Sigma_true <- 1.0 * E_true + diag(nugget)

set.seed(123)
nu_demo <- 500
S_drift <- rWishart(1, df = nu_demo, Sigma = Sigma_true / nu_demo)[, , 1]

## ----drift-fit----------------------------------------------------------------
g_drift <- wishart_drift_covariates(site_variance, model = "wishart_covariance")

fit_drift <- terradish(
  S_drift ~ forestcover + altitude,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = g_drift,
  nu                = nu_demo
)

fit_scalar <- terradish(
  S_drift ~ forestcover + altitude,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = wishart_covariance,
  nu                = nu_demo
)

summary(fit_drift)

## ----drift-interpret----------------------------------------------------------
# Recovered log-diagonal slope (true value -0.6)
fit_drift$fit$phi["gamma_var1", 1]

# Compare only because the response, sites, graph, likelihood, and nu match.
aic_table(
  list(fit_drift, fit_scalar),
  mod_names = c("diagonal covariate", "scalar nugget")
)

## ----check-gw-response--------------------------------------------------------
S_dist <- dist_from_cov(S_cov)
gw_check <- check_distance_response(S_dist)
gw_check[c("admissible", "relative_min_eigenvalue", "n_negative")]

## ----ade4-cross-check, eval = FALSE-------------------------------------------
# # D is already squared, so take its square root for ade4::is.euclid().
# ade4::is.euclid(as.dist(sqrt(D)))

## ----compare-gw---------------------------------------------------------------
fit_gw <- terradish(
  S_dist ~ forestcover + altitude,
  data              = surface,
  conductance_model = loglinear_conductance,
  measurement_model = generalized_wishart,
  nu                = 1000
)

# Side-by-side parameter estimates
cat("wishart_covariance estimates:\n")
print(round(coef(fit_wc), 3))

cat("\ngeneralized_wishart estimates:\n")
print(round(coef(fit_gw), 3))

## ----full-workflow, eval = FALSE----------------------------------------------
# library(terradish)
# library(terra)
# 
# # 1. Load landscape data and build conductance surface
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
# # 2. Prepare the genetic covariance matrix from raw genotype data
# #    (replace with your actual data and effective nu)
# 
# # Option A: biallelic SNPs (Y = derived allele counts, N = haploid totals)
# S <- cov_from_biallelic(Y, N)
# 
# # Option B: individual-level numeric features grouped by population
# S <- cov_from_genetic_data(dosage_matrix, groups = pop_vector,
#                             normalize = "features")
# 
# # Option C: microsatellite allele calls
# S <- cov_from_genetic_data(allele_df, groups = pop_vector,
#                             input = "allele_calls",
#                             loci  = locus_id_vector)
# 
# nu <- effective_marker_df  # justify, report, and evaluate sensitivity
# 
# # 3. Fit with wishart_covariance
# fit <- terradish(
#   S ~ forestcover + altitude,
#   data              = surface,
#   conductance_model = loglinear_conductance,
#   measurement_model = wishart_covariance,
#   nu                = nu
# )
# 
# # 4. Summarize and visualize
# summary(fit)
# coef(fit)
# 
# plot(fit, type = "fit")
# plot(fit, type = "surface", data = surface)
# plot(fit, type = "marginal", data = surface)   # conductance-scale associations

