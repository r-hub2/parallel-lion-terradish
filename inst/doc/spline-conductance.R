## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width  = 7,
  fig.height = 5,
  warning   = FALSE,
  message   = FALSE
)


## ----basis-demo, fig.cap = "***Natural spline (ns) and B-spline (bs) basis functions for a standardized covariate.*** Each colored line is one basis column. The fitted conductance association is a weighted combination of all basis functions in a panel.", fig.height = 5----
library(splines)
library(terradish)

x_seq <- seq(-2, 2, length.out = 300)

# Natural spline basis with df = 4
B_ns <- ns(x_seq, df = 4)

# B-spline basis with k = 4 columns, degree 3 (cubic)
B_bs <- bs(x_seq, df = 4, degree = 3)

par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

matplot(x_seq, B_ns, type = "l", lty = 1, lwd = 2,
        col = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
        xlab = "Standardized covariate", ylab = "Basis value",
        main = "Natural spline (ns, df = 4)",
        ylim = c(-0.4, 1))
abline(h = 0, col = "grey70", lty = 2)
legend("topright", legend = paste0("basis ", 1:4),
       col = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
       lty = 1, lwd = 2, bty = "n", cex = 0.8)

matplot(x_seq, B_bs, type = "l", lty = 1, lwd = 2,
        col = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
        xlab = "Standardized covariate", ylab = "Basis value",
        main = "B-spline (bs, df = 4, degree = 3)",
        ylim = c(-0.4, 1))
abline(h = 0, col = "grey70", lty = 2)
legend("topright", legend = paste0("basis ", 1:4),
       col = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"),
       lty = 1, lwd = 2, bty = "n", cex = 0.8)

## ----load-data----------------------------------------------------------------
library(terradish)
library(terra)

data(melip)
melip.altitude    <- terra::unwrap(melip.altitude)
melip.forestcover <- terra::unwrap(melip.forestcover)
melip.coords      <- terra::unwrap(melip.coords)

cat("Number of sampling sites:", nrow(melip.coords), "\n")
cat("Dimension of Fst matrix: ", dim(melip.Fst), "\n")

## ----prepare-covariates-------------------------------------------------------
covariates <- c(melip.altitude, melip.forestcover)
names(covariates) <- c("altitude", "forestcover")

# Coarsen the demonstration graph so the many candidate fits remain suitable
# for a package vignette. Use a scientifically justified resolution in an
# applied analysis.
covariates <- terra::aggregate(covariates, fact = 3, na.rm = TRUE)

# Center and scale each layer to mean 0, sd 1
covariates <- scale_covariates(covariates)

# Build the conductance surface graph
# saveStack = TRUE is required for conductance() and plot(type = "surface")
surface <- conductance_surface(
  covariates,
  melip.coords,
  directions = 8,
  saveStack  = TRUE
)

cat("Graph vertices:", nrow(surface$x), "\n")
head(surface$x)

## ----assess-spline, eval = FALSE----------------------------------------------
# # Assessment for the altitude natural-spline model (df = 4)
# # Uses the same surface built above (saveStack = TRUE required for solver probes)
# assessment_spline <- terradish_assess_settings(
#   melip.Fst ~ forestcover + s(altitude, df = 4),
#   data              = surface,
#   conductance_model = smooth_loglinear_conductance,
#   measurement_model = mlpe,
#   probe_maxit       = 2,
#   verbose           = TRUE
# )
# 
# assessment_spline

## ----assess-spline-use, eval = FALSE------------------------------------------
# # Extract recommended optimizer and control once
# rec <- assessment_spline$recommended
# 
# fit_ns_assessed <- terradish(
#   melip.Fst ~ forestcover + s(altitude, df = 4),
#   data              = surface,
#   conductance_model = smooth_loglinear_conductance,
#   measurement_model = mlpe,
#   optimizer         = rec$optimizer,
#   control           = rec$control,
#   solver            = rec$solver,
#   solver_control    = rec$solver_control
# )

## ----fit-models---------------------------------------------------------------
assessment_spline_light <- terradish_assess_settings(
  melip.Fst ~ forestcover + s(altitude, df = 4),
  data              = surface,
  conductance_model = smooth_loglinear_conductance,
  measurement_model = mlpe,
  probe_maxit       = 1,
  verbose           = FALSE
)
rec <- assessment_spline_light$recommended

# 1. Log-linear baseline (one parameter per covariate)
fit_linear <- terradish(
  melip.Fst ~ altitude + forestcover,
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

# 2. Natural spline on altitude, linear forest cover (df = 4)
fit_ns_alt <- terradish(
  melip.Fst ~ forestcover + s(altitude, df = 4),
  data              = surface,
  conductance_model = smooth_loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

# 3. Natural splines on both covariates (df = 4 each)
fit_ns_both <- terradish(
  melip.Fst ~ s(altitude, df = 4) + s(forestcover, df = 4),
  data              = surface,
  conductance_model = smooth_loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

# 4. Cubic B-spline on altitude, linear forest cover (k = 4 columns, degree 3)
fit_bs_alt <- terradish(
  melip.Fst ~ forestcover + s(altitude, k = 4, basis = "bs", degree = 3),
  data              = surface,
  conductance_model = smooth_loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

## ----inspect-design-----------------------------------------------------------
# Build the smooth model object for the altitude-spline specification
smooth_spec <- smooth_loglinear_conductance(
  ~ forestcover + s(altitude, df = 4),
  surface$x
)

# Starting parameter vector (one entry per design column)
cat("Parameter names:\n")
names(attr(smooth_spec, "default"))

# Smooth term details
cat("\nSmooth term metadata:\n")
str(attr(smooth_spec, "smooth_loglinear_info"))

## ----show-coefs---------------------------------------------------------------
cat("--- Linear baseline ---\n")
print(round(coef(fit_linear), 4))

cat("\n--- Altitude ns spline ---\n")
print(round(coef(fit_ns_alt), 4))

cat("\n--- Both ns splines ---\n")
print(round(coef(fit_ns_both), 4))

cat("\n--- Altitude bs spline ---\n")
print(round(coef(fit_bs_alt), 4))

## ----aic-table----------------------------------------------------------------
aic_table(
  list(fit_linear, fit_ns_alt, fit_bs_alt, fit_ns_both),
  mod_names = c(
    "Linear",
    "Altitude ns (df=4)",
    "Altitude bs (k=4)",
    "Both ns (df=4 each)"
  ),
  AICc = TRUE
)

## ----marginal-ns-vs-bs, fig.cap = "***Altitude marginal association on the conductance scale.*** Left: natural spline (ns). Right: B-spline (bs). Both use df/k = 4 and degree 3.", fig.height = 4.5----
library(ggplot2)

p_ns <- plot(fit_ns_alt, type = "marginal", data = surface, n = 80,
             marginal_covariates = "altitude") +
  ggtitle("Altitude: ns (df = 4)")

p_bs <- plot(fit_bs_alt, type = "marginal", data = surface, n = 80,
             marginal_covariates = "altitude") +
  ggtitle("Altitude: bs (k = 4)")

# patchwork is optional; fall back to one panel per page without it
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  p_ns | p_bs
} else {
  print(p_ns)
  print(p_bs)
}

## ----marginal-all, fig.cap = "***Altitude marginal association across all four candidate models.*** The linear model constrains the association to a straight line (on the log-conductance scale). Both spline models can represent non-monotone shapes.", fig.height = 6----
p_lin   <- plot(fit_linear,  type = "marginal", data = surface, n = 80,
                marginal_covariates = "altitude") +
  ggtitle("Linear baseline")

p_ns_a  <- plot(fit_ns_alt,  type = "marginal", data = surface, n = 80,
                marginal_covariates = "altitude") +
  ggtitle("Altitude ns (df = 4)")

p_ns_b  <- plot(fit_ns_both, type = "marginal", data = surface, n = 80,
                marginal_covariates = "altitude") +
  ggtitle("Both ns (df = 4): altitude panel")

p_bs_a  <- plot(fit_bs_alt,  type = "marginal", data = surface, n = 80,
                marginal_covariates = "altitude") +
  ggtitle("Altitude bs (k = 4)")

(p_lin | p_ns_a) / (p_ns_b | p_bs_a)

## ----marginal-response, fig.cap = "***Marginal response-scale associations from the altitude natural-spline model.*** Shaded bands are approximate pointwise 95% predictive intervals. Left: natural spline. Right: B-spline.", fig.height = 4.5----
p_resp_ns <- plot(fit_ns_alt, type = "marginal_response", data = surface,
                  n = 80, marginal_covariates = "altitude") +
  ggtitle("ns (df = 4): response scale")

p_resp_bs <- plot(fit_bs_alt, type = "marginal_response", data = surface,
                  n = 80, marginal_covariates = "altitude") +
  ggtitle("bs (k = 4): response scale")

p_resp_ns | p_resp_bs

## ----spline-support-clamp, eval = FALSE---------------------------------------
# plot(
#   fit_ns_alt,
#   type = "marginal",
#   data = surface,
#   n = 80,
#   support = "focal",
#   support_probs = c(0.01, 0.99),
#   clamp_covariates = "altitude"
# )
# 
# cond_ns_alt_focal <- conductance(
#   surface,
#   fit_ns_alt,
#   support = "focal",
#   support_probs = c(0.01, 0.99),
#   clamp_covariates = c("altitude", "forestcover")
# )

## ----obs-vs-fitted, fig.cap = "***Observed vs. fitted pairwise Fst for all candidate models.*** Each point is a pair of sampling sites. The diagonal (dashed) represents a perfect fit.", fig.height = 5----
p_fit_linear <- plot(fit_linear,  type = "fit", main = "Linear baseline")
p_fit_ns_alt <- plot(fit_ns_alt,  type = "fit", main = "Altitude ns (df = 4)")
p_fit_ns_both <- plot(fit_ns_both, type = "fit", main = "Both ns (df = 4)")
p_fit_bs_alt <- plot(fit_bs_alt,  type = "fit", main = "Altitude bs (k = 4)")

(p_fit_linear | p_fit_ns_alt) / (p_fit_ns_both | p_fit_bs_alt)

## ----surface-ns-alt, fig.cap = "***Estimated conductance surface from the altitude natural-spline model.*** The three panels show the point estimate and 95% CI bounds. Red points are sampling sites.", fig.height = 4----
cond_ns_alt <- conductance(surface, fit_ns_alt)

plot(cond_ns_alt,
     main = c("ns(altitude): estimate",
              "ns(altitude): lower 95%",
              "ns(altitude): upper 95%"))
points(melip.coords, pch = 19, cex = 0.6, col = "red")

## ----surface-compare, fig.cap = "***Point-estimate conductance surfaces: linear model (left) vs. altitude natural-spline model (right).*** Differences in surface texture reflect the non-linear altitude association represented by the spline.", fig.height = 4----
cond_linear <- conductance(surface, fit_linear)

par(mfrow = c(1, 2))
plot(cond_linear[[1]], main = "Linear baseline")
points(melip.coords, pch = 19, cex = 0.6, col = "red")

plot(cond_ns_alt[[1]], main = "Altitude ns (df = 4)")
points(melip.coords, pch = 19, cex = 0.6, col = "red")

## ----ns-df-comparison, fig.cap = "***Natural spline altitude association for df = 2, 3, 4, and 6.*** Higher df allows more shape complexity. Compare the smoothness and ecological plausibility of each curve before selecting a value.", fig.height = 6----
fit_ns2 <- terradish(
  melip.Fst ~ forestcover + s(altitude, df = 2),
  data              = surface,
  conductance_model = smooth_loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

fit_ns3 <- terradish(
  melip.Fst ~ forestcover + s(altitude, df = 3),
  data              = surface,
  conductance_model = smooth_loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

fit_ns6 <- terradish(
  melip.Fst ~ forestcover + s(altitude, df = 6),
  data              = surface,
  conductance_model = smooth_loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

# Use marginal_covariates = "altitude" to pull only the altitude panel from
# each fit, then tile the four ggplot objects with patchwork.
p2 <- plot(fit_ns2,    type = "marginal", data = surface, n = 80,
           marginal_covariates = "altitude") +
  ggtitle("ns, df = 2 (quadratic-like)")

p3 <- plot(fit_ns3,    type = "marginal", data = surface, n = 80,
           marginal_covariates = "altitude") +
  ggtitle("ns, df = 3")

p4 <- plot(fit_ns_alt, type = "marginal", data = surface, n = 80,
           marginal_covariates = "altitude") +
  ggtitle("ns, df = 4")

p6 <- plot(fit_ns6,    type = "marginal", data = surface, n = 80,
           marginal_covariates = "altitude") +
  ggtitle("ns, df = 6")

if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  (p2 | p3) / (p4 | p6)
} else {
  for (p in list(p2, p3, p4, p6)) print(p)
}

## ----ns-df-aic----------------------------------------------------------------
aic_table(
  list(fit_linear, fit_ns2, fit_ns3, fit_ns_alt, fit_ns6),
  mod_names = c("Linear", "ns df=2", "ns df=3", "ns df=4", "ns df=6"),
  AICc = TRUE
)

## ----bs-degree-comparison, fig.cap = "***B-spline altitude association for degree 1, 2, and 3 (k = 4 basis columns each).*** Degree 1 is piecewise linear; degree 3 is cubic. Smoother curves do not always fit better.", fig.height = 4----
fit_bs1 <- terradish(
  melip.Fst ~ forestcover + s(altitude, k = 4, basis = "bs", degree = 1),
  data              = surface,
  conductance_model = smooth_loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

fit_bs2 <- terradish(
  melip.Fst ~ forestcover + s(altitude, k = 4, basis = "bs", degree = 2),
  data              = surface,
  conductance_model = smooth_loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

pd1 <- plot(fit_bs1,    type = "marginal", data = surface, n = 80,
            marginal_covariates = "altitude") +
  ggtitle("bs, degree = 1 (piecewise linear)")

pd2 <- plot(fit_bs2,    type = "marginal", data = surface, n = 80,
            marginal_covariates = "altitude") +
  ggtitle("bs, degree = 2 (quadratic)")

pd3 <- plot(fit_bs_alt, type = "marginal", data = surface, n = 80,
            marginal_covariates = "altitude") +
  ggtitle("bs, degree = 3 (cubic)")

pd1 | pd2 | pd3

## ----mixed-terms--------------------------------------------------------------
# Forest cover linear, altitude as natural spline
fit_mixed <- terradish(
  melip.Fst ~ forestcover + s(altitude, df = 3),
  data              = surface,
  conductance_model = smooth_loglinear_conductance,
  measurement_model = mlpe,
  optimizer         = rec$optimizer,
  control           = rec$control,
  solver            = rec$solver,
  solver_control    = rec$solver_control,
  approximation     = rec$approximation,
  approximation_control = rec$approximation_control
)

coef(fit_mixed)

## ----mixed-marginal, fig.cap = "***Marginal associations from the mixed model (linear forest cover + altitude natural spline, df = 3).*** Left: the linear forest-cover association (log-conductance scale). Right: the non-linear altitude association (smooth curve).", fig.height = 4----
p_fc  <- plot(fit_mixed, type = "marginal", data = surface, n = 80,
              marginal_covariates = "forestcover") +
  ggtitle("Forest cover (linear term)")

p_alt <- plot(fit_mixed, type = "marginal", data = surface, n = 80,
              marginal_covariates = "altitude") +
  ggtitle("Altitude: ns (df = 3)")

p_fc | p_alt

## ----final-aic----------------------------------------------------------------
aic_table(
  list(fit_linear, fit_ns2, fit_ns3, fit_ns_alt, fit_ns6,
       fit_bs_alt, fit_ns_both, fit_mixed),
  mod_names = c(
    "Linear",
    "Altitude ns (df=2)",
    "Altitude ns (df=3)",
    "Altitude ns (df=4)",
    "Altitude ns (df=6)",
    "Altitude bs (k=4, deg=3)",
    "Both ns (df=4 each)",
    "Altitude ns (df=3) + linear forestcover"
  ),
  AICc = TRUE
)

## ----workflow, eval = FALSE---------------------------------------------------
# # Step 1: Fit a linear baseline: this is your reference
# assessment <- terradish_assess_settings(
#   S ~ covar2 + s(covar1, df = 3),
#   data              = surface,
#   conductance_model = smooth_loglinear_conductance,
#   measurement_model = mlpe,
#   probe_maxit       = 1,
#   verbose           = FALSE
# )
# rec <- assessment$recommended
# 
# fit_linear <- terradish(
#   S ~ covar1 + covar2,
#   data              = surface,
#   conductance_model = loglinear_conductance,
#   measurement_model = mlpe,
#   optimizer         = rec$optimizer,
#   control           = rec$control,
#   solver            = rec$solver,
#   solver_control    = rec$solver_control,
#   approximation     = rec$approximation,
#   approximation_control = rec$approximation_control
# )
# 
# # Step 2: Inspect the linear marginal associations.
# # Which covariate might plausibly have a non-linear association?
# # Use marginal_covariates to view one covariate at a time.
# plot(fit_linear, type = "marginal", data = surface,
#      marginal_covariates = "covar1")
# 
# # Step 3: Add one spline term with small df first
# fit_spline <- terradish(
#   S ~ covar2 + s(covar1, df = 3),
#   data              = surface,
#   conductance_model = smooth_loglinear_conductance,
#   measurement_model = mlpe,
#   optimizer         = rec$optimizer,
#   control           = rec$control,
#   solver            = rec$solver,
#   solver_control    = rec$solver_control,
#   approximation     = rec$approximation,
#   approximation_control = rec$approximation_control
# )
# 
# # Step 4: Compare with AICc
# aic_table(list(fit_linear, fit_spline),
#           mod_names = c("Linear", "Spline df=3"),
#           AICc = TRUE)
# 
# # Step 5: Inspect the spline marginal association
# plot(fit_spline, type = "marginal", data = surface)
# 
# # Step 6: If AICc improved and the curve is plausible, try more flexibility
# fit_spline4 <- terradish(
#   S ~ covar2 + s(covar1, df = 4),
#   data              = surface,
#   conductance_model = smooth_loglinear_conductance,
#   measurement_model = mlpe,
#   optimizer         = rec$optimizer,
#   control           = rec$control,
#   solver            = rec$solver,
#   solver_control    = rec$solver_control,
#   approximation     = rec$approximation,
#   approximation_control = rec$approximation_control
# )
# 
# aic_table(list(fit_linear, fit_spline, fit_spline4),
#           mod_names = c("Linear", "df=3", "df=4"),
#           AICc = TRUE)
# 
# # Step 7: Visualize the selected model
# best_fit <- fit_spline  # or whichever won
# plot(best_fit, type = "fit",             data = surface)
# plot(best_fit, type = "surface",         data = surface)
# plot(best_fit, type = "marginal",        data = surface)
# plot(best_fit, type = "marginal_response", data = surface)

## ----overfit-check, eval = FALSE----------------------------------------------
# # Quick eigenvalue check: are estimated parameters degenerate?
# hess <- fit_ns_alt$fit$hessian
# ev <- eigen(hess, symmetric = TRUE, only.values = TRUE)$values
# min(ev)   # Should be clearly positive; near-zero → numerical problems
# 
# # Compare ns and bs marginal shapes side by side
# par(mfrow = c(1, 2))
# plot(fit_ns_alt, type = "marginal", data = surface, main = "ns")
# plot(fit_bs_alt, type = "marginal", data = surface, main = "bs")

## ----recovery-path------------------------------------------------------------
system.file("examples", "spline-recovery-melip.R", package = "terradish")

## ----recovery-script, eval = FALSE--------------------------------------------
# source(system.file("examples", "spline-recovery-melip.R",
#                    package = "terradish"))
# 
# # Simulate a hump-shaped altitude response and check recovery
# recovery <- run_spline_recovery_example()
# recovery$model_comparison   # AICc table
# recovery$parameter_recovery # true vs estimated theta
# 
# # Full melip worked example
# melip_spline <- fit_melip_spline_example()
# melip_spline$model_comparison

## ----quickref, eval = FALSE---------------------------------------------------
# library(terradish)
# library(terra)
# 
# # ---- 1. Prepare the landscape ----
# data(melip)
# melip.altitude    <- terra::unwrap(melip.altitude)
# melip.forestcover <- terra::unwrap(melip.forestcover)
# melip.coords      <- terra::unwrap(melip.coords)
# 
# covariates <- c(melip.altitude, melip.forestcover)
# names(covariates) <- c("altitude", "forestcover")
# covariates <- scale_covariates(covariates)   # Always scale before fitting
# surface <- conductance_surface(covariates, melip.coords,
#                                directions = 8, saveStack = TRUE)
# 
# # ---- 2. Assess once and reuse recommended settings ----
# assessment <- terradish_assess_settings(
#   melip.Fst ~ forestcover + s(altitude, df = 4),
#   data              = surface,
#   conductance_model = smooth_loglinear_conductance,
#   measurement_model = mlpe,
#   probe_maxit       = 1,
#   verbose           = FALSE
# )
# rec <- assessment$recommended
# 
# # ---- 3. Fit candidate models ----
# # Linear baseline
# fit_lin <- terradish(
#   melip.Fst ~ altitude + forestcover,
#   data = surface, conductance_model = loglinear_conductance,
#   measurement_model = mlpe,
#   optimizer = rec$optimizer,
#   control = rec$control,
#   solver = rec$solver,
#   solver_control = rec$solver_control,
#   approximation = rec$approximation,
#   approximation_control = rec$approximation_control
# )
# 
# # Natural spline on altitude (df = 4)
# fit_ns <- terradish(
#   melip.Fst ~ forestcover + s(altitude, df = 4),
#   data = surface, conductance_model = smooth_loglinear_conductance,
#   measurement_model = mlpe,
#   optimizer = rec$optimizer,
#   control = rec$control,
#   solver = rec$solver,
#   solver_control = rec$solver_control,
#   approximation = rec$approximation,
#   approximation_control = rec$approximation_control
# )
# 
# # B-spline on altitude (k = 4, cubic)
# fit_bs <- terradish(
#   melip.Fst ~ forestcover + s(altitude, k = 4, basis = "bs", degree = 3),
#   data = surface, conductance_model = smooth_loglinear_conductance,
#   measurement_model = mlpe,
#   optimizer = rec$optimizer,
#   control = rec$control,
#   solver = rec$solver,
#   solver_control = rec$solver_control,
#   approximation = rec$approximation,
#   approximation_control = rec$approximation_control
# )
# 
# # ---- 4. Compare ----
# aic_table(list(fit_lin, fit_ns, fit_bs),
#           mod_names = c("Linear", "ns df=4", "bs k=4"),
#           AICc = TRUE)
# 
# # ---- 5. Visualize the best model ----
# best <- fit_ns   # replace with whichever model won
# 
# plot(best, type = "fit")                       # Observed vs fitted
# plot(best, type = "surface",  data = surface)  # Conductance raster + CI
# 
# # Marginal plots: all covariates together, or one at a time
# plot(best, type = "marginal", data = surface)
# plot(best, type = "marginal", data = surface,
#      marginal_covariates = "altitude")          # altitude panel only
# 
# # Side-by-side comparison of two models' altitude associations. patchwork is
# # optional; install it for the `|` and `/` layout operators used here.
# library(patchwork)
# p1 <- plot(fit_lin,  type = "marginal", data = surface,
#            marginal_covariates = "altitude") + ggtitle("Linear")
# p2 <- plot(fit_ns,   type = "marginal", data = surface,
#            marginal_covariates = "altitude") + ggtitle("ns df=4")
# p1 | p2
# 
# plot(best, type = "marginal_response",         # Response-scale associations
#      data = surface, marginal_covariates = "altitude")
# 
# # Extract the conductance surface as a SpatRaster
# cond_raster <- conductance(surface, best)
# # Layer names: "est", "lower95", "upper95"
# plot(cond_raster)

