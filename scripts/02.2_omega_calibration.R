## ---------------------------------------------------------------------------
## Stress testing the robustness weight omega, rather than fixing it by
## convention, against a design prior that turns out to be wrong.
##
## The threshold is fixed at p = 0.95 throughout, matching sign.crit and every
## other design in this thesis. Under the design prior pi_MAP used for the
## sample size of Section 3.1.1, omega = 0 sits exactly at nominal (avgT1E =
## 0.05, the exact control result applied to a genuinely aligned design), and
## every omega > 0 sits strictly below it: a robust analysis prior is wider
## than pi_MAP, which inflates the posterior variance and so raises the
## decision threshold k = q + z_p sqrt(V_delta), making success harder to
## declare under the null exactly as it does under the alternative.
##
## Three ways to stress test omega against a design prior that turns out to be
## wrong are compared, mirroring the two types of prior conflict already
## isolated in the single Normal closed form of Figure 2.3 (location shift,
## scale shift), plus the discrete alternative of Best et al. (2025):
##   SWITCH  evaluate avgT1E under one of the named design priors of
##           Section 3.1 (skeptical, misspecified) instead of pi_MAP: a
##           different, discrete scenario, not a perturbation of pi_MAP.
##   SHIFT   slide pi_MAP itself by a location drift d: a continuous, one
##           parameter family that isolates a pure change in location,
##           holding shape and spread fixed. Attacks avgT1E through the bias
##           term (1-W)(mu_A,C - mu_D,C).
##   SCALE   stretch pi_MAP's spread by a factor s, mean and weights fixed.
##           Attacks avgT1E through the variance term (1-W)^2 sigma_D,C^2
##           instead of the bias term, so it is a structurally different
##           failure mode from SHIFT even when it produces similar avgT1E.
## ---------------------------------------------------------------------------

suppressMessages({
  library(RBesT)
  library(parallel)
})

source("scripts/00_functions.R")
source("scripts/00_shared_setup.R")

n_C <- 20
n_T <- 40
P_FIXED <- 0.95        # the threshold used everywhere else in the thesis
OUTDIR  <- "Output/02.2_omega_calibration"
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

## Shift a Normal mixture in location only: every component mean moves by d,
## the weights and the component spreads are untouched.
shift_mix <- function(mix, d) { mix[2, ] <- mix[2, ] + d; mix }

## Scale a Normal mixture's spread by a factor s: every component sd is
## multiplied by s, the means and weights are untouched.
scale_mix <- function(mix, s) { mix[3, ] <- mix[3, ] * s; mix }

## avgT1E rate for a given robustness weight and design prior, at the fixed
## threshold.
avg_t1e_at <- function(omega, pd, p = P_FIXED) {
  pa <- if (omega == 0) p_MAP else robustify(p_MAP, omega, mean = mu_A)
  crit <- decision2S(pc = p, qc = 0, lower.tail = TRUE)
  f <- avgoc2S.normMix(prior1 = p_vague, prior2 = pa, n1 = n_T, n2 = n_C,
                       decision = crit, delta = 0, design_prior2 = pd)
  as.numeric(f(0, pd))
}

skept <- design_priors[["skeptical (-90)"]]
missp <- design_priors[["misspecified (-10)"]]
sw_ws <- c(0, 0.2, 0.35, 0.5)
switch_shift <- do.call(rbind, mclapply(sw_ws, function(w) data.frame(
  omega        = w,
  t1e_map      = avg_t1e_at(w, p_MAP),
  switch_skept = avg_t1e_at(w, skept),
  switch_missp = avg_t1e_at(w, missp),
  shift_d15    = avg_t1e_at(w, shift_mix(p_MAP, -15)),
  shift_d30    = avg_t1e_at(w, shift_mix(p_MAP, -30)),
  scale_s2     = avg_t1e_at(w, scale_mix(p_MAP, 2)),
  scale_s3     = avg_t1e_at(w, scale_mix(p_MAP, 3))
), mc.cores = length(sw_ws)))
saveRDS(switch_shift, file.path(OUTDIR, "switch_vs_shift.RDS"))

## ---- report ---------------------------------------------------------------
cat("Fixed design, n_C = ", n_C, ", n_T = ", n_T, ", p = ", P_FIXED, "\n\n", sep = "")
print(switch_shift, row.names = FALSE)
