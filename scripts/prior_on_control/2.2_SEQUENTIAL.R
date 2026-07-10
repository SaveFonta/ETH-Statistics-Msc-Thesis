# =============================================================================
# Sequential Design — verification of the Monte Carlo engine
# =============================================================================
# This file produces NO thesis results. Its only purpose is to check that the
# Monte Carlo engine reproduces the exact analytical sequential operating
# characteristics, so that 2.3_SEQUENTIAL.MC.R can be trusted.
#
#   reference  : oc2S_seq.dual.normMix  (00.functions.R, exact, conditional)
#   under test : oc2_seq_mc.normMix     (00.functions.MC.R)
#
# Only the CONDITIONAL operating characteristics are compared. That is enough:
# the averaged metrics are obtained by drawing theta_C from the design prior,
# which is exact Monte Carlo integration by construction, so once the conditional
# engine is right the averaged one is too.
#
# Each analytical evaluation costs a few seconds, so the grid is deliberately
# small. Expect a couple of minutes.
#
# Run from the Thesis/ root:
#   Rscript scripts/prior_on_control/2.2_SEQUENTIAL.R
# =============================================================================

suppressMessages({
  library(RBesT)
  library(dplyr)
})
source("scripts/prior_on_control/00.functions.R")      # oc2S_seq.dual.normMix
source("scripts/prior_on_control/00.functions.MC.R")   # oc2_seq_mc.normMix


# ---------------------------------------------------------------------------
# Design parameters (identical to 2.3_SEQUENTIAL.MC.R)
# ---------------------------------------------------------------------------
sigma      <- 88
delta_MCID <- 60

n1_seq <- c(20, 40)   # treatment arm, cumulative
n2_seq <- c(10, 20)   # control arm, cumulative
n_stg1 <- n1_seq[1] + n2_seq[1]
n_stg2 <- n1_seq[2] + n2_seq[2]

n_sim <- 200000       # MC replicates per evaluation


# ---------------------------------------------------------------------------
# Priors
# ---------------------------------------------------------------------------
p_MAP <- mixnorm(
  c(0.4848, -52.457, 21.154), c(0.4598, -47.465, 7.843), c(0.0554, -50.355, 48.164),
  sigma = sigma, param = "ms"
)
p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")


# ---------------------------------------------------------------------------
# Decision criteria (identical to 2.3)
# ---------------------------------------------------------------------------
# Efficacy, dual criterion: Pr(delta > 0) > 0.95 AND Pr(delta > 50) > 0.50
sign.crit <- decision2S(pc = c(0.95, 0.50), qc = c(0, -50), lower.tail = TRUE)
# Futility, stage 1 only: Pr(delta < 40 | data) >= 0.90
fut.crit  <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)

# The analytical function accepts a bare decision2S for the final stage.
decisions_an <- list(list(success = sign.crit, futility = fut.crit), sign.crit)
# The MC engine wants the success/futility slots at every stage.
decisions_mc <- list(list(success = sign.crit, futility = fut.crit),
                     list(success = sign.crit, futility = NULL))


# ---------------------------------------------------------------------------
# Tolerances: everything is compared against 4 Monte Carlo standard errors.
# N takes only the two values n_stg1 and n_stg2, so its variance is Bernoulli.
# ---------------------------------------------------------------------------
se_prop <- function(p) sqrt(p * (1 - p) / n_sim)
se_EN   <- function(p_stop1) sqrt(p_stop1 * (1 - p_stop1) / n_sim) * abs(n_stg2 - n_stg1)

report <- function(label, analytic, mc, tol) {
  d  <- mc - analytic
  ok <- abs(d) < tol
  cat("  ", formatC(label, width = -18),
      " analytic =", formatC(analytic, format = "f", digits = 5, width = 9),
      "  MC =",      formatC(mc,       format = "f", digits = 5, width = 9),
      "  diff =",    formatC(d,        format = "f", digits = 5, width = 9),
      " ", if (ok) "PASS" else "**FAIL**", "\n")
  ok
}


# ---------------------------------------------------------------------------
# Conditional OC: analytical versus MC
# ---------------------------------------------------------------------------
cat("\n===============================================================\n")
cat("Conditional OC, analytical versus MC (n_sim =", n_sim, ")\n")
cat("tolerance: 4 Monte Carlo standard errors\n")
cat("===============================================================\n")

theta_c_check   <- c(-70, -50, -30)
analysis_priors <- list("MAP" = p_MAP, "Vague" = p_vague)
all_ok <- TRUE

for (pname in names(analysis_priors)) {

  oc_an <- oc2S_seq.dual.normMix(
    prior1 = p_vague, prior2 = analysis_priors[[pname]],
    n1 = n1_seq, n2 = n2_seq, decisions = decisions_an
  )

  for (tc in theta_c_check) {
    # delta = 0 gives the type I error, delta = delta_MCID gives the power.
    for (d in c(0, delta_MCID)) {

      a <- oc_an(tc - d, tc)                       # theta_T = theta_C - delta
      m <- oc2_seq_mc.normMix(
        theta_1 = tc - d, theta_2 = tc,
        prior_1 = p_vague, prior_2 = analysis_priors[[pname]],
        n1_seq  = n1_seq,  n2_seq  = n2_seq,
        decisions_list = decisions_mc, n_sim = n_sim, seed = 1
      )
      ps <- m$Per_Stage
      p_stop1 <- a$P_Stop_Eff_Stg1 + a$P_Stop_Fut_Stg1

      cat("\n", pname, " | theta_C =", tc, " | delta =", d, "\n")
      all_ok <- report("P_Stop_Eff_Stg1", a$P_Stop_Eff_Stg1, ps$P_Succ[1],
                       4 * se_prop(a$P_Stop_Eff_Stg1)) && all_ok
      all_ok <- report("P_Stop_Fut_Stg1", a$P_Stop_Fut_Stg1, ps$P_Fut[1],
                       4 * se_prop(a$P_Stop_Fut_Stg1)) && all_ok
      all_ok <- report("Total_Power",     a$Total_Power,     m$Overall[["Power"]],
                       4 * se_prop(a$Total_Power)) && all_ok
      all_ok <- report("EN (Trt + Ctrl)", a$EN_Trt + a$EN_Pbo,
                       m$Overall[["EN_t"]] + m$Overall[["EN_c"]],
                       4 * se_EN(p_stop1)) && all_ok
    }
  }
}

cat("\n\n===============================================================\n")
if (all_ok) {
  cat("ALL CHECKS PASSED: the MC engine matches the analytical reference.\n")
} else {
  cat("SOME CHECKS FAILED: inspect the lines marked **FAIL** above.\n")
}
cat("===============================================================\n")
