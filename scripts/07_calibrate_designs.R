# =============================================================================
# STEP 1 - Calibrate the designs to a common avgPower and avgT1E
# =============================================================================
# Every design (fixed, one interim, two interims) is calibrated so that
#   avgT1E   = target_t1e                (via the common threshold p)
#   avgPower = target_power at delta_MCID (via the sample size)
# both under the aligned prior. After this step the designs are exchangeable on
# their error rates: they differ only in the stopping structure and in the
# sample size that structure requires.
#
#   calibrate_threshold : given a schedule, find p such that avgT1E = target.
#   calibrate_design    : jointly find the smallest schedule (proportional
#                         scaling) and p such that both targets hold. The
#                         threshold is recalibrated at every candidate size, so
#                         the T1E target holds exactly at the returned n.
#
# The calibrated designs are saved to Output/07_calibrate_designs/designs.RDS.
# The EUII comparison across delta is done later, in 08_calibrated_euii_comparison.R.
# =============================================================================

library(dplyr)
source("scripts/00_shared_setup.R")          # sigma, p_MAP, p_vague, delta_MCID
source("scripts/00_functions.R")

out_dir <- "Output/07_calibrate_designs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
target_t1e   <- 0.05
target_power <- 0.90

n_sim_calib <- 5e4    # sims per evaluation inside the calibration search
n_sim_check <- 2e5    # sims for the final verification of each design
mc_seed     <- 123


# ---------------------------------------------------------------------------
# Threshold calibration alone, at a fixed schedule
# ---------------------------------------------------------------------------
# Two looks at half the sample. The classical Pocock threshold for two equally
# spaced looks is 1 - 0.0294 = 0.9706 without borrowing. With borrowing the
# interim carries a higher effective information fraction, the looks are more
# correlated, and slightly less correction is needed.
cat("=============================================================\n")
cat("Threshold calibration at a fixed schedule (demo)\n")
cat("=============================================================\n")

cal_p <- calibrate_threshold(target_t1e,
                             n1_seq = c(20, 40), n2_seq = c(10, 20),
                             prior_cntrl = p_MAP, prior_treat = p_vague,
                             criterion = "SC", n_sim = n_sim_calib, seed = mc_seed)
cat("  calibrated p    :", round(cal_p$p, 5), "\n")
cat("  achieved avgT1E :", round(cal_p$achieved_t1e, 5), "\n")
cat("  Pocock (vague, two looks) would be 0.9706\n\n")



cal_p2 <- calibrate_threshold(target_t1e,
                             n1_seq = 60, n2_seq = 30,
                             prior_cntrl = p_vague, prior_treat = p_vague,
                             criterion = "SC", n_sim = n_sim_calib, seed = mc_seed)


# ---------------------------------------------------------------------------
# Joint calibration per design: same avgT1E, same avgPower at delta_MCID
# ---------------------------------------------------------------------------
# The base schedules define the look fractions and the 2:1 allocation.
cat("=============================================================\n")
cat("Joint (n, p) calibration per design\n")
cat("=============================================================\n")

base_schedules <- list(
  "Fixed"        = list(n1 = 40,            n2 = 20),
  "One interim"  = list(n1 = c(20, 40),     n2 = c(10, 20)),
  "Two interims" = list(n1 = c(13, 27, 40), n2 = c(7, 13, 20))
)

designs_cal <- lapply(names(base_schedules), function(nm) {
  cat("\n---", nm, "---\n")
  bs <- base_schedules[[nm]]
  calibrate_design(target_t1e = target_t1e, target_power = target_power,
                   delta_power = delta_MCID,
                   n1_base = bs$n1, n2_base = bs$n2,
                   prior_cntrl = p_MAP, prior_treat = p_vague,
                   criterion = "SC",
                   n_sim = n_sim_calib, seed = mc_seed)
})
names(designs_cal) <- names(base_schedules)


# ---------------------------------------------------------------------------
# Verification at higher precision, fresh seed
# ---------------------------------------------------------------------------
cat("\n=============================================================\n")
cat("Verification of the calibrated designs (fresh seed, n_sim =", n_sim_check, ")\n")
cat("=============================================================\n\n")
cat(" design        n_max   p        avgT1E   avgPower  E[N](", delta_MCID, ")\n")

for (nm in names(designs_cal)) {
  d <- designs_cal[[nm]]
  r <- avgoc2_seq_mc.normMix(
    prior_1        = p_vague,
    prior_2        = p_MAP,
    n1_seq         = d$n1_seq,
    n2_seq         = d$n2_seq,
    decisions_list = d$decisions,
    delta          = c(0, delta_MCID),
    design_prior_c = p_MAP,
    n_sim          = n_sim_check,
    seed           = mc_seed + 1     # fresh seed: honest check of the calibration
  )
  ov <- r[[paste0("delta.", delta_MCID)]]$Overall
  cat(" ", format(nm, width = 12),
      format(max(d$n1_seq) + max(d$n2_seq), width = 5),
      round(d$p, 4), " ",
      round(r[["delta.0"]]$Overall[["Power"]], 4), "  ",
      round(ov[["Power"]], 4), "   ",
      round(ov[["EN_t"]] + ov[["EN_c"]], 1), "\n")
}

saveRDS(list(designs_cal   = designs_cal,
             target_t1e    = target_t1e,
             target_power  = target_power,
             delta_MCID    = delta_MCID,
             n_sim_calib   = n_sim_calib,
             mc_seed       = mc_seed),
        file = file.path(out_dir, "designs.RDS"))
cat("\nCalibrated designs saved to", file.path(out_dir, "designs.RDS"), "\n")


