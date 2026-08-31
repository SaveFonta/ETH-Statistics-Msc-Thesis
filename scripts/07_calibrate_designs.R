# =============================================================================
# Calibrate the SC designs to a common avgT1E and avgPower
# =============================================================================
# avgT1E = target_t1e (via the common threshold p) and avgPower = target_power
# at delta_MCID (via the sample size). Only for SC.
#
#   calibrate_threshold : given a schedule, find p such that avgT1E = target.
#   calibrate_design    : jointly find n and p (T1E + power).
#
# The calibrated designs are saved to Output/07_calibrate_designs/designs.RDS.
# The EUII comparison across delta is done later, in 08_calibrated_euii_comparison.R.
# =============================================================================

source("scripts/00_shared_setup.R")          # sigma, p_MAP, p_vague, delta_MCID, DV
source("scripts/00_functions.R")             

out_dir <- "Output/07_calibrate_designs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
target_t1e   <- 0.05
target_power <- 0.90

n_sim_calib <- 1e6    # sims per evaluation inside the calibration search
n_sim_check <- 2e6    # sims for the final verification of each design
mc_seed     <- 123

base_schedules <- list(
  "Fixed"        = list(n1 = 40,            n2 = 20),
  "One interim"  = list(n1 = c(20, 40),     n2 = c(10, 20)),
  "Two interims" = list(n1 = c(13, 27, 40), n2 = c(7, 13, 20))
)
cores <- length(base_schedules)  


# ---------------------------------------------------------------------------
# Joint calibration, one per number of looks
# ---------------------------------------------------------------------------
cat("===== Start Calibration ====")

combos <- data.frame(Looks = names(base_schedules), Criterion = "SC",
                     stringsAsFactors = FALSE)

if (.Platform$OS.type != "unix") cores <- 1L   # mclapply cannot fork on Windows
cat(nrow(combos), "combinations, cores =", cores, "\n")

designs_cal <- parallel::mclapply(combos$Looks, function(looks) {
  bs <- base_schedules[[looks]]
  calibrate_design(target_t1e = target_t1e, target_power = target_power,
                   delta_power = delta_MCID,
                   n1_base = bs$n1, n2_base = bs$n2,
                   prior_cntrl = p_MAP, prior_treat = p_vague,
                   criterion = "SC", DV = DV,
                   n_sim = n_sim_calib, seed = mc_seed)
}, mc.cores = cores)
names(designs_cal) <- paste(combos$Looks, combos$Criterion, sep = " | ")

bad <- vapply(designs_cal, inherits, logical(1), "try-error")
if (any(bad)) {
  cat("Calibration failed for:\n")
  cat(paste0("   ", names(designs_cal)[bad], "\n"), sep = "")
  stop("calibrate_design failed for ", sum(bad), " of ", nrow(combos), " combinations.")
}


# ---------------------------------------------------------------------------
# Verification at higher precision, fresh seed
# ---------------------------------------------------------------------------
cat("\n=============================================================\n")
cat("Verification of the calibrated designs (fresh seed, n_sim =", n_sim_check, ")\n")
cat("=============================================================\n\n")
cat(" design                    n_max   p        avgT1E   avgPower  E[N](", delta_MCID, ")\n")

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
  cat(" ", format(nm, width = 24),
      format(max(d$n1_seq) + max(d$n2_seq), width = 5),
      round(d$p, 4), " ",
      round(r[["delta.0"]]$Overall[["Power"]], 4), "  ",
      round(ov[["Power"]], 4), "   ",
      round(ov[["EN_t"]] + ov[["EN_c"]], 1), "\n")

  # The search picks the smallest n whose *estimated* power clears the target
  # during calibration (n_sim_calib); this fresh, higher-precision check is
  # what actually confirms it. Flag it if that confirmation doesn't hold up -
  # either the search got lucky off calibration noise (below target here), or
  # it's close enough to target that it plausibly could have (within ~2 MC
  # standard errors of this fresh estimate).
  se <- sqrt(ov[["Power"]] * (1 - ov[["Power"]]) / n_sim_check)
  if (ov[["Power"]] < target_power) {
    warning(nm, ": re-verified avgPower (", round(ov[["Power"]], 4),
            ") is BELOW target_power (", target_power,
            ") - the calibration search likely accepted calibration-noise luck.")
  } else if (ov[["Power"]] - target_power < 2 * se) {
    warning(nm, ": re-verified avgPower (", round(ov[["Power"]], 4),
            ") is within 2 MC standard errors of target_power (", target_power,
            ") - margin is thin, consider a larger n_sim_calib.")
  }
}

saveRDS(list(designs_cal   = designs_cal,
             combos        = combos,
             target_t1e    = target_t1e,
             target_power  = target_power,
             delta_MCID    = delta_MCID,
             n_sim_calib   = n_sim_calib,
             mc_seed       = mc_seed),
        file = file.path(out_dir, "designs.RDS"))
cat("\nCalibrated designs saved to", file.path(out_dir, "designs.RDS"), "\n")
