# =============================================================================
# STEP 1 - Calibrate different boundary shapes to a common avgPower and avgT1E
# =============================================================================
# Fixes the number of looks (one interim + final, K = 2) and instead varies
# the *shape* of the boundary across that fixed schedule: Flat (Pocock-style
# constant threshold, calibrate_design's default), O'Brien-Fleming,
# Haybittle-Peto, a power-family spending shape, and the Shi & Yin (2019)
# interpolation style. All builders live in 00_functions_examples.R.
#
# Every shape is calibrated so that
#   avgT1E   = target_t1e                (via the shape's own parameter)
#   avgPower = target_power at delta_MCID (via the sample size)
# so the shapes are exchangeable on their error rates: they differ only in
# how they spend that error across the two looks, and in the sample size that
# choice requires.
#
# Holding K fixed here isolates the shape's effect: 07_calibrate_designs.R
# instead holds the shape fixed (flat) and varies K, so the two comparisons
# don't get confounded with one another.
#
# The calibrated designs are saved to Output/09_calibrate_boundary_shapes/designs.RDS.
# The EUII comparison across delta is done later, in 10_boundary_shapes_euii_comparison.R.
# =============================================================================

library(dplyr)
source("scripts/00_shared_setup.R")            # sigma, p_MAP, p_vague, delta_MCID, DV
source("scripts/00_functions_examples.R")      # sources 00_functions.R too, plus the builders

out_dir <- "Output/09_calibrate_boundary_shapes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
target_t1e   <- 0.05
target_power <- 0.90

n_sim_calib <- 5e4    # sims per evaluation inside the calibration search
n_sim_check <- 2e5    # sims for the final verification of each design
mc_seed     <- 123
cores       <- 6      # each shape's calibration is independent, run them in parallel

base_schedule <- list(n1 = c(20, 40), n2 = c(10, 20))   # one interim, same as 07's "One interim"

# Shi & Yin's parameter u lives on a different scale than the other builders'
# posterior probability p, see 00_functions_examples.R section 7. Two
# adjustments on top of the generic builder_sy_1int (which assumes a raw
# nominal alpha = 0.025, half our target_t1e):
#   - build the nominal boundary at alpha = target_t1e directly, so u = 0
#     already lands close to target instead of needlessly overshooting;
#   - even so, the nominal frequentist Z-boundary under-covers once translated
#     into a posterior-probability threshold applied to a design that borrows
#     via an informative (MAP) control prior (the "aligned prior" T1E result
#     is a single-look, non-borrowing-adjusted identity, see the calibration
#     discussion), so u = 0 is still too strict here; the search interval
#     must allow u < 0 (loosening below the nominal boundary) to attain
#     target_t1e = 0.05, and must stop short of u = 1 (threshold = 1 is a
#     degenerate, unsolvable decision boundary).
builder_sy_matched <- make_builder_sy_obf(c(1 / 2, 1), alpha = target_t1e)

shapes <- list(
  "Flat"           = list(builder = NULL,               p_interval = c(0.6, 0.999)),
  "OBF"            = list(builder = builder_obf_1int,    p_interval = c(0.6, 0.999)),
  "Haybittle-Peto" = list(builder = builder_hp_1int,     p_interval = c(0.6, 0.999)),
  "Power-family"   = list(builder = builder_pow_1int,    p_interval = c(0.6, 0.999)),
  "Shi & Yin"      = list(builder = builder_sy_matched,  p_interval = c(-5, 0.999))
)
shape_names <- names(shapes)


# ---------------------------------------------------------------------------
# Joint calibration, every shape
# ---------------------------------------------------------------------------
cat("=============================================================\n")
cat("Joint (n, p) calibration, every boundary shape, K = 2 fixed\n")
cat("=============================================================\n")
cat(length(shapes), "shapes, cores =", if (.Platform$OS.type != "unix") 1L else cores, "\n")

if (.Platform$OS.type != "unix") cores <- 1L   # mclapply cannot fork on Windows

designs_cal <- parallel::mclapply(shape_names, function(nm) {
  sh <- shapes[[nm]]
  calibrate_design(target_t1e = target_t1e, target_power = target_power,
                   delta_power = delta_MCID,
                   n1_base = base_schedule$n1, n2_base = base_schedule$n2,
                   prior_cntrl = p_MAP, prior_treat = p_vague,
                   criterion = "SC", p_interval = sh$p_interval,
                   decisions_builder = sh$builder,
                   n_sim = n_sim_calib, seed = mc_seed)
}, mc.cores = cores)
names(designs_cal) <- shape_names

bad <- vapply(designs_cal, inherits, logical(1), "try-error")
if (any(bad)) {
  cat("Calibration failed for:\n")
  cat(paste0("   ", shape_names[bad], "\n"), sep = "")
  stop("calibrate_design failed for ", sum(bad), " of ", length(shapes), " shapes.")
}


# ---------------------------------------------------------------------------
# Verification at higher precision, fresh seed
# ---------------------------------------------------------------------------
cat("\n=============================================================\n")
cat("Verification of the calibrated designs (fresh seed, n_sim =", n_sim_check, ")\n")
cat("=============================================================\n\n")
cat(" shape             n_max   p        avgT1E   avgPower  E[N](", delta_MCID, ")\n")

for (nm in shape_names) {
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
  cat(" ", format(nm, width = 16),
      format(max(d$n1_seq) + max(d$n2_seq), width = 5),
      round(d$p, 4), " ",
      round(r[["delta.0"]]$Overall[["Power"]], 4), "  ",
      round(ov[["Power"]], 4), "   ",
      round(ov[["EN_t"]] + ov[["EN_c"]], 1), "\n")

  # Same check as 07_calibrate_designs.R: the search picks the smallest n
  # whose *estimated* power clears target during calibration; this fresh,
  # higher-precision pass is what actually confirms it.
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
             shape_names   = shape_names,
             target_t1e    = target_t1e,
             target_power  = target_power,
             delta_MCID    = delta_MCID,
             n_sim_calib   = n_sim_calib,
             mc_seed       = mc_seed),
        file = file.path(out_dir, "designs.RDS"))
cat("\nCalibrated designs saved to", file.path(out_dir, "designs.RDS"), "\n")
