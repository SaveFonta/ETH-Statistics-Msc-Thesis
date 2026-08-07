# =============================================================================
# STEP 1 - Calibrate the designs to a common avgPower (and, for SC, avgT1E)
# =============================================================================
# Every (looks x criterion) combination is calibrated so that avgPower =
# target_power at delta_MCID (via the sample size). For SC, avgT1E is also
# pinned to target_t1e (via the common threshold p), same as before.
#
# DC cannot be pinned to the same avgT1E target: its second condition
# Pr(delta > DV | data) > p_DV is evaluated under the null (delta = 0), and as
# n grows the posterior concentrates around that true null, so the chance of
# it spuriously wandering past DV collapses toward 0 - checked directly, the
# avgT1E ceiling (at the most liberal p) is already below target_t1e = 0.05
# by n2 = 20, and keeps shrinking as n grows further. So for DC, p is just
# fixed at the most liberal boundary (p_interval[1]) and only n is searched
# for power; avgT1E is left to be whatever that naturally gives, and is
# expected to land well below SC's.
#
#   calibrate_threshold        : given a schedule, find p such that avgT1E = target (SC).
#   calibrate_design           : jointly find n and p for SC (T1E + power).
#   calibrate_design_natural_t1e : same outer n-search as calibrate_design, but
#                                  p fixed, no T1E calibration (DC).
#
# The calibrated designs are saved to Output/07_calibrate_designs/designs.RDS.
# The EUII comparison across delta is done later, in 08_calibrated_euii_comparison.R.
# =============================================================================

library(dplyr)
source("scripts/00_shared_setup.R")          # sigma, p_MAP, p_vague, delta_MCID, DV
source("scripts/00_functions.R")

out_dir <- "Output/07_calibrate_designs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# DC helper: same outer (bracket + bisect) sample-size search as
# calibrate_design, but p is fixed at p_interval[1] throughout, no inner
# calibrate_threshold call.
# ---------------------------------------------------------------------------
calibrate_design_natural_t1e <- function(target_power, delta_power, n1_base, n2_base,
                                         prior_cntrl, prior_treat, design_prior = prior_cntrl,
                                         criterion, DV, p_DV = 0.5, p_interval = c(0.6, 0.999),
                                         sigma_1 = 88, sigma_2 = 88,
                                         n_sim = 1e5, seed = 123,
                                         n2_max_cap = 64 * max(n2_base)) {
  allo_ratio <- n1_base[length(n1_base)] / n2_base[length(n2_base)]
  info_ratio <- n2_base / n2_base[length(n2_base)]
  p_fixed <- p_interval[1]

  schedule <- function(n_2max) {
    n2 <- pmax(1, round(info_ratio * n_2max))
    n2 <- cummax(n2)
    if (any(diff(n2) <= 0)) n2 <- n2[1] + cumsum(c(0, pmax(1, diff(n2))))
    n1 <- pmax(1, round(allo_ratio * n2))
    list(n1 = n1, n2 = n2)
  }

  eval_at <- function(n_2max) {
    ns <- schedule(n_2max)
    crit <- if (criterion == "SC") {
      RBesT::decision2S(pc = p_fixed, qc = 0, lower.tail = TRUE)
    } else {
      RBesT::decision2S(pc = c(p_fixed, p_DV), qc = c(0, -DV), lower.tail = TRUE)
    }
    decisions <- lapply(seq_along(ns$n1), function(k) list(success = crit, futility = NULL))
    r <- avgoc2_seq_mc.normMix(
      prior_1 = prior_treat, prior_2 = prior_cntrl,
      n1_seq = ns$n1, n2_seq = ns$n2, decisions_list = decisions,
      delta = c(0, delta_power), design_prior_c = design_prior,
      sigma_1 = sigma_1, sigma_2 = sigma_2, n_sim = n_sim, seed = seed
    )
    list(power = r[[paste0("delta.", delta_power)]]$Overall[["Power"]],
         t1e   = r[["delta.0"]]$Overall[["Power"]],
         ns = ns, decisions = decisions)
  }

  n_max_curr <- n2_base[length(n2_base)]
  e_curr <- eval_at(n_max_curr)

  if (e_curr$power >= target_power) {
    hi <- n_max_curr; e_hi <- e_curr
    lo <- max(2, floor(n_max_curr / 2)); e_lo <- eval_at(lo)
    while (e_lo$power >= target_power && lo > 2) {
      hi <- lo; e_hi <- e_lo
      lo <- max(2, floor(lo / 2)); e_lo <- eval_at(lo)
    }
    if (e_lo$power >= target_power) {
      stop("Target power is reached even at the minimum sample size of 2.")
    }
  } else {
    lo <- n_max_curr
    hi <- 2 * n_max_curr; e_hi <- eval_at(hi)
    while (e_hi$power < target_power) {
      lo <- hi; hi <- 2 * hi
      if (hi > n2_max_cap) stop("power target not reached below n2_max_cap = ", n2_max_cap)
      e_hi <- eval_at(hi)
    }
  }

  while (hi - lo > 1) {
    mid <- floor((lo + hi) / 2)
    e_mid <- eval_at(mid)
    if (e_mid$power >= target_power) { hi <- mid; e_hi <- e_mid } else lo <- mid
  }

  list(n1_seq = e_hi$ns$n1, n2_seq = e_hi$ns$n2, p = p_fixed,
       achieved_t1e = e_hi$t1e, achieved_power = e_hi$power,
       decisions = e_hi$decisions)
}


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
target_t1e   <- 0.05
target_power <- 0.90

n_sim_calib <- 5e4    # sims per evaluation inside the calibration search
n_sim_check <- 2e5    # sims for the final verification of each design
mc_seed     <- 123

base_schedules <- list(
  "Fixed"        = list(n1 = 40,            n2 = 20),
  "One interim"  = list(n1 = c(20, 40),     n2 = c(10, 20)),
  "Two interims" = list(n1 = c(13, 27, 40), n2 = c(7, 13, 20))
)
criteria <- c("SC", "DC")
cores <- 6   # each (looks x criterion) calibration is independent, run them in parallel


# ---------------------------------------------------------------------------
# Joint calibration, every looks x criterion combination
# ---------------------------------------------------------------------------
# Same avgPower at delta_MCID for all 3 x 2 = 6 combinations; same avgT1E too,
# but only within SC (DC's avgT1E is left at its natural, much lower level).
cat("=============================================================\n")
cat("Joint (n, p) calibration, every looks x criterion combination\n")
cat("=============================================================\n")

combos <- expand.grid(Looks = names(base_schedules), Criterion = criteria,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

if (.Platform$OS.type != "unix") cores <- 1L   # mclapply cannot fork on Windows
cat(nrow(combos), "combinations, cores =", cores, "\n")

designs_cal <- parallel::mclapply(seq_len(nrow(combos)), function(i) {
  looks <- combos$Looks[i]; crit <- combos$Criterion[i]
  bs <- base_schedules[[looks]]
  if (crit == "SC") {
    calibrate_design(target_t1e = target_t1e, target_power = target_power,
                     delta_power = delta_MCID,
                     n1_base = bs$n1, n2_base = bs$n2,
                     prior_cntrl = p_MAP, prior_treat = p_vague,
                     criterion = crit, DV = DV,
                     n_sim = n_sim_calib, seed = mc_seed)
  } else {
    calibrate_design_natural_t1e(target_power = target_power, delta_power = delta_MCID,
                                 n1_base = bs$n1, n2_base = bs$n2,
                                 prior_cntrl = p_MAP, prior_treat = p_vague,
                                 criterion = crit, DV = DV,
                                 n_sim = n_sim_calib, seed = mc_seed)
  }
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
