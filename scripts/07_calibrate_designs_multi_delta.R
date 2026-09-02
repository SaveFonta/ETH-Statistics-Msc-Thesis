# =============================================================================
# Multi delta calibration - number of interim analyses
# =============================================================================
# (n, p) is recalibrated at every delta_target in delta_targets, for each
# number of looks, and the EUII is read at the design's own calibration point.
# =============================================================================

library(dplyr)
library(ggplot2)
source("scripts/00_shared_setup.R")          # sigma, p_MAP, p_vague, DV, my_theme
source("scripts/00_functions.R")             # calibrate_design, compute_euii

out_dir <- "Output/07_calibrate_designs_multi_delta"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
target_t1e    <- 0.05
target_power  <- 0.90
delta_targets <- c(50, 60, 70, 80, 90)

n_sim_calib <- 1e6    # sims per evaluation inside the calibration search
n_sim_check <- 2e6    # sims for the verification of each calibrated design
mc_seed     <- 123
prior_H1    <- c(0.01, 0.1, 0.5)

# Only starting points for calibrate_design()'s bracket and bisect search over
# n, so the same three are reused at every delta_target: the search finds the
# right n on its own, whatever the target.
base_schedules <- list(
  "Fixed"        = list(n1 = 40,            n2 = 20),
  "One interim"  = list(n1 = c(20, 40),     n2 = c(10, 20)),
  "Two interims" = list(n1 = c(13, 27, 40), n2 = c(7, 13, 20))
)

combos <- expand.grid(Looks       = names(base_schedules),
                      DeltaTarget = delta_targets,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
combos$Label <- paste(combos$Looks, combos$DeltaTarget, sep = " | d")

cores <- min(32L, nrow(combos))
if (.Platform$OS.type != "unix") cores <- 1L   # mclapply not works on Windows

n_A <- as.numeric(ess(p_MAP))   # ~40.5

cat("===== Multi delta calibration - number of looks =====\n")
cat(nrow(combos), "combinations,", length(delta_targets), "targets x",
    length(base_schedules), "look counts, cores =", cores, "\n")
cat("targets:", paste(delta_targets, collapse = ", "), "\n")
cat("ESS of the analysis prior, n_A =", round(n_A, 1), "\n\n")


# ---------------------------------------------------------------------------
# Calibration, one cell per (look count x delta target)
# ---------------------------------------------------------------------------
designs_cal <- parallel::mclapply(seq_len(nrow(combos)), function(i) {
  bs <- base_schedules[[combos$Looks[i]]]
  calibrate_design(target_t1e = target_t1e, target_power = target_power,
                   delta_power = combos$DeltaTarget[i],
                   n1_base = bs$n1, n2_base = bs$n2,
                   prior_cntrl = p_MAP, prior_treat = p_vague,
                   criterion = "SC", DV = DV,
                   n_sim = n_sim_calib, seed = mc_seed)
}, mc.cores = cores)
names(designs_cal) <- combos$Label

bad <- vapply(designs_cal, inherits, logical(1), "try-error") # shouldn't happen
if (any(bad)) {
  cat("Calibration failed for:\n")
  cat(paste0("   ", combos$Label[bad], "\n"), sep = "")
  stop("calibrate_design failed for ", sum(bad), " of ", nrow(combos), " combinations.")
}


# ---------------------------------------------------------------------------
# Verification at new seed and higher precision, and the EUII at the
# design's own calibration point, both from the same simulation call
# ---------------------------------------------------------------------------
cat("===== Verification at a new seed, n_sim =", n_sim_check, "=====\n\n")

verified <- parallel::mclapply(seq_len(nrow(combos)), function(i) {
  d  <- designs_cal[[i]]
  dt <- combos$DeltaTarget[i]

  r <- avgoc2_seq_mc.normMix(
    prior_1        = p_vague,
    prior_2        = p_MAP,
    n1_seq         = d$n1_seq,
    n2_seq         = d$n2_seq,
    decisions_list = d$decisions,
    delta          = c(0, dt),
    design_prior_c = p_MAP,
    n_sim          = n_sim_check,
    seed           = mc_seed + 1     # fresh seed
  )

  ov  <- r[[paste0("delta.", dt)]]$Overall
  # EUII at every prior_H1 in the grid, so the plot can show the same
  # Pr(H1)-sensitivity band as the other EUII figures: EUII is the point
  # value at Pr(H1) = 0.5 (used throughout the prose and tables), EUII_lo/hi
  # the min/max across the full prior_H1 grid (used only for the ribbon).
  eu_all  <- compute_euii(r, prior_H1 = prior_H1)
  eu_vals <- vapply(eu_all, function(e) e$EUII[e$Delta == dt], numeric(1))
  # Effective sample size on each branch, across the full prior_H1 grid, for
  # the third panel of the EUII figures: (E[1/N | success])^-1 and
  # (E[1/N | failure])^-1. EN_success/EN_failure are the point values at
  # Pr(H1) = 0.5 (used throughout the prose and tables), the _lo/_hi the
  # min/max across the grid (used only for the ribbon).
  en_success_vals <- vapply(eu_all, function(e) 1 / e$E_invN_sig[e$Delta == dt], numeric(1))
  en_failure_vals <- vapply(eu_all, function(e) 1 / e$E_invN_nonsig[e$Delta == dt], numeric(1))

  list(
    n_max      = max(d$n1_seq) + max(d$n2_seq),
    n_C        = max(d$n2_seq),
    p          = d$p,
    avgT1E     = r[["delta.0"]]$Overall[["Power"]],
    avgPower   = ov[["Power"]],
    EN_success    = en_success_vals[["0.5"]],
    EN_success_lo = min(en_success_vals),
    EN_success_hi = max(en_success_vals),
    EN_failure    = en_failure_vals[["0.5"]],
    EN_failure_lo = min(en_failure_vals),
    EN_failure_hi = max(en_failure_vals),
    EN       = ov[["EN_t"]] + ov[["EN_c"]],
    EUII     = eu_vals[["0.5"]],
    EUII_lo  = min(eu_vals),
    EUII_hi  = max(eu_vals)
  )
}, mc.cores = cores)

bad_v <- vapply(verified, inherits, logical(1), "try-error")
if (any(bad_v)) {
  cat("Verification failed for:\n")
  cat(paste0("   ", combos$Label[bad_v], "\n"), sep = "")
  stop("verification failed for ", sum(bad_v), " of ", nrow(combos), " combinations.")
}


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
summary_df <- combos |>
  mutate(n_max    = vapply(verified, function(v) v$n_max,    numeric(1)),
         n_C      = vapply(verified, function(v) v$n_C,      numeric(1)),
         p        = vapply(verified, function(v) v$p,        numeric(1)),
         avgT1E   = vapply(verified, function(v) v$avgT1E,   numeric(1)),
         avgPower = vapply(verified, function(v) v$avgPower, numeric(1)),
         EN       = vapply(verified, function(v) v$EN,       numeric(1)),
         EN_success    = vapply(verified, function(v) v$EN_success,    numeric(1)),
         EN_success_lo = vapply(verified, function(v) v$EN_success_lo, numeric(1)),
         EN_success_hi = vapply(verified, function(v) v$EN_success_hi, numeric(1)),
         EN_failure    = vapply(verified, function(v) v$EN_failure,    numeric(1)),
         EN_failure_lo = vapply(verified, function(v) v$EN_failure_lo, numeric(1)),
         EN_failure_hi = vapply(verified, function(v) v$EN_failure_hi, numeric(1)),
         EUII     = vapply(verified, function(v) v$EUII,     numeric(1)),
         EUII_lo  = vapply(verified, function(v) v$EUII_lo,  numeric(1)),
         EUII_hi  = vapply(verified, function(v) v$EUII_hi,  numeric(1)),
         n_A_over_nC = n_A / n_C) |>
  group_by(DeltaTarget) |>
  mutate(EUII_ratio_to_fixed = EUII / EUII[Looks == "Fixed"]) |>
  ungroup() |>
  mutate(Looks = factor(Looks, levels = names(base_schedules))) |>
  arrange(DeltaTarget, Looks)


# ---------------------------------------------------------------------------
# Print, and flag any cell whose power margin is thin 
# ---------------------------------------------------------------------------
cat("Calibrated to avgT1E =", target_t1e, "and avgPower =", target_power,
    "at each design's own delta target.\n")
cat("EUII is quoted at that same target, at Pr(H1) = 0.5.\n\n")

print(as.data.frame(summary_df |>
        select(DeltaTarget, Looks, n_max, p, avgT1E, avgPower, EN, EUII,
               n_A_over_nC, EUII_ratio_to_fixed) |>
        mutate(across(c(p, avgT1E, avgPower, EUII, EUII_ratio_to_fixed),
                      \(x) round(x, 4)),
               across(c(EN, n_A_over_nC), \(x) round(x, 2)))),
      row.names = FALSE)

for (i in seq_len(nrow(summary_df))) {
  pw <- summary_df$avgPower[i]
  se <- sqrt(pw * (1 - pw) / n_sim_check)
  lb <- paste0(summary_df$Looks[i], " at delta = ", summary_df$DeltaTarget[i])
  if (pw < target_power) {
    warning(lb, ": re-verified avgPower (", round(pw, 4),
            ") is BELOW target_power (", target_power, ")")
  } else if (pw - target_power < 2 * se) {
    warning(lb, ": re-verified avgPower (", round(pw, 4),
            ") is within 2 MC standard errors of target_power (", target_power,
            ") - margin is thin")
  }
}


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------
looks_cols <- c("Fixed" = "#000000", "One interim" = "#0072B2", "Two interims" = "#D55E00")

p_euii <- ggplot(summary_df, aes(x = DeltaTarget, y = EUII,
                                 color = Looks, group = Looks)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.4) +
  scale_color_manual(values = looks_cols) +
  labs(x = expression(delta[target]), y = "EUII at the calibration point",
       title = paste0("EUII of designs calibrated at their own target (avgT1E = ",
                      target_t1e, ", avgPower = ", target_power, ")")) +
  my_theme

p_ratio <- ggplot(summary_df, aes(x = DeltaTarget, y = EUII_ratio_to_fixed,
                                  color = Looks, group = Looks)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.4) +
  scale_color_manual(values = looks_cols) +
  labs(x = expression(delta[target]), y = "EUII relative to the fixed design",
       title = "Does the sequential advantage depend on the effect size designed for?") +
  my_theme

# n_A is fixed while n_C shrinks as the target grows, so the borrowing regime
# itself changes across the grid. The single delta comparison cannot show this.
p_borrow <- ggplot(summary_df, aes(x = DeltaTarget, y = n_A_over_nC,
                                   color = Looks, group = Looks)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.4) +
  scale_color_manual(values = looks_cols) +
  labs(x = expression(delta[target]),
       y = expression(n[A] / n[C]),
       title = "Borrowing regime: historical ESS relative to the concurrent control arm") +
  my_theme

ggsave(file.path(out_dir, "euii_vs_target.pdf"),        p_euii,   width = 8, height = 5)
ggsave(file.path(out_dir, "euii_ratio_vs_target.pdf"),  p_ratio,  width = 8, height = 5)
ggsave(file.path(out_dir, "borrowing_vs_target.pdf"),   p_borrow, width = 8, height = 5)

saveRDS(list(designs_cal   = designs_cal,
             summary_df    = summary_df,
             combos        = combos,
             delta_targets = delta_targets,
             target_t1e    = target_t1e,
             target_power  = target_power,
             n_sim_calib   = n_sim_calib,
             n_sim_check   = n_sim_check,
             prior_H1      = prior_H1,
             n_A           = n_A,
             mc_seed       = mc_seed),
        file = file.path(out_dir, "results.RDS"))
cat("\nResults saved to", file.path(out_dir, "results.RDS"), "\n")
