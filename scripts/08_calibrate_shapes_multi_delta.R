# =============================================================================
# Multi delta calibration - boundary shape
# =============================================================================
# Companion to 07_calibrate_designs_multi_delta.R. There the boundary is held
# constant and the number of looks varies; here the number of looks is held at
# K = 3 (two interims) and the shape of
# the boundary varies instead, so the two comparisons stay unconfounded.
#
# As in script 07, (n, p) is recalibrated at every delta_target and the EUII
# is read at the design's own calibration point, out of the same simulation
# call that produces avgT1E and avgPower.
# =============================================================================

library(dplyr)
library(ggplot2)
source("scripts/00_shared_setup.R")            # sigma, p_MAP, p_vague, DV, my_theme
source("scripts/00_functions.R")               # calibrate_design, compute_euii
source("scripts/00_boundary_builders.R")       # OBF, Haybittle-Peto, Wang Tsiatis

out_dir <- "Output/08_calibrate_shapes_multi_delta"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
target_t1e    <- 0.05
target_power  <- 0.90
delta_targets <- c(50, 60, 70, 80, 90)

n_sim_calib <- 1e6
n_sim_check <- 2e6
mc_seed     <- 123
prior_H1    <- c(0.01, 0.1, 0.5)

# Two interims, matching the original trial protocol. Cumulative totals at the three looks are
# 20, 40 and 60, so the information fractions are 1/3, 2/3 and 1, which is what
# the builders assume. Only a starting point for the search over n.
base_schedule <- list(n1 = c(13, 27, 40), n2 = c(7, 13, 20))

shapes <- list(
  "Pocock"         = list(builder = NULL,                p_interval = c(0.6, 0.999)),
  "OBF"            = list(builder = builder_obf_2int,    p_interval = c(0.6, 0.999)),
  "Haybittle-Peto" = list(builder = builder_hp_2int,     p_interval = c(0.6, 0.999)),
  "Wang Tsiatis"   = list(builder = builder_wt_2int,     p_interval = c(0.6, 0.999))
)
shape_names <- names(shapes)

combos <- expand.grid(Shape       = shape_names,
                      DeltaTarget = delta_targets,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
combos$Label <- paste(combos$Shape, combos$DeltaTarget, sep = " | d")

cores <- min(32L, nrow(combos))
if (.Platform$OS.type != "unix") cores <- 1L

n_A <- as.numeric(ess(p_MAP))

cat("===== Multi delta calibration, axis 2 (boundary shape, K = 3) =====\n")
cat(nrow(combos), "combinations,", length(delta_targets), "targets x",
    length(shapes), "shapes, cores =", cores, "\n")
cat("targets:", paste(delta_targets, collapse = ", "), "\n\n")


# ---------------------------------------------------------------------------
# Calibration, one cell per (shape x delta target)
# ---------------------------------------------------------------------------
designs_cal <- parallel::mclapply(seq_len(nrow(combos)), function(i) {
  sh <- shapes[[combos$Shape[i]]]
  calibrate_design(target_t1e = target_t1e, target_power = target_power,
                   delta_power = combos$DeltaTarget[i],
                   n1_base = base_schedule$n1, n2_base = base_schedule$n2,
                   prior_cntrl = p_MAP, prior_treat = p_vague,
                   criterion = "SC", p_interval = sh$p_interval,
                   decisions_builder = sh$builder,
                   n_sim = n_sim_calib, seed = mc_seed)
}, mc.cores = cores)
names(designs_cal) <- combos$Label

bad <- vapply(designs_cal, inherits, logical(1), "try-error")
if (any(bad)) {
  cat("Calibration failed for:\n")
  cat(paste0("   ", combos$Label[bad], "\n"), sep = "")
  stop("calibrate_design failed for ", sum(bad), " of ", nrow(combos), " combinations.")
}


# ---------------------------------------------------------------------------
# Verification and EUII, both from the same simulation call
# ---------------------------------------------------------------------------
cat("===== Verification at a fresh seed, n_sim =", n_sim_check, "=====\n\n")

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
    seed           = mc_seed + 1
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
# Assemble the summary table
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
  mutate(EUII_ratio_to_pocock = EUII / EUII[Shape == "Pocock"]) |>
  ungroup() |>
  mutate(Shape = factor(Shape, levels = shape_names)) |>
  arrange(DeltaTarget, Shape)


# ---------------------------------------------------------------------------
# Print, and flag any thin power margin
# ---------------------------------------------------------------------------
cat("Calibrated to avgT1E =", target_t1e, "and avgPower =", target_power,
    "at each design's own delta target, K = 3 throughout.\n")
cat("EUII is quoted at that same target, at Pr(H1) = 0.5.\n\n")

print(as.data.frame(summary_df |>
        select(DeltaTarget, Shape, n_max, p, avgT1E, avgPower, EN, EUII,
               EUII_ratio_to_pocock) |>
        mutate(across(c(p, avgT1E, avgPower, EUII, EUII_ratio_to_pocock),
                      \(x) round(x, 4)),
               EN = round(EN, 2))),
      row.names = FALSE)

for (i in seq_len(nrow(summary_df))) {
  pw <- summary_df$avgPower[i]
  se <- sqrt(pw * (1 - pw) / n_sim_check)
  lb <- paste0(summary_df$Shape[i], " at delta = ", summary_df$DeltaTarget[i])
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
shape_cols <- c("Pocock" = "#000000", "OBF" = "#0072B2", "Haybittle-Peto" = "#D55E00",
                "Wang Tsiatis" = "#009E73")

p_euii <- ggplot(summary_df, aes(x = DeltaTarget, y = EUII,
                                 color = Shape, group = Shape)) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.4) +
  scale_color_manual(values = shape_cols) +
  labs(x = expression(delta[target]), y = "EUII at the calibration point",
       title = paste0("EUII by boundary shape, each calibrated at its own target (K = 3)")) +
  my_theme

p_ratio <- ggplot(summary_df, aes(x = DeltaTarget, y = EUII_ratio_to_pocock,
                                  color = Shape, group = Shape)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.4) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.4) +
  scale_color_manual(values = shape_cols) +
  labs(x = expression(delta[target]), y = "EUII relative to Pocock",
       title = "Does the ranking of boundary shapes depend on the effect size designed for?") +
  my_theme

ggsave(file.path(out_dir, "euii_vs_target.pdf"),       p_euii,  width = 8, height = 5)
ggsave(file.path(out_dir, "euii_ratio_vs_target.pdf"), p_ratio, width = 8, height = 5)

saveRDS(list(designs_cal   = designs_cal,
             summary_df    = summary_df,
             combos        = combos,
             delta_targets = delta_targets,
             shape_names   = shape_names,
             base_schedule = base_schedule,
             target_t1e    = target_t1e,
             target_power  = target_power,
             n_sim_calib   = n_sim_calib,
             n_sim_check   = n_sim_check,
             prior_H1      = prior_H1,
             n_A           = n_A,
             mc_seed       = mc_seed),
        file = file.path(out_dir, "results.RDS"))
cat("\nResults saved to", file.path(out_dir, "results.RDS"), "\n")
