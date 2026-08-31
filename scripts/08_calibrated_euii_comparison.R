# =============================================================================
# STEP 2 - EUII across delta for the calibrated designs
# =============================================================================
# Loads the designs calibrated in 07_calibrate_designs.R (all sharing the
# same avgT1E, and the same avgPower at delta_MCID) and compares them across
# the whole range of treatment effects. Since the error rates coincide by
# construction, the EUII differences are attributable to the stopping
# structure (number of looks) alone, priced per patient.
#
# SC only. 07_calibrate_designs.R also calibrates a DC family, but only to
# reach avgPower at delta_MCID: the DC's avgT1E cannot be brought to the same
# target (see the comment there), so its error rates never coincide with the
# SC's, and its EUII across delta is not meaningful (used in the thesis only
# through the DC's calibrated n_max, read directly off designs.RDS - not
# through anything computed here). Simulating its full delta grid here would
# cost about half this script's runtime for output nothing downstream reads.
#
# Run from the Thesis/ root, AFTER 07_calibrate_designs.R:
#   Rscript scripts/08_calibrated_euii_comparison.R
# =============================================================================

library(dplyr)
library(ggplot2)
source("scripts/00_shared_setup.R")          # sigma, p_MAP, p_vague, my_theme, ph1_ref
source("scripts/00_functions.R")

out_dir <- "Output/08_calibrated_euii_comparison"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
n_sim_eval   <- 1e6
delta_values <- c(0, seq(10, 100, by = 10))   # delta = 0 required by compute_euii
prior_H1     <- c(0.01, 0.1, 0.5)

cal <- readRDS("Output/07_calibrate_designs/designs.RDS")
design_names <- grep(" \\| SC$", names(cal$designs_cal), value = TRUE)
designs_cal  <- cal$designs_cal[design_names]
cat("Loaded", length(designs_cal), "calibrated SC designs (avgT1E =", cal$target_t1e,
    ", avgPower =", cal$target_power, "at delta =", cal$delta_MCID, ")\n")


# ---------------------------------------------------------------------------
# Evaluate every calibrated design on the full delta grid
# ---------------------------------------------------------------------------
res_cal <- lapply(designs_cal, function(d) {
  avgoc2_seq_mc.normMix(
    prior_1        = p_vague,
    prior_2        = p_MAP,
    n1_seq         = d$n1_seq,
    n2_seq         = d$n2_seq,
    decisions_list = d$decisions,
    delta          = delta_values,
    design_prior_c = p_MAP,
    n_sim          = n_sim_eval,
    seed           = cal$mc_seed
  )
})

# EUII for every prior_H1, line drawn at ph1_ref with a band over the grid.
# Design = "Looks | SC", Looks split back out for colour.
df_euii <- lapply(design_names, function(nm) {
  e <- compute_euii(res_cal[[nm]], prior_H1 = prior_H1)
  bind_rows(lapply(names(e), function(g) {
    e[[g]] |> mutate(prior_H1 = as.numeric(g), Design = nm)
  }))
}) |> bind_rows() |>
  mutate(Design = factor(Design, levels = design_names),
         Looks  = factor(sub(" \\|.*", "", Design), levels = unique(cal$combos$Looks)),
         Criterion = "SC")

band <- value_band(df_euii, "EUII", c("Design", "Looks", "Criterion", "Delta"))


# ---------------------------------------------------------------------------
# Summary at delta_MCID and the EUII figure
# ---------------------------------------------------------------------------
cat("\nEUII at delta =", cal$delta_MCID, "(prior_H1 = ", ph1_ref, "):\n")
for (nm in design_names) {
  eu <- band |> filter(Design == nm, Delta == cal$delta_MCID) |> pull(mid)
  ov <- res_cal[[nm]][[paste0("delta.", cal$delta_MCID)]]$Overall
  cat(" ", format(nm, width = 24),
      " EUII =", round(eu, 4),
      "  E[N] =", round(ov[["EN_t"]] + ov[["EN_c"]], 1), "\n")
}

looks_cols <- c("Fixed" = "#000000", "One interim" = "#0072B2", "Two interims" = "#D55E00")

p_euii <- ggplot(band) +
  geom_ribbon(aes(x = Delta, ymin = lo, ymax = hi, fill = Looks), alpha = 0.30) +
  geom_line(aes(x = Delta, y = mid, color = Looks), linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", linewidth = 0.5) +
  scale_color_manual(values = looks_cols) +
  scale_fill_manual(values = looks_cols, guide = "none") +
  labs(x = expression(delta),
       y = "EUII",
       title = paste0("EUII at equal avgT1E (", cal$target_t1e,
                      ") and equal avgPower (", cal$target_power,
                      " at delta = ", cal$delta_MCID, "), SC by number of looks")) +
  my_theme

print(p_euii)

saveRDS(list(res_cal = res_cal, df_euii = df_euii, band = band),
        file = file.path(out_dir, "euii.RDS"))
cat("\nResults saved to", file.path(out_dir, "euii.RDS"), "\n")
