# =============================================================================
# STEP 2 - EUII across delta for the calibrated boundary shapes
# =============================================================================
# Loads the shapes calibrated in 09_calibrate_boundary_shapes.R (all sharing
# the same avgT1E, and the same avgPower at delta_MCID, at a fixed K = 2
# schedule) and compares them across the whole range of treatment effects.
# Since the error rates coincide by construction, the EUII differences are
# attributable to the boundary shape alone, priced per patient.
#
# Run from the Thesis/ root, AFTER 09_calibrate_boundary_shapes.R:
#   Rscript scripts/10_boundary_shapes_euii_comparison.R
# =============================================================================

library(dplyr)
library(ggplot2)
source("scripts/00_shared_setup.R")          # sigma, p_MAP, p_vague, my_theme, ph1_ref
source("scripts/00_functions.R")

out_dir <- "Output/10_boundary_shapes_euii_comparison"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
n_sim_eval   <- 2e5
delta_values <- c(0, seq(10, 100, by = 10))   # delta = 0 required by compute_euii
prior_H1     <- c(0.01, 0.1, 0.5)

cal <- readRDS("Output/09_calibrate_boundary_shapes/designs.RDS")
designs_cal <- cal$designs_cal
shape_names <- cal$shape_names
cat("Loaded", length(designs_cal), "calibrated shapes (avgT1E =", cal$target_t1e,
    ", avgPower =", cal$target_power, "at delta =", cal$delta_MCID, ")\n")


# ---------------------------------------------------------------------------
# Evaluate every calibrated shape on the full delta grid
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
df_euii <- lapply(shape_names, function(nm) {
  e <- compute_euii(res_cal[[nm]], prior_H1 = prior_H1)
  bind_rows(lapply(names(e), function(g) {
    e[[g]] |> mutate(prior_H1 = as.numeric(g), Shape = nm)
  }))
}) |> bind_rows() |>
  mutate(Shape = factor(Shape, levels = shape_names))

band <- value_band(df_euii, "EUII", c("Shape", "Delta"))


# ---------------------------------------------------------------------------
# Summary at delta_MCID and the EUII figure
# ---------------------------------------------------------------------------
cat("\nEUII at delta =", cal$delta_MCID, "(prior_H1 = ", ph1_ref, "):\n")
for (nm in shape_names) {
  eu <- band |> filter(Shape == nm, Delta == cal$delta_MCID) |> pull(mid)
  ov <- res_cal[[nm]][[paste0("delta.", cal$delta_MCID)]]$Overall
  cat(" ", format(nm, width = 16),
      " EUII =", round(eu, 4),
      "  E[N] =", round(ov[["EN_t"]] + ov[["EN_c"]], 1), "\n")
}

shape_cols <- c("Flat" = "#000000", "OBF" = "#0072B2", "Haybittle-Peto" = "#D55E00",
                "Power-family" = "#009E73", "Shi & Yin" = "#CC79A7")

p_euii <- ggplot(band) +
  geom_ribbon(aes(x = Delta, ymin = lo, ymax = hi, fill = Shape), alpha = 0.30) +
  geom_line(aes(x = Delta, y = mid, color = Shape), linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = shape_cols) +
  scale_fill_manual(values = shape_cols, guide = "none") +
  labs(x = expression(delta),
       y = "EUII",
       title = paste0("EUII at equal avgT1E (", cal$target_t1e,
                      ") and equal avgPower (", cal$target_power,
                      " at delta = ", cal$delta_MCID, "), K = 2 fixed")) +
  my_theme

print(p_euii)

saveRDS(list(res_cal = res_cal, df_euii = df_euii, band = band),
        file = file.path(out_dir, "euii.RDS"))
cat("\nResults saved to", file.path(out_dir, "euii.RDS"), "\n")
