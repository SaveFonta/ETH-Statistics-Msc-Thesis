# =============================================================================
# STEP 2 - EUII across delta for the calibrated designs
# =============================================================================
# Loads the designs calibrated in 07_calibrate_designs.R (all sharing the
# same avgT1E and the same avgPower at delta_MCID) and compares them across the
# whole range of treatment effects. Since the error rates coincide by
# construction, the EUII differences are attributable to the stopping structure
# alone, priced per patient.
#
# Run from the Thesis/ root, AFTER 2.5:
#   Rscript scripts/08_calibrated_euii_comparison.R
# =============================================================================

library(RBesT)
library(dplyr)
library(ggplot2)
source("scripts/00_functions_monte_carlo.R")


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
sigma <- 88

n_sim_eval   <- 2e5
delta_values <- c(0, seq(10, 100, by = 10))   # delta = 0 required by compute_euii
prior_H1     <- c(0.01, 0.1, 0.5)
ph1_ref      <- 0.5

p_MAP <- mixnorm(
  c(0.4848, -52.457, 21.154), c(0.4598, -47.465, 7.843), c(0.0554, -50.355, 48.164),
  sigma = sigma, param = "ms"
)
p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")

my_theme <- theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())

cal <- readRDS("Output/CALIBRATED_designs.RDS")
designs_cal <- cal$designs_cal
design_names <- names(designs_cal)
cat("Loaded", length(designs_cal), "calibrated designs (avgT1E =", cal$target_t1e,
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

# EUII for every prior_H1, line drawn at ph1_ref with a band over the grid
df_euii <- lapply(design_names, function(nm) {
  e <- compute_euii(res_cal[[nm]], prior_H1 = prior_H1)
  bind_rows(lapply(names(e), function(g) {
    e[[g]] |> mutate(prior_H1 = as.numeric(g), Design = nm)
  }))
}) |> bind_rows() |>
  mutate(Design = factor(Design, levels = design_names))

band <- df_euii |>
  group_by(Design, Delta) |>
  summarise(lo  = min(EUII),
            hi  = max(EUII),
            mid = EUII[prior_H1 == ph1_ref],
            .groups = "drop")


# ---------------------------------------------------------------------------
# Summary at delta_MCID and the EUII figure
# ---------------------------------------------------------------------------
cat("\nEUII at delta =", cal$delta_MCID, "(prior_H1 = ", ph1_ref, "):\n")
for (nm in design_names) {
  eu <- band |> filter(Design == nm, Delta == cal$delta_MCID) |> pull(mid)
  ov <- res_cal[[nm]][[paste0("delta.", cal$delta_MCID)]]$Overall
  cat(" ", format(nm, width = 12),
      " EUII =", round(eu, 4),
      "  E[N] =", round(ov[["EN_t"]] + ov[["EN_c"]], 1), "\n")
}

design_cols <- c("Fixed" = "#000000", "One interim" = "#0072B2",
                 "Two interims" = "#D55E00")

p_euii <- ggplot(band) +
  geom_ribbon(aes(x = Delta, ymin = lo, ymax = hi, fill = Design), alpha = 0.30) +
  geom_line(aes(x = Delta, y = mid, color = Design), linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = design_cols) +
  scale_fill_manual(values = design_cols, guide = "none") +
  labs(x = expression(delta),
       y = "EUII",
       title = paste0("EUII at equal avgT1E (", cal$target_t1e,
                      ") and equal avgPower (", cal$target_power,
                      " at delta = ", cal$delta_MCID, ")")) +
  my_theme

print(p_euii)

saveRDS(list(res_cal = res_cal, df_euii = df_euii, band = band),
        file = "Output/CALIBRATED_euii.RDS")
cat("\nResults saved to Output/CALIBRATED_euii.RDS\n")
