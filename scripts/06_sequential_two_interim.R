# =============================================================================
# Sequential Design with TWO interim analyses — MC evaluation of AVG OC and EUII
#
# Same analysis as 05_sequential_one_interim.R (shares every summary/plotting
# helper via 00_functions.R), but the trial has three looks, timed
# as in the original protocol (Hueber et al. 2012): total enrollment of 17, 41
# and 60 patients at the two interims and the final analysis. With the 2:1
# allocation this gives cumulative arm sizes that are only approximately 2:1 at
# the interims (11:6 and 27:14).
# Futility can stop the trial at both interims, efficacy at every look.
# Results saved to Output/06_sequential_two_interim/avg.RDS and .../euii.RDS.
#
# Structure: (1) data / simulation, (2) saveRDS, (3) plots and tables, each
# labelled with the manuscript.Rnw chunk it feeds.
# =============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)   # plot_layout(guides = "collect") for a single shared legend
library(parallel)
source("scripts/00_shared_setup.R")          # sigma, priors, decision criteria, theme, crit_cols...
source("scripts/00_functions.R")

out_dir <- "Output/06_sequential_two_interim"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Schedule specific to the two-interim design
# ---------------------------------------------------------------------------
n1_seq <- c(11, 27, 40)   # treatment arm cumulative sizes (three looks)
n2_seq <- c(6, 14, 20)    # control arm cumulative sizes, stage totals 17, 41, 60

delta_values <- c(0, seq(10, 100, by = 10))
mc_seed      <- 123    # shared by every sweep (paired comparison)

# Prior probability of H1 in the EUII. It does not enter the simulation, only the
# processing with compute_euii, so the whole grid comes from a single sweep.
prior_H1 <- c(0.01, 0.1, 0.5)


# ---------------------------------------------------------------------------
# Decision criteria per stage (success/futility criteria themselves come from
# 00_shared_setup.R: sign.crit, dual.crit, fut.crit)
# ---------------------------------------------------------------------------
decisions_fut <- list(
  list(success = dual.crit, futility = fut.crit),   # interim 1
  list(success = dual.crit, futility = fut.crit),   # interim 2
  list(success = dual.crit, futility = NULL)        # final look
)
decisions_no_fut <- list(
  list(success = dual.crit, futility = NULL),
  list(success = dual.crit, futility = NULL),
  list(success = dual.crit, futility = NULL)
)

decisions_fut.single <- list(
  list(success = sign.crit, futility = fut.crit),   # interim 1
  list(success = sign.crit, futility = fut.crit),   # interim 2
  list(success = sign.crit, futility = NULL)        # final look
)
decisions_no_fut.single <- list(
  list(success = sign.crit, futility = NULL),
  list(success = sign.crit, futility = NULL),
  list(success = sign.crit, futility = NULL)
)


# =============================================================================
# DATA — Run one MC sweep per decision criterion
# =============================================================================
# Four designs are compared. All use the SAME seed and the SAME n_sim, so theta_C
# and the stage 1 data are common to all four and the comparison is paired.
n_sim_avg <- 1e6

designs <- list(
  "SC"            = decisions_no_fut.single,
  "SC + Futility" = decisions_fut.single,
  "DC"            = decisions_no_fut,
  "DC + Futility" = decisions_fut
)

sweeps <- lapply(design_labels, function(nm) {
  cat("\n=== ", nm, " ===")
  run_sequential_sweep(designs[[nm]], combos, analysis_priors, design_priors,
                       prior_1 = p_vague, n1_seq = n1_seq, n2_seq = n2_seq,
                       delta_values = delta_values,
                       n_sim_avg = n_sim_avg, mc_seed = mc_seed)
})
names(sweeps) <- design_labels


# ---------------------------------------------------------------------------
# DATA — Stack the four sweeps, tagging each with its Criterion
# ---------------------------------------------------------------------------
df_oc <- bind_rows(lapply(design_labels, function(nm) {
  summarise_sweep(sweeps[[nm]], combos, extract_sweep_total) |> mutate(Criterion = nm)
})) |>
  group_by(Criterion, Analysis_Prior, Design_Prior) |>
  mutate(avgT1E = Power[delta == 0]) |>
  ungroup()

df_stage <- bind_rows(lapply(design_labels, function(nm) {
  summarise_sweep(sweeps[[nm]], combos, extract_sweep_per_stage) |> mutate(Criterion = nm)
})) |>
  mutate(delta = as.numeric(delta))

df_oc <- df_oc |> mutate(
  Criterion      = factor(Criterion, levels = design_labels),
  Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
  Design_Prior   = factor(Design_Prior,   levels = dprior_names)
)

df_stage <- df_stage |> mutate(
  Criterion      = factor(Criterion, levels = design_labels),
  Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
  Design_Prior   = factor(Design_Prior,   levels = dprior_names)
)


# =============================================================================
# DATA — Sequential EUII, correctly specified prior
# =============================================================================
# hERE, WE ONLY USE THE correctly specified analysis prior (so MAP for both analysis and design prior)

df_euii <- bind_rows(lapply(design_labels, function(nm) {
  sweep_euii_table(sweeps[[nm]], combos, prior_H1, omega_map, prior_names, dprior_names) |>
    mutate(Criterion = nm)
})) |>
  filter(Analysis_Prior == "MAP", Design_Prior == "MAP") |>
  mutate(Criterion = factor(Criterion, levels = design_labels),
         Base      = if_else(grepl("^SC", as.character(Criterion)), "SC", "DC"))

## Fixed baselines
fixed_base <- read_fixed_baselines()

df_euii <- df_euii |>
  left_join(fixed_base, by = c("Delta", "Base")) |>
  mutate(EUII_ratio = EUII / EUII_fixed)


# =============================================================================
# SAVE
# =============================================================================
saveRDS(list(df_oc = df_oc, df_stage = df_stage, sweeps = sweeps),
        file = file.path(out_dir, "avg.RDS"))
cat("\nMC average and per stage results saved to", file.path(out_dir, "avg.RDS"), "\n")

saveRDS(list(df_euii = df_euii, prior_H1 = prior_H1),
        file = file.path(out_dir, "euii.RDS"))
cat("MC EUII results saved to", file.path(out_dir, "euii.RDS"), "\n")


# =============================================================================
# PLOTS AND TABLES
# =============================================================================

# --- manuscript.Rnw chunk: plot_seq2_oc -----------------------------------
# Average T1E and average Power, one panel per criterion
t1e <- df_oc |> filter(delta == 0)
p_t1e <- oc_facet_layout(
  ggplot(t1e,
         aes(x = Design_Prior, y = Power,
             color = Analysis_Prior, group = Analysis_Prior,
             linetype = Analysis_Prior == "Vague"))) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.5) +
  labs(x = "Design Prior", y = "Average T1E")

print(p_t1e)

p_pow <- oc_facet_layout(
  ggplot(filter(df_oc, delta == delta_MCID),
         aes(x = Design_Prior, y = Power,
             color = Analysis_Prior, group = Analysis_Prior,
             linetype = Analysis_Prior == "Vague"))) +
  labs(x = "Design Prior",
       y = bquote("Average Power"))

print(p_pow)

# Expected sample size. Futility leaves the average T1E and the average Power
# basically unchanged, because it only stops trials that were going to fail
# anyway. --> from here we don't see the advantage of futility analyis. Luckily WE
# HAVE DEVELOPED THE EUII


# --- manuscript.Rnw chunk: plot_seq2_stage --------------------------------
# What the interims actually do: how often the trial stops early for efficacy when
# the treatment works, and how often it stops early for futility under the null.
# With two interims there is one row per interim and per metric, four rows in total.
df_stage1 <- bind_rows(
  df_stage |> filter(Stage <= 2, delta == delta_MCID) |>
    transmute(Criterion, Analysis_Prior, Design_Prior, Prob = P_Succ,
              Metric = paste0("Efficacy stop, interim ", Stage,
                              " (delta = ", delta_MCID, ")")),
  df_stage |> filter(Stage <= 2, delta == 0) |>
    transmute(Criterion, Analysis_Prior, Design_Prior, Prob = P_Fut,
              Metric = paste0("Futility stop, interim ", Stage, " (delta = 0)"))
)

p_stage <- ggplot(df_stage1,
                  aes(x = Design_Prior, y = Prob,
                      color = Analysis_Prior, group = Analysis_Prior,
                      linetype = Analysis_Prior == "Vague")) +
  geom_point(size = 2.2) +
  geom_line(linewidth = 0.9) +
  facet_grid(Metric ~ Criterion) +
  scale_linetype_vague() +
  scale_color_prior() +
  labs(x = "Design Prior", y = "Probability of stopping at the interim") +
  my_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p_stage)

# NOTE: tab_seq2_summary in the manuscript builds its table directly from
# Output/06_sequential_two_interim/avg.RDS, so there is no corresponding
# table-building block here.


# --- Own diagnostic only, no manuscript chunk reads this -------------------
# Absolute EUII of the two-interim designs against the fixed-design references.
# The manuscript's two-interim EUII figure is plot_seq2_vs_one below (ratio to
# the one-interim designs) instead of this absolute view.
crit_levels <- c(design_labels, "Fixed SC", "Fixed DC")

# The two fixed designs carry no prior_H1: in a fixed design N is constant. They
# are single reference lines, with no band.
df_fixed_lines <- fixed_base |>
  transmute(Delta, EUII = EUII_fixed,
            Criterion = factor(paste("Fixed", Base), levels = crit_levels))

df_euii_leveled <- df_euii |> mutate(Criterion = factor(as.character(Criterion), levels = crit_levels))
df_band_abs   <- value_band(df_euii_leveled, "EUII",       c("Criterion", "Delta"))
df_band_ratio <- value_band(df_euii_leveled, "EUII_ratio", c("Criterion", "Delta"))

# top panel: absolute EUII, with the two fixed designs as references
p_euii_abs <- ggplot() +
  geom_ribbon(data = df_band_abs,
              aes(x = Delta, ymin = lo, ymax = hi, fill = Criterion), alpha = 0.30) +
  geom_line(data = df_band_abs,
            aes(x = Delta, y = mid, color = Criterion, linetype = Criterion),
            linewidth = 1.2) +
  geom_line(data = df_fixed_lines,
            aes(x = Delta, y = EUII, color = Criterion, linetype = Criterion),
            linewidth = 1.0) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", linewidth = 0.4) +
  scale_color_manual(values = crit_cols, drop = FALSE) +
  scale_linetype_manual(values = crit_ltys, drop = FALSE) +
  scale_fill_manual(values = crit_cols, guide = "none") +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  coord_cartesian(xlim = c(DV, max(delta_values))) +
  labs(x = NULL, y = "EUII") +
  my_theme +
  theme(axis.text.x = element_blank())

# bottom panel: ratio to the corresponding fixed design
p_euii_ratio <- ggplot(df_band_ratio) +
  geom_ribbon(aes(x = Delta, ymin = lo, ymax = hi, fill = Criterion), alpha = 0.30) +
  geom_line(aes(x = Delta, y = mid, color = Criterion, linetype = Criterion),
            linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = crit_cols) +
  scale_linetype_manual(values = crit_ltys) +
  scale_fill_manual(values = crit_cols, guide = "none") +
  guides(color = "none", linetype = "none") +
  coord_cartesian(xlim = c(DV, max(delta_values))) +
  labs(x = expression(delta), y = "EUII ratio") +
  my_theme

# patchwork stacks the panels and collects the guides, so they share one legend
print((p_euii_abs / p_euii_ratio) +
        patchwork::plot_layout(guides = "collect") &
        theme(legend.position = "bottom"))


# --- manuscript.Rnw chunk: plot_seq2_vs_one (bottom panel: eff. sample size) ---
# Effective sample size, 1 / E[1/N | outcome]. As in the EUII panels, the line
# is drawn at prior_H1 = 0.5 and there is a ribbon for the others.
n_min <- n1_seq[1] + n2_seq[1]                              #  first interim total, 20
n_mid <- n1_seq[2] + n2_seq[2]                              #  second interim total, 40
n_max <- n1_seq[length(n1_seq)] + n2_seq[length(n2_seq)]    #  final total, 60

sample_size <- df_euii |>
  filter(grepl("Futility", Criterion)) |>
  select(Delta, prior_H1, Criterion, E_invN_sig, E_invN_nonsig) |>
  tidyr::pivot_longer(cols = c(E_invN_sig, E_invN_nonsig),
                      names_to = "E_inv_type", values_to = "E_inv_val") |>
  mutate(
    N_eff     = 1 / E_inv_val,
    Outcome   = factor(if_else(E_inv_type == "E_invN_sig", "Success", "Non success"),
                       levels = c("Success", "Non success")),
    Criterion = factor(as.character(Criterion), levels = design_labels)
  )
band_N <- value_band(sample_size, "N_eff", c("Criterion", "Outcome", "Delta"))

p_N_eff <- ggplot(band_N) +
  geom_ribbon(aes(x = Delta, ymin = lo, ymax = hi, fill = Criterion), alpha = 0.30) +
  geom_line(aes(x = Delta, y = mid, color = Criterion), linewidth = 1.2) +
  geom_hline(yintercept = n_max, linetype = "dashed", color = "black",  linewidth = 0.4) +
  geom_hline(yintercept = n_mid, linetype = "dotted", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = n_min, linetype = "dotted", color = "grey40", linewidth = 0.4) +
  annotate("text", x = DV + 2, y = n_max - 1.3, hjust = 0, size = 3, color = "grey30",
           label = paste0("without futility (= ", n_max, ")")) +
  annotate("text", x = DV + 2, y = n_mid + 1.1, hjust = 0, size = 3, color = "grey50",
           label = paste0("second interim (= ", n_mid, ")")) +
  annotate("text", x = DV + 2, y = n_min + 1.1, hjust = 0, size = 3, color = "grey50",
           label = paste0("first interim (= ", n_min, ")")) +
  facet_wrap(~Outcome) +
  scale_color_manual(values = crit_cols) +
  scale_fill_manual(values = crit_cols, guide = "none") +
  coord_cartesian(xlim = c(DV, max(delta_values)), ylim = c(n_min, n_max)) +
  labs(x = expression(delta),
       y = expression("Effective sample size  " * (E * "[1/N | outcome]")^-1)) +
  my_theme

print(p_N_eff)


# --- manuscript.Rnw chunk: plot_seq2_vs_one (top panel: ratio to one interim) ---
# The natural question of this script: does the second interim add information per
# patient beyond the first? The one interim results of 05_sequential_one_interim.R are read
# from Output/05_sequential_one_interim/euii.RDS and the ratio EUII(two interims) / EUII(one interim)
# is computed per design and prior_H1. Both runs share the same seed, so theta_C is
# common and the comparison is paired.

df_euii_1int <- readRDS("Output/05_sequential_one_interim/euii.RDS")$df_euii |>
  select(Delta, prior_H1, Criterion, EUII_1int = EUII)

df_ratio_12 <- df_euii |>
  select(Delta, prior_H1, Criterion, EUII) |>
  mutate(Criterion = as.character(Criterion)) |>
  left_join(df_euii_1int |> mutate(Criterion = as.character(Criterion)),
            by = c("Delta", "prior_H1", "Criterion")) |>
  mutate(ratio_12 = EUII / EUII_1int)

band_12 <- df_ratio_12 |>
  mutate(Criterion = factor(Criterion, levels = design_labels)) |>
  value_band("ratio_12", c("Criterion", "Delta"))

p_ratio_12 <- ggplot(band_12) +
  geom_ribbon(aes(x = Delta, ymin = lo, ymax = hi, fill = Criterion), alpha = 0.30) +
  geom_line(aes(x = Delta, y = mid, color = Criterion), linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = crit_cols) +
  scale_fill_manual(values = crit_cols, guide = "none") +
  coord_cartesian(xlim = c(DV, max(delta_values))) +
  labs(x = expression(delta),
       y = "EUII (two interims) / EUII (one interim)") +
  my_theme

print(p_ratio_12)
