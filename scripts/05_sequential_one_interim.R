# =============================================================================
# Sequential Design One interim — Monte Carlo evaluation of AVERAGE OC and EUII
# =============================================================================
# Structure: (1) data / simulation, (2) saveRDS, (3) plots and tables, each
# labelled with the manuscript.Rnw chunk it feeds.
# =============================================================================

library(dplyr)
library(ggplot2)
library(patchwork)   # plot_layout(guides = "collect") for a single shared legend
library(parallel)
source("scripts/00_shared_setup.R")
source("scripts/00_functions.R")

out_dir <- "Output/05_sequential_one_interim"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Schedule
# ---------------------------------------------------------------------------
n1_seq <- c(20, 40)   # trt arm cumulative sizes
n2_seq <- c(10, 20)   # cntrl arm cumulative sizes

delta_values <- c(0, seq(10, 100, by = 10))
mc_seed      <- 123    # shared by every sweep (paired comparison)

# Prior probability of H1 in the EUII. It does not enter the simulation, only the
# processing with compute_euii.
prior_H1 <- c(0.01, 0.1, 0.5)


# ---------------------------------------------------------------------------
# Decision criteria per stage
# ---------------------------------------------------------------------------
decisions_fut <- list(
  list(success = dual.crit, futility = fut.crit),
  list(success = dual.crit, futility = NULL)
)
decisions_no_fut <- list(
  list(success = dual.crit, futility = NULL),
  list(success = dual.crit, futility = NULL)
)

decisions_fut.single <- list(
  list(success = sign.crit, futility = fut.crit),
  list(success = sign.crit, futility = NULL)
)
decisions_no_fut.single <- list(
  list(success = sign.crit, futility = NULL),
  list(success = sign.crit, futility = NULL)
)


# =============================================================================
# DATA — Run one MC sweep per decision criterion
# =============================================================================
n_sim_avg <- 1e6

designs <- list(
  "SC"            = decisions_no_fut.single,
  "SC + Futility" = decisions_fut.single,
  "DC"            = decisions_no_fut,
  "DC + Futility" = decisions_fut
)

# sweep over all decision lists (i.e. designs)
sweeps <- lapply(design_labels, function(nm) {
  cat("\n=== ", nm, " ===")
  run_sequential_sweep(decisions_list = designs[[nm]], combos, analysis_priors, design_priors,
                       prior_1 = p_vague, n1_seq = n1_seq, n2_seq = n2_seq,
                       delta_values = delta_values,
                       n_sim_avg = n_sim_avg, mc_seed = mc_seed)
})
names(sweeps) <- design_labels
names(sweeps$SC[[1]])
sweeps$SC[["MAP | Dirac (-50)"]]  |> names()

# ---------------------------------------------------------------------------
# DATA — Stack the four sweeps, adding each with its Criterion
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
# DATA — EUII for correctly specified prior
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

# --- manuscript.Rnw chunks: seq_oc_setup, plot_seq_oc --------------------
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


# --- manuscript.Rnw chunk: plot_seq_stage ---------------------------------
# What the interim actually does: how often it stops early for efficacy when the
# treatment works, and how often it stops early for futility under the null.
df_stage1 <- bind_rows(
  df_stage |> filter(Stage == 1, delta == delta_MCID) |>
    transmute(Criterion, Analysis_Prior, Design_Prior, Prob = P_Succ,
              Metric = paste0("Early efficacy stop (delta = ", delta_MCID, ")")),
  df_stage |> filter(Stage == 1, delta == 0) |>
    transmute(Criterion, Analysis_Prior, Design_Prior, Prob = P_Fut,
              Metric = "Early futility stop (delta = 0)")
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


# --- manuscript.Rnw chunk: tab_seq_summary --------------------------------
tab_stage <- df_stage |>
  filter(Analysis_Prior == "MAP", Design_Prior == "MAP", delta %in% c(0, 60)) |>
  group_by(Criterion, delta) |>
  summarise(
    P_succ_interim = P_Succ[Stage == 1],
    P_fut_interim  = P_Fut[Stage == 1],
    P_succ_final   = P_Succ[Stage == max(Stage)],
        .groups = "drop"
  ) |>
  mutate(Total_succ = P_succ_interim + P_succ_final)

tab_designs <- tab_stage |>
  left_join(df_oc |>
              filter(Analysis_Prior == "MAP", Design_Prior == "MAP") |>
              select(Criterion, delta, EN),
            by = c("Criterion", "delta")) |>
  arrange(delta, Criterion)

knitr::kable(tab_designs, format = "pipe", digits = 3)


# --- manuscript.Rnw chunk: plot_seq_euii (top/middle panels) -------------
# Left: EUII of the four sequential designs, with the two fixed designs dashed.
# Right: ratio of each sequential design to ITS OWN fixed counterpart, so SC and
# SC + Futility are divided by the fixed SC, and DC and DC + Futility by the fixed DC.
# colour encodes the DESIGN (crit_cols/crit_ltys from 00_shared_setup.R)
crit_levels <- c(design_labels, "Fixed SC", "Fixed DC")

# The two fixed designs carry no prior_H1: in a fixed design N is constant. 
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


# --- manuscript.Rnw chunk: plot_seq_euii (bottom panel: effective sample size) ---
# 1 / E[1/N | outcome]. As in the EUII panels, the line is drawn at prior_H1 = 0.5
# and there is a ribbon for the others.
n_min <- n1_seq[1] + n2_seq[1]                              #  30
n_max <- n1_seq[length(n1_seq)] + n2_seq[length(n2_seq)]    #  60

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
  geom_hline(yintercept = n_min, linetype = "dotted", color = "grey40", linewidth = 0.4) +
  annotate("text", x = DV + 2, y = n_max - 1.3, hjust = 0, size = 3, color = "grey30",
           label = paste0("without futility (= ", n_max, ")")) +
  annotate("text", x = DV + 2, y = n_min + 1.1, hjust = 0, size = 3, color = "grey50",
           label = paste0("interim (= ", n_min, ")")) +
  facet_wrap(~Outcome) +
  scale_color_manual(values = crit_cols) +
  scale_fill_manual(values = crit_cols, guide = "none") +
  coord_cartesian(xlim = c(DV, max(delta_values)), ylim = c(n_min, n_max)) +
  labs(x = expression(delta),
       y = expression("Effective sample size  " * (E * "[1/N | outcome]")^-1)) +
  my_theme

print(p_N_eff)
