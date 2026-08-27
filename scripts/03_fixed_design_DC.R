# =============================================================================
# Fixed Design with Dual Criterion — Conditional and Average OC + EUII
# =============================================================================
#
# it is the same as 02_fixed_design_SC.R exactly (same priors, sample sizes, plots, structure)
# but uses the dual efficacy criterion from Gsponer:
#   Pr(delta > 0  | data) > 0.95  AND  Pr(delta > 50 | data) > 0.50
# instead of the single criterion Pr(delta > 0 | data) > 0.95.
#
# decision2S accepts vector pc/qc for multi-criterion AND logic, and this
# works directly with oc2S and avgoc2S.normMix — no custom functions needed.
#
# Structure: (1) data / simulation, (2) saveRDS, (3) plots and tables, each
# labelled with the manuscript.Rnw chunk it feeds.
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
source("scripts/00_functions.R")
source("scripts/00_shared_setup.R")   # sigma, p_MAP, p_vague, analysis_priors, design_priors,
                                       # prior_names, dprior_names, dual.crit, my_theme, mu_A, delta_MCID


# ---------------------------------------------------------------------------
# Design parameters specific to the fixed design
# ---------------------------------------------------------------------------
n_T <- 40
n_C <- 20

out_dir <- "Output/03_fixed_design_DC"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# =============================================================================
# DATA — Part 1: Conditional T1E and Power
# =============================================================================
theta_c <- seq(-110, 10, by = 2)

cat("Building conditional OC functions...\n")

oc_funs <- lapply(analysis_priors, function(prior_C) {
  oc2S(prior1 = p_vague, prior2 = prior_C,
       n1 = n_T, n2 = n_C, decision = dual.crit)
})

df_cond <- bind_rows(lapply(prior_names, function(pname) {
  f <- oc_funs[[pname]]
  data.frame(
    Prior   = pname,
    theta_C = theta_c,
    T1E     = f(theta_c,              theta_c),
    Power   = f(theta_c - delta_MCID, theta_c)
  )
})) |> mutate(Prior = factor(Prior, levels = prior_names))


# =============================================================================
# DATA — Part 2: Average T1E and Power
# =============================================================================
delta_values <- seq(0, 100, by = 10)

cat("Computing average OC...\n")

rows_avg <- lapply(prior_names, function(aname) {
  avg_fun <- avgoc2S.normMix(
    prior1        = p_vague,
    prior2        = analysis_priors[[aname]],
    n1            = n_T,
    n2            = n_C,
    decision      = dual.crit,
    delta         = 0,
    design_prior2 = design_priors[[1]]
  )
  lapply(dprior_names, function(dname) {
    dp      <- design_priors[[dname]]
    t1e_val <- avg_fun(delta_new = 0, design_prior2_new = dp)
    bind_rows(lapply(delta_values, function(dv) {
      data.frame(
        Analysis_Prior = aname,
        Design_Prior   = dname,
        delta          = dv,
        avgT1E         = t1e_val,
        avgPower       = avg_fun(delta_new = -dv, design_prior2_new = dp)
      )
    }))
  })
})

df_avg <- bind_rows(unlist(rows_avg, recursive = FALSE)) |>
  mutate(
    Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
    Design_Prior   = factor(Design_Prior,   levels = dprior_names)
  )


# =============================================================================
# DATA — Part 3: EUII
# =============================================================================
n_trial <- n_T + n_C   # omega_map comes from 00_shared_setup.R

df_euii <- df_avg |>
  mutate(
    omega    = omega_map[as.character(Analysis_Prior)],
    # log space: qlogis(p) = logit(p) = log(p / (1 - p)), so
    # log_DOR = logit(avgPower) - logit(avgT1E) avoids DOR overflowing to Inf
    log_DOR  = qlogis(avgPower) - qlogis(avgT1E),
    log_EUII = log_DOR / n_trial,
    EUII     = exp(log_EUII)
  )


# =============================================================================
# DATA — Part 4: EUII Ratio: Dual / Single-criterion
# =============================================================================
fixed_data    <- readRDS("Output/02_fixed_design_SC/euii.RDS")
df_euii_fixed <- fixed_data$df_euii |>
  mutate(
    Analysis_Prior = as.character(Analysis_Prior),
    Design_Prior   = as.character(Design_Prior)
  ) |>
  select(Analysis_Prior, Design_Prior, delta, EUII_fixed = EUII)

df_ratio <- df_euii |>
  mutate(
    Analysis_Prior = as.character(Analysis_Prior),
    Design_Prior   = as.character(Design_Prior)
  ) |>
  select(Analysis_Prior, Design_Prior, delta, omega, EUII_dual = EUII) |>
  left_join(df_euii_fixed, by = c("Analysis_Prior", "Design_Prior", "delta")) |>
  mutate(EUII_ratio = EUII_dual / EUII_fixed)

# Safety check: a mismatched join key would silently propagate NA into the ratio.
if (any(is.na(df_ratio$EUII_ratio))) {
  stop("EUII ratio join produced NA - check that Analysis/Design prior labels match between DUAL and FIXED outputs.")
}

# Restore factor ordering for plotting
df_ratio <- df_ratio |>
  mutate(
    Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
    Design_Prior   = factor(Design_Prior,   levels = dprior_names)
  )


# =============================================================================
# SAVE
# =============================================================================
saveRDS(list(df_cond = df_cond), file = file.path(out_dir, "cond.RDS"))
cat("Conditional results saved to", file.path(out_dir, "cond.RDS"), "\n")

saveRDS(list(df_avg = df_avg), file = file.path(out_dir, "avg.RDS"))
cat("Average results saved to", file.path(out_dir, "avg.RDS"), "\n")

saveRDS(list(df_euii = df_euii), file = file.path(out_dir, "euii.RDS"))
cat("EUII results saved to", file.path(out_dir, "euii.RDS"), "\n")

saveRDS(list(df_ratio = df_ratio), file = file.path(out_dir, "vs_SC_euii.RDS"))
cat("EUII ratio saved to", file.path(out_dir, "vs_SC_euii.RDS"), "\n")


# =============================================================================
# PLOTS AND TABLES
# =============================================================================

# --- manuscript.Rnw chunk: plot_cond_both -------------------------------
p_cond_T1E <- ggplot(df_cond,
                     aes(x = theta_C, y = T1E,
                         color = Prior, linetype = Prior == "vague")) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  geom_vline(xintercept = mu_A, linetype = "dotted",
             color = "gray", linewidth = 0.8) +
  scale_linetype_vague() +
  scale_color_prior() +
  labs(x = expression("True Control Mean (" * theta[C] * ")"),
       y = "conditional type I error") +
  my_theme

p_cond_Power <- ggplot(df_cond,
                       aes(x = theta_C, y = Power,
                           color = Prior, linetype = Prior == "vague")) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = mu_A, linetype = "dotted",
             color = "gray", linewidth = 0.8) +
  scale_linetype_vague() +
  scale_color_prior() +
  labs(x = expression("True Control Mean (" * theta[C] * ")"),
       y = bquote("conditional power (" * delta * " = 60)")) +
  my_theme

grid.arrange(p_cond_T1E, p_cond_Power, ncol = 2)


# --- manuscript.Rnw chunks: plot_avg_both, plot_power_delta_both --------
p_avg_T1E <- ggplot(df_avg,
                    aes(x = Design_Prior, y = avgT1E,
                        color = Analysis_Prior, group = Analysis_Prior,
                        linetype = Analysis_Prior == "vague")) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  scale_linetype_vague() +
  scale_color_prior() +
  labs(x = "Design Prior", y = "average type I error") +
  my_theme

p_avg_Power <- ggplot(filter(df_avg, delta == delta_MCID),
                      aes(x = Design_Prior, y = avgPower,
                          color = Analysis_Prior, group = Analysis_Prior,
                          linetype = Analysis_Prior == "vague")) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1.0) +
  scale_linetype_vague() +
  scale_color_prior() +
  labs(x = "Design Prior",
       y = bquote("average power (" * delta * " = 60)")) +
  my_theme

grid.arrange(p_avg_T1E, p_avg_Power, ncol = 2)


# --- manuscript.Rnw chunk: plot_euii_both (top/middle rows) --------------
df_euii_60 <- filter(df_euii, delta == 60)

euii_vague_60 <- df_euii_60 |>
  filter(Analysis_Prior == "vague") |>
  pull(EUII) |>
  unique()

p_EUII_omega <- ggplot(df_euii_60,
                       aes(x = omega, y = EUII,
                           color = Design_Prior, group = Design_Prior)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = euii_vague_60, linetype = "dotted",
             color = "black", linewidth = 0.6) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_color_prior() +
  labs(x = expression("Robustness weight (" * omega * ")"),
       y = expression("EUII (" * delta * " = 60)")) +
  my_theme

print(p_EUII_omega)

p_EUII_delta <- ggplot(df_euii,
                       aes(x = delta, y = EUII,
                           color = Analysis_Prior,
                           linetype = Analysis_Prior == "vague")) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~Design_Prior, nrow = 1) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  scale_linetype_vague() +
  scale_color_prior() +
  labs(x = expression(delta), y = "EUII") +
  my_theme

print(p_EUII_delta)


# --- manuscript.Rnw chunk: plot_euii_both (bottom row: ratio DC/SC) -----
# --- Plot A: ratio vs omega at delta = 60 (one line per design prior) ---
df_ratio_60 <- filter(df_ratio, delta == delta_MCID)

p_ratio_omega <- ggplot(df_ratio_60,
                        aes(x = omega, y = EUII_ratio,
                            color = Design_Prior, group = Design_Prior)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_color_prior() +
  labs(x = expression("Robustness weight (" * omega * ")"),
       y = bquote("EUII"["Dual"] / "EUII"["Sign"] ~ "(" * delta * " = 60)")) +
  my_theme

print(p_ratio_omega)

# --- Plot B: ratio vs delta (facet by design prior, color by analysis prior) ---
p_ratio_delta <- ggplot(df_ratio,
                        aes(x = delta, y = EUII_ratio,
                            color = Analysis_Prior,
                            linetype = Analysis_Prior == "vague")) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~Design_Prior, nrow = 1) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  scale_linetype_vague() +
  scale_color_prior() +
  labs(x = expression(delta),
       y = bquote("EUII"["Dual"] / "EUII"["Sign"])) +
  my_theme

print(p_ratio_delta)
