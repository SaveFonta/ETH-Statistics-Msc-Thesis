# =============================================================================
# Sequential Design (one interim) — Conditional and Average OC + EUII
# =============================================================================
# Run from the Thesis/ root directory:
#   setwd("/home/menelao/Thesis")
#   source("scripts/prior_on_control/2.2_SEQUENTIAL.R")
#
# Part 1 – Conditional OC  (oc2S_seq.dual.normMix, all analysis priors)
# Part 2 – Average OC      (avgoc2S_seq2.normMix, all prior combinations)
# Part 3 – Sequential EUII (Held et al. 2025, eq. 15)
#
# Results saved to Output/SEQ_cond.RDS, Output/SEQ_avg.RDS, Output/SEQ_euii.RDS
#
# NOTE: Part 2 is computationally intensive (nested numerical integration).
# =============================================================================

library(RBesT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
source("scripts/prior_on_control/00.functions.R")


# ---------------------------------------------------------------------------
# Design parameters  (same as fixed design)
# ---------------------------------------------------------------------------
sigma      <- 88
delta_MCID <- 50
mu_A       <- -50

# Sequential sample sizes (cumulative per arm)
# Interim: n_T = 20, n_C = 10  (halfway, same 2:1 allocation)
# Final:   n_T = 40, n_C = 20  (same total as fixed design)
n1_seq <- c(20, 40)   # treatment arm cumulative sizes
n2_seq <- c(10, 20)   # control arm cumulative sizes
n_interim <- n1_seq[1] + n2_seq[1]   # 30
n_final   <- n1_seq[2] + n2_seq[2]   # 60


# ---------------------------------------------------------------------------
# Analysis priors  (identical to 2.1_FIXED.R)
# ---------------------------------------------------------------------------
p_MAP <- mixnorm(
  c(0.4848, -52.457, 21.154), c(0.4598, -47.465, 7.843), c(0.0554, -50.355, 48.164),
  sigma = sigma, param = "ms"
)
p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")

analysis_priors <- list(
  "MAP"            = p_MAP,
  "Robust (w=0.2)" = robustify(p_MAP, 0.2, mean = -50),
  "Robust (w=0.4)" = robustify(p_MAP, 0.4, mean = -50),
  "Robust (w=0.6)" = robustify(p_MAP, 0.6, mean = -50),
  "Robust (w=0.8)" = robustify(p_MAP, 0.8, mean = -50),
  "Vague"          = p_vague
)
prior_names <- names(analysis_priors)


# ---------------------------------------------------------------------------
# Design priors  (identical to 2.1_FIXED.R)
# ---------------------------------------------------------------------------
design_priors <- list(
  "Dirac (-50)"        = mixnorm(c(1, -50, 1e-16), sigma = sigma, param = "ms"),
  "MAP"                = mixnorm(c(0.4848, -52.457, 21.154), c(0.4598, -47.465, 7.843), c(0.0554, -50.355, 48.164),
                                 sigma = sigma, param = "ms"),
  "Skeptical (-90)"    = mixnorm(c(1, -90, 17.6), sigma = sigma, param = "ms"),
  "Misspecified (-10)" = mixnorm(c(1, -10, 25), sigma = sigma, param = "ms")
)
dprior_names <- names(design_priors)


# ---------------------------------------------------------------------------
# Decision criteria  (Gsponer et al. 2014 / manuscript §Decision Rules)
# ---------------------------------------------------------------------------
# Efficacy (both stages): dual criterion
#   Pr(δ > 0  | data) > 0.95  AND  Pr(δ > 50 | data) > 0.50
# qc is in the scale of θ_T − θ_C = −δ, so δ-thresholds flip sign.
sign.crit <- decision2S(pc = 0.95, qc =   0, lower.tail = TRUE) &
             decision2S(pc = 0.50, qc = -50, lower.tail = TRUE)

# Futility (stage 1 only): Pr(δ < 40 | data) ≥ 0.90
#   ↔  Pr(δ > 40 | data) < 0.10
#   ↔  Pr(θ_T − θ_C < −40 | data) < 0.10
fut.crit  <- decision2S(pc = 0.10, qc = -40, lower.tail = FALSE)

# decisions[[1]]: stage 1 — efficacy + futility
# decisions[[2]]: stage 2 — efficacy only
decisions <- list(
  list(success = sign.crit, futility = fut.crit),
  sign.crit
)


# ---------------------------------------------------------------------------
# Shared theme
# ---------------------------------------------------------------------------
my_theme <- theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())


# =============================================================================
# PART 1 — Conditional OC
# =============================================================================
# oc2S_seq.dual.normMix returns a flat data.frame per (theta_T, theta_C):
#   P_Stop_Eff_Stg1, P_Stop_Fut_Stg1, P_Continue_Stg1, P_Success_Stg2,
#   Total_Power, EN_Trt, EN_Pbo,
#   EN_Trt_Rej, EN_Pbo_Rej, EN_Trt_NoRej, EN_Pbo_NoRej

theta_c <- seq(-110, 10, by = 2)

cat("Building conditional sequential OC functions...\n")

oc_funs_seq <- lapply(analysis_priors, function(prior_C) {
  oc2S_seq.dual.normMix(
    prior1    = p_vague,
    prior2    = prior_C,
    n1        = n1_seq,
    n2        = n2_seq,
    decisions = decisions
  )
})

df_cond_seq <- bind_rows(lapply(prior_names, function(pname) {
  f <- oc_funs_seq[[pname]]
  bind_rows(lapply(theta_c, function(tc) {
    r_t1e   <- f(tc,              tc)   # delta = 0
    r_power <- f(tc - delta_MCID, tc)   # delta = MCID
    data.frame(
      Prior            = pname,
      theta_C          = tc,
      # T1E quantities (delta = 0)
      T1E              = r_t1e$Total_Power,
      P_Eff1_T1E       = r_t1e$P_Stop_Eff_Stg1,
      P_Fut1_T1E       = r_t1e$P_Stop_Fut_Stg1,
      EN_T1E           = r_t1e$EN_Trt + r_t1e$EN_Pbo,
      # Power quantities (delta = MCID)
      Power            = r_power$Total_Power,
      P_Eff1_Power     = r_power$P_Stop_Eff_Stg1,
      P_Fut1_Power     = r_power$P_Stop_Fut_Stg1,
      EN_Power         = r_power$EN_Trt + r_power$EN_Pbo,
      EN_Succ          = r_power$EN_Trt_Rej + r_power$EN_Pbo_Rej,
      EN_Fail          = r_power$EN_Trt_NoRej + r_power$EN_Pbo_NoRej
    )
  }))
})) |> mutate(Prior = factor(Prior, levels = prior_names))

# --- Left: Conditional T1E ---
p_seq_T1E <- ggplot(df_cond_seq,
                    aes(x = theta_C, y = T1E,
                        color = Prior, linetype = Prior == "Vague")) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  geom_vline(xintercept = mu_A, linetype = "dotted",
             color = "gray", linewidth = 0.8) +
  scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                        guide = "none") +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = expression("True Control Mean (" * theta[C] * ")"),
       y = "Conditional Type I Error") +
  my_theme

# --- Right: Conditional Power ---
p_seq_Power <- ggplot(df_cond_seq,
                      aes(x = theta_C, y = Power,
                          color = Prior, linetype = Prior == "Vague")) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = mu_A, linetype = "dotted",
             color = "gray", linewidth = 0.8) +
  scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                        guide = "none") +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = expression("True Control Mean (" * theta[C] * ")"),
       y = bquote("Conditional Power (" * delta * " = 50)")) +
  my_theme

# --- Centre: Conditional E[N] at MCID (early stopping benefit) ---
p_seq_EN <- ggplot(df_cond_seq,
                   aes(x = theta_C, y = EN_Power,
                       color = Prior, linetype = Prior == "Vague")) +
  geom_line(linewidth = 1.2) +
  geom_hline(yintercept = n_final, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  geom_vline(xintercept = mu_A, linetype = "dotted",
             color = "gray", linewidth = 0.8) +
  scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                        guide = "none") +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = expression("True Control Mean (" * theta[C] * ")"),
       y = bquote("E[N] (" * delta * " = 50)")) +
  my_theme

grid.arrange(p_seq_T1E, p_seq_Power, p_seq_EN, ncol = 3)

saveRDS(list(df_cond_seq = df_cond_seq), file = "Output/SEQ_cond.RDS")
cat("Sequential conditional results saved to Output/SEQ_cond.RDS\n")


# =============================================================================
# PART 2 — Average OC
# =============================================================================
# avgoc2S_seq2.normMix returns, per (delta_new, design_prior):
#   avg_Power, avg_Fut, avg_EN,
#   avg_EN_Succ, avg_EN_Fail,   (E_D[N^+ | H_i], E_D[N^- | H_i])
#   avg_inv_EN_Succ, avg_inv_EN_Fail  (E_D[1/N^+ | H_i], E_D[1/N^- | H_i])
#
# Convention: delta_new = theta_T - theta_C = -delta (lower is better).
# avgT1E = avg_Power at delta_new = 0.

delta_values <- seq(0, 100, by = 10)

cat("Computing average sequential OC — 6 x 4 analysis x design prior combinations...\n")
cat("(This is computationally intensive; expect several minutes.)\n")

rows_avg_seq <- lapply(prior_names, function(aname) {
  cat(" Analysis prior:", aname, "\n")
  avg_fun <- avgoc2S_seq2.normMix(
    prior1        = p_vague,
    prior2        = analysis_priors[[aname]],
    n1            = n1_seq,
    n2            = n2_seq,
    decisions     = decisions,
    delta         = 0,
    design_prior2 = design_priors[[1]]
  )
  lapply(dprior_names, function(dname) {
    dp <- design_priors[[dname]]
    bind_rows(lapply(delta_values, function(dv) {
      res <- avg_fun(delta_new = -dv, design_prior2_new = dp)
      data.frame(
        Analysis_Prior  = aname,
        Design_Prior    = dname,
        delta           = dv,
        avg_Power       = res$avg_Power,
        avg_Fut         = res$avg_Fut,
        avg_EN          = res$avg_EN,
        avg_EN_Succ     = res$avg_EN_Succ,
        avg_EN_Fail     = res$avg_EN_Fail,
        avg_inv_EN_Succ = res$avg_inv_EN_Succ,
        avg_inv_EN_Fail = res$avg_inv_EN_Fail
      )
    }))
  })
})

df_avg_seq <- bind_rows(unlist(rows_avg_seq, recursive = FALSE)) |>
  mutate(
    Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
    Design_Prior   = factor(Design_Prior,   levels = dprior_names)
  )

# avgT1E = avg_Power at delta = 0 (same for all delta rows within a prior pair)
df_avg_seq <- df_avg_seq |>
  group_by(Analysis_Prior, Design_Prior) |>
  mutate(
    avgT1E          = avg_Power[delta == 0],
    avgT1E_inv_Succ = avg_inv_EN_Succ[delta == 0],
    avgT1E_inv_Fail = avg_inv_EN_Fail[delta == 0]
  ) |>
  ungroup()

# --- Left panel: Average T1E ---
p_avg_seq_T1E <- ggplot(
  filter(df_avg_seq, delta == delta_MCID),
  aes(x = Design_Prior, y = avgT1E,
      color = Analysis_Prior, group = Analysis_Prior,
      linetype = Analysis_Prior == "Vague")) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1.0) +
  geom_hline(yintercept = 0.05, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                        guide = "none") +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = "Design Prior", y = "Average Type I Error") +
  my_theme

# --- Right panel: Average Power at MCID ---
p_avg_seq_Power <- ggplot(
  filter(df_avg_seq, delta == delta_MCID),
  aes(x = Design_Prior, y = avg_Power,
      color = Analysis_Prior, group = Analysis_Prior,
      linetype = Analysis_Prior == "Vague")) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1.0) +
  scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                        guide = "none") +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = "Design Prior",
       y = bquote("Average Power (" * delta * " = 50)")) +
  my_theme

grid.arrange(p_avg_seq_T1E, p_avg_seq_Power, ncol = 2)

saveRDS(list(df_avg_seq = df_avg_seq), file = "Output/SEQ_avg.RDS")
cat("Sequential average results saved to Output/SEQ_avg.RDS\n")


# =============================================================================
# PART 3 — Sequential EUII  (Held et al. 2025, eq. 15)
# =============================================================================
# EUII_seq = (LR+)^{E(1/N+)} / (LR-)^{E(1/N-)}
#
# where LR+ = avgPower / avgT1E,  LR- = (1-avgPower) / (1-avgT1E)
#
# E(1/N+) = E_D[1/N+|H0]*Pr(H0|+) + E_D[1/N+|H1]*Pr(H1|+)
# E(1/N-) = E_D[1/N-|H0]*Pr(H0|-) + E_D[1/N-|H1]*Pr(H1|-)
#
# Pr(H1|+)  = LR+ * prior_H1 / (LR+ * prior_H1 + (1-prior_H1))
# Pr(H1|-)  = LR- * prior_H1 / (LR- * prior_H1 + (1-prior_H1))
#
# E_D[1/N+|H0] = avg_inv_EN_Succ at delta = 0
# E_D[1/N+|H1] = avg_inv_EN_Succ at delta = MCID
# E_D[1/N-|H0] = avg_inv_EN_Fail at delta = 0
# E_D[1/N-|H1] = avg_inv_EN_Fail at delta = MCID

prior_H1 <- 0.5   # prior probability of the alternative hypothesis

omega_map <- c(
  "MAP"            = 0.0,
  "Robust (w=0.2)" = 0.2,
  "Robust (w=0.4)" = 0.4,
  "Robust (w=0.6)" = 0.6,
  "Robust (w=0.8)" = 0.8,
  "Vague"          = 1.0
)

df_euii_seq <- df_avg_seq |>
  filter(delta == delta_MCID) |>
  mutate(
    omega = omega_map[as.character(Analysis_Prior)],

    # Likelihood ratios
    LR_pos = pmax(avg_Power / pmax(avgT1E,         1e-8), 1e-8),
    LR_neg = pmax((1 - avg_Power) / pmax(1 - avgT1E, 1e-8), 1e-8),

    # Posterior Pr(H1 | outcome)
    Pr_H1_pos = LR_pos * prior_H1 / (LR_pos * prior_H1 + (1 - prior_H1)),
    Pr_H0_pos = 1 - Pr_H1_pos,
    Pr_H1_neg = LR_neg * prior_H1 / (LR_neg * prior_H1 + (1 - prior_H1)),
    Pr_H0_neg = 1 - Pr_H1_neg,

    # Mixed expected inverse sample sizes
    E_inv_N_pos = avg_inv_EN_Succ * Pr_H1_pos + avgT1E_inv_Succ * Pr_H0_pos,
    E_inv_N_neg = avg_inv_EN_Fail * Pr_H1_neg + avgT1E_inv_Fail * Pr_H0_neg,

    # Sequential EUII
    log_EUII_seq = E_inv_N_pos * log(LR_pos) - E_inv_N_neg * log(LR_neg),
    EUII_seq     = exp(log_EUII_seq)
  )

# --- Plot A: EUII_seq vs omega at delta = MCID (one line per design prior) ---
euii_vague_seq <- df_euii_seq |>
  filter(Analysis_Prior == "Vague") |>
  pull(EUII_seq) |>
  unique()

p_EUII_seq_omega <- ggplot(df_euii_seq,
                           aes(x = omega, y = EUII_seq,
                               color = Design_Prior, group = Design_Prior)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = euii_vague_seq, linetype = "dotted",
             color = "black", linewidth = 0.6) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = expression("Robustness weight (" * omega * ")"),
       y = expression("EUII"[seq] * " (" * delta * " = 50)")) +
  my_theme

print(p_EUII_seq_omega)

# --- Plot B: EUII_seq vs delta (all analysis x design prior, re-evaluated) ---
# Re-compute for all delta values (not just MCID)
df_euii_seq_full <- df_avg_seq |>
  mutate(
    omega = omega_map[as.character(Analysis_Prior)],
    LR_pos = pmax(avg_Power / pmax(avgT1E,         1e-8), 1e-8),
    LR_neg = pmax((1 - avg_Power) / pmax(1 - avgT1E, 1e-8), 1e-8),
    Pr_H1_pos = LR_pos * prior_H1 / (LR_pos * prior_H1 + (1 - prior_H1)),
    Pr_H0_pos = 1 - Pr_H1_pos,
    Pr_H1_neg = LR_neg * prior_H1 / (LR_neg * prior_H1 + (1 - prior_H1)),
    Pr_H0_neg = 1 - Pr_H1_neg,
    E_inv_N_pos  = avg_inv_EN_Succ * Pr_H1_pos + avgT1E_inv_Succ * Pr_H0_pos,
    E_inv_N_neg  = avg_inv_EN_Fail * Pr_H1_neg + avgT1E_inv_Fail * Pr_H0_neg,
    log_EUII_seq = E_inv_N_pos * log(LR_pos) - E_inv_N_neg * log(LR_neg),
    EUII_seq     = exp(log_EUII_seq)
  )

p_EUII_seq_delta <- ggplot(df_euii_seq_full,
                           aes(x = delta, y = EUII_seq,
                               color = Analysis_Prior,
                               linetype = Analysis_Prior == "Vague")) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~Design_Prior, nrow = 1) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                        guide = "none") +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = expression(delta),
       y = expression("EUII"[seq])) +
  my_theme

print(p_EUII_seq_delta)

saveRDS(list(df_euii_seq = df_euii_seq, df_euii_seq_full = df_euii_seq_full),
        file = "Output/SEQ_euii.RDS")
cat("Sequential EUII results saved to Output/SEQ_euii.RDS\n")
