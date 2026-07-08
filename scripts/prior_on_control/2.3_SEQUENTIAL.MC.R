# =============================================================================
# Sequential Design — Monte Carlo verification of OC and EUII
# =============================================================================
# Run from the Thesis/ root directory:
#   setwd("/home/menelao/Thesis")
#   source("scripts/prior_on_control/2.3_SEQUENTIAL.MC.R")
#
# Part 1 – Conditional OC  (oc2_seq_mc.normMix, all analysis priors)
# Part 2 – Average OC      (avgoc2_seq_mc.normMix, all prior combinations)
# Part 3 – Sequential EUII (compute_euii, same prior_H1 = 0.5 as 2.2)
#
# Results saved to Output/SEQ_MC_cond.RDS, Output/SEQ_MC_avg.RDS,
#                  Output/SEQ_MC_euii.RDS
#
# NOTE: All three parts are MC-based (n_sim = 100 000 per evaluation).
#       Part 2 runs 6 x 4 x 11 = 264 MC evaluations — allow ~20–30 min.
# =============================================================================

library(RBesT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
source("scripts/prior_on_control/00.functions.MC.R")


# ---------------------------------------------------------------------------
# Design parameters  (identical to 2.2_SEQUENTIAL.R)
# ---------------------------------------------------------------------------
sigma      <- 88
delta_MCID <- 60
mu_A       <- -50

n1_seq <- c(20, 40)   # treatment arm cumulative sizes
n2_seq <- c(10, 20)   # control arm cumulative sizes
n_interim <- n1_seq[1] + n2_seq[1]   # 30
n_final   <- n1_seq[2] + n2_seq[2]   # 60

n_sim_cond <- 200000   # sims per conditional OC point (Part 1)
n_sim_avg  <- 100000   # sims per average OC evaluation (Part 2)


# ---------------------------------------------------------------------------
# Analysis priors  
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
# Design priors  
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
# Decision criteria  
# ---------------------------------------------------------------------------
# Efficacy (both stages): 
# θ_T − θ_C = −delta
#   Pr(delta > 0  | data) > 0.95  AND  Pr(delta > 50 | data) > 0.50
# qc is in the scale of 
sign.crit <- list(
  decision2S(pc = 0.95, qc =   0, lower.tail = TRUE),  
  decision2S(pc = 0.50, qc = -50, lower.tail = TRUE)   
)

# Futility (stage 1 only): Pr(theta_T - theta_C > -40 | data) ≥ 0.90
#
fut.crit  <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)

decisions_mc <- list(
  list(success = sign.crit, futility = fut.crit),  # stage 1
  list(success = sign.crit, futility = NULL)        # stage 2 (final)
)


# ---------------------------------------------------------------------------
# Shared theme
# ---------------------------------------------------------------------------
my_theme <- theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())


# =============================================================================
# PART 1 — Conditional OC (MC)
# =============================================================================
# Each point requires two MC evaluations (delta=0 and delta=MCID).

theta_c_mc <- seq(-110, 10, by = 5)


df_cond_seq_mc <- bind_rows(lapply(prior_names, function(pname) {
  cat(" Analysis prior:", pname, "\n")
  bind_rows(lapply(seq_along(theta_c_mc), function(i) {
    tc <- theta_c_mc[i]
    r_t1e <- oc2_seq_mc.normMix(
      theta_1        = tc,
      theta_2        = tc,
      prior_1        = p_vague,
      prior_2        = analysis_priors[[pname]],
      n1_seq         = n1_seq,
      n2_seq         = n2_seq,
      decisions_list = decisions_mc,
      n_sim          = n_sim_cond,
      seed           = 42 + i
    )
    r_pwr <- oc2_seq_mc.normMix(
      theta_1        = tc - delta_MCID,
      theta_2        = tc,
      prior_1        = p_vague,
      prior_2        = analysis_priors[[pname]],
      n1_seq         = n1_seq,
      n2_seq         = n2_seq,
      decisions_list = decisions_mc,
      n_sim          = n_sim_cond,
      seed           = 142 + i
    )
    data.frame(
      Prior        = pname,
      theta_C      = tc,
      T1E          = r_t1e$Overall[["Power"]],
      P_Eff1_T1E   = r_t1e$Per_Stage$P_Succ[1],
      P_Fut1_T1E   = r_t1e$Per_Stage$P_Fut[1],
      Power        = r_pwr$Overall[["Power"]],
      P_Eff1_Power = r_pwr$Per_Stage$P_Succ[1],
      P_Fut1_Power = r_pwr$Per_Stage$P_Fut[1],
      EN_Power     = r_pwr$Overall[["EN_t"]] + r_pwr$Overall[["EN_c"]]
    )
  }))
})) |> mutate(Prior = factor(Prior, levels = prior_names))

p_mc_T1E <- ggplot(df_cond_seq_mc,
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
       y = "Conditional Type I Error (MC)") +
  my_theme

p_mc_Power <- ggplot(df_cond_seq_mc,
                     aes(x = theta_C, y = Power,
                         color = Prior, linetype = Prior == "Vague")) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = mu_A, linetype = "dotted",
             color = "gray", linewidth = 0.8) +
  scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                        guide = "none") +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = expression("True Control Mean (" * theta[C] * ")"),
       y = bquote("Conditional Power (MC, " * delta * " = 60)")) +
  my_theme

grid.arrange(p_mc_T1E, p_mc_Power, ncol = 2)

saveRDS(list(df_cond_seq_mc = df_cond_seq_mc), file = "Output/SEQ_MC_cond.RDS")
cat("MC conditional results saved to Output/SEQ_MC_cond.RDS\n")


# =============================================================================
# PART 2 — Average OC (MC)
# =============================================================================
# avgoc2_seq_mc.normMix is called once per (analysis_prior, design_prior) pair
# with all delta values.  The returned list is named "delta.0", "delta.10", etc.

delta_values <- c(0, seq(10, 100, by = 10))

cat("\nPart 2: Average OC via MC (n_sim =", n_sim_avg, ") ...\n")
cat("Running", length(prior_names), "x", length(dprior_names), "combinations.\n")

# Store raw MC results (needed by compute_euii)
mc_results_all <- list()

rows_avg_mc <- lapply(prior_names, function(aname) {
  cat(" Analysis prior:", aname, "\n")
  lapply(dprior_names, function(dname) {
    cat("   Design prior:", dname, "\n")
    res_mc <- avgoc2_seq_mc.normMix(
      prior_1        = p_vague,
      prior_2        = analysis_priors[[aname]],
      n1_seq         = n1_seq,
      n2_seq         = n2_seq,
      decisions_list = decisions_mc,
      delta          = delta_values,
      design_prior_c = design_priors[[dname]],
      n_sim          = n_sim_avg,
      seed           = 123
    )
    # Store for EUII computation
    mc_results_all[[paste(aname, dname, sep = "||")]] <<- res_mc

    bind_rows(lapply(names(res_mc), function(nm) {
      ov <- res_mc[[nm]]$Overall
      ps <- res_mc[[nm]]$Per_Stage
      dv <- as.numeric(sub("delta\\.", "", nm))
      pwr <- ov[["Power"]]
      # E_D[1/N^+ | H_i] * avgPower (numerator) = sum_k P_Succ_k / n_k
      num_inv_EN_Succ <- sum(ps$P_Succ / (ps$N_Trt + ps$N_Ctrl))
      # E_D[1/N^- | H_i] * (1-avgPower) (numerator)
      K <- nrow(ps)
      nonsig_prob <- c(ps$P_Fut[-K], 1 - ps$Cum_P_Succ[K] - ps$Cum_P_Fut[K])
      num_inv_EN_Fail <- sum(nonsig_prob / (ps$N_Trt + ps$N_Ctrl))
      data.frame(
        Analysis_Prior  = aname,
        Design_Prior    = dname,
        delta           = dv,
        avg_Power       = pwr,
        avg_Fut         = sum(ps$P_Fut),
        avg_EN          = ov[["EN_t"]] + ov[["EN_c"]],
        avg_EN_Succ     = if (pwr   > 1e-8) (ov[["EN_t_Succ"]] + ov[["EN_c_Succ"]]) else NA_real_,
        avg_EN_Fail     = if (1-pwr > 1e-8) (ov[["EN_t_Fail"]] + ov[["EN_c_Fail"]]) else NA_real_,
        avg_inv_EN_Succ = if (pwr   > 1e-8) num_inv_EN_Succ / pwr   else NA_real_,
        avg_inv_EN_Fail = if (1-pwr > 1e-8) num_inv_EN_Fail / (1 - pwr) else NA_real_
      )
    }))
  })
})

df_avg_seq_mc <- bind_rows(unlist(rows_avg_mc, recursive = FALSE)) |>
  mutate(
    Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
    Design_Prior   = factor(Design_Prior,   levels = dprior_names)
  ) |>
  group_by(Analysis_Prior, Design_Prior) |>
  mutate(
    avgT1E          = avg_Power[delta == 0],
    avgT1E_inv_Succ = avg_inv_EN_Succ[delta == 0],
    avgT1E_inv_Fail = avg_inv_EN_Fail[delta == 0]
  ) |>
  ungroup()

p_avg_mc_T1E <- ggplot(
  filter(df_avg_seq_mc, delta == delta_MCID),
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
  labs(x = "Design Prior", y = "Average Type I Error (MC)") +
  my_theme

p_avg_mc_Power <- ggplot(
  filter(df_avg_seq_mc, delta == delta_MCID),
  aes(x = Design_Prior, y = avg_Power,
      color = Analysis_Prior, group = Analysis_Prior,
      linetype = Analysis_Prior == "Vague")) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1.0) +
  scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                        guide = "none") +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = "Design Prior",
       y = bquote("Average Power (MC, " * delta * " = 60)")) +
  my_theme

grid.arrange(p_avg_mc_T1E, p_avg_mc_Power, ncol = 2)

saveRDS(list(df_avg_seq_mc = df_avg_seq_mc, mc_results_all = mc_results_all),
        file = "Output/SEQ_MC_avg.RDS")
cat("MC average results saved to Output/SEQ_MC_avg.RDS\n")


# =============================================================================
# PART 3 — Sequential EUII via compute_euii
# =============================================================================
# compute_euii correctly uses E[1/N] (not 1/E[N]).
# Called once per (analysis_prior, design_prior) on the stored MC result list.

prior_H1 <- 0.5

cat("\nPart 3: Sequential EUII via compute_euii ...\n")

omega_map <- c(
  "MAP"            = 0.0,
  "Robust (w=0.2)" = 0.2,
  "Robust (w=0.4)" = 0.4,
  "Robust (w=0.6)" = 0.6,
  "Robust (w=0.8)" = 0.8,
  "Vague"          = 1.0
)

euii_rows <- lapply(prior_names, function(aname) {
  lapply(dprior_names, function(dname) {
    key    <- paste(aname, dname, sep = "||")
    res_mc <- mc_results_all[[key]]
    euii   <- compute_euii(res_mc, prior_H1 = prior_H1)[[as.character(prior_H1)]]
    euii |>
      mutate(
        Analysis_Prior = aname,
        Design_Prior   = dname,
        omega          = omega_map[[aname]]
      )
  })
})

df_euii_seq_mc <- bind_rows(unlist(euii_rows, recursive = FALSE)) |>
  mutate(
    Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
    Design_Prior   = factor(Design_Prior,   levels = dprior_names)
  )

# --- Plot A: EUII vs omega at delta = MCID ---
df_euii_mc_50 <- filter(df_euii_seq_mc, Delta == delta_MCID)

euii_vague_mc <- df_euii_mc_50 |>
  filter(Analysis_Prior == "Vague") |>
  pull(EUII) |>
  unique()

p_EUII_mc_omega <- ggplot(df_euii_mc_50,
                          aes(x = omega, y = EUII,
                              color = Design_Prior, group = Design_Prior)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = euii_vague_mc, linetype = "dotted",
             color = "black", linewidth = 0.6) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = expression("Robustness weight (" * omega * ")"),
       y = expression("EUII"[seq] * " (MC, " * delta * " = 60)")) +
  my_theme

print(p_EUII_mc_omega)

# --- Plot B: EUII vs delta ---
p_EUII_mc_delta <- ggplot(df_euii_seq_mc,
                           aes(x = Delta, y = EUII,
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
       y = expression("EUII"[seq] * " (MC)")) +
  my_theme

print(p_EUII_mc_delta)

saveRDS(list(df_euii_seq_mc = df_euii_seq_mc), file = "Output/SEQ_MC_euii.RDS")
cat("MC EUII results saved to Output/SEQ_MC_euii.RDS\n")


# =============================================================================
# PART 4 — EUII Ratio: Sequential vs Fixed Dual Criterion
# =============================================================================
# Plots EUII_seq / EUII_dual_fixed for two designs side-by-side:
#   "With Futility"    — decisions_mc (Part 2/3, already computed)
#   "Without Futility" — efficacy-only stopping at both stages
#
# Re-uses mc_results_all from Part 2; runs a second MC sweep for no-futility.
# Loads Output/DUAL_euii.RDS for the fixed dual criterion baseline.
# =============================================================================

# --- 4a: Decision list without futility stopping ----------------------------
decisions_mc_no_fut <- list(
  list(success = sign.crit, futility = NULL),   # stage 1: efficacy only
  list(success = sign.crit, futility = NULL)    # stage 2: final
)

cat("\nPart 4: Average MC without futility stopping (n_sim =", n_sim_avg, ") ...\n")

mc_results_no_fut <- list()

lapply(prior_names, function(aname) {
  cat(" Analysis prior:", aname, "\n")
  lapply(dprior_names, function(dname) {
    cat("   Design prior:", dname, "\n")
    res_mc <- avgoc2_seq_mc.normMix(
      prior_1        = p_vague,
      prior_2        = analysis_priors[[aname]],
      n1_seq         = n1_seq,
      n2_seq         = n2_seq,
      decisions_list = decisions_mc_no_fut,
      delta          = delta_values,
      design_prior_c = design_priors[[dname]],
      n_sim          = n_sim_avg,
      seed           = 456
    )
    mc_results_no_fut[[paste(aname, dname, sep = "||")]] <<- res_mc
    NULL
  })
  NULL
})

# --- 4b: EUII for no-futility design ----------------------------------------
euii_rows_no_fut <- lapply(prior_names, function(aname) {
  lapply(dprior_names, function(dname) {
    key  <- paste(aname, dname, sep = "||")
    euii <- compute_euii(mc_results_no_fut[[key]],
                         prior_H1 = prior_H1)[[as.character(prior_H1)]]
    euii |> mutate(Analysis_Prior = aname,
                   Design_Prior   = dname,
                   omega          = omega_map[[aname]])
  })
})

df_euii_no_fut <- bind_rows(unlist(euii_rows_no_fut, recursive = FALSE)) |>
  mutate(
    Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
    Design_Prior   = factor(Design_Prior,   levels = dprior_names)
  )

# --- 4c: Load fixed dual criterion EUII baseline ----------------------------
fixed_dual         <- readRDS("Output/DUAL_euii.RDS")
df_euii_fixed_dual <- fixed_dual$df_euii |>
  mutate(
    Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
    Design_Prior   = factor(Design_Prior,   levels = dprior_names)
  ) |>
  select(Analysis_Prior, Design_Prior, delta, omega, EUII_fixed = EUII)

# --- 4d: Build combined ratio dataframe -------------------------------------
make_ratio_df <- function(df_euii_seq, stopping_label) {
  df_euii_seq |>
    rename(delta = Delta) |>
    select(Analysis_Prior, Design_Prior, delta, omega, EUII_seq = EUII) |>
    left_join(df_euii_fixed_dual,
              by = c("Analysis_Prior", "Design_Prior", "delta", "omega")) |>
    mutate(EUII_ratio = EUII_seq / EUII_fixed,
           Stopping   = stopping_label)
}

df_ratio_combined <- bind_rows(
  make_ratio_df(df_euii_seq_mc, "With Futility"),
  make_ratio_df(df_euii_no_fut, "Without Futility")
) |>
  mutate(Stopping = factor(Stopping,
                           levels = c("With Futility", "Without Futility")))

# --- 4e: Plot A — ratio vs omega at delta = delta_MCID ----------------------
p_ratio_omega <- ggplot(
  filter(df_ratio_combined, delta == delta_MCID),
  aes(x = omega, y = EUII_ratio,
      color = Design_Prior,
      group = interaction(Design_Prior, Stopping),
      linetype = Stopping)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_linetype_manual(
    values = c("With Futility" = "solid", "Without Futility" = "dashed")) +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = expression("Robustness weight (" * omega * ")"),
       y = bquote("EUII"["Seq"] / "EUII"["Dual"] ~ "(" * delta * " = 60)"),
       linetype = "Futility stopping") +
  my_theme

print(p_ratio_omega)

# --- 4f: Plot B — ratio vs delta (facet by design prior) --------------------
p_ratio_delta <- ggplot(
  df_ratio_combined,
  aes(x = delta, y = EUII_ratio,
      color = Analysis_Prior,
      group = interaction(Analysis_Prior, Stopping),
      linetype = Stopping)) +
  geom_line(linewidth = 1.2) +
  facet_wrap(~Design_Prior, nrow = 1) +
  geom_hline(yintercept = 1, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  scale_linetype_manual(
    values = c("With Futility" = "solid", "Without Futility" = "dashed")) +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = expression(delta),
       y = bquote("EUII"["Seq"] / "EUII"["Dual"]),
       linetype = "Futility stopping") +
  my_theme

print(p_ratio_delta)

saveRDS(list(df_ratio_combined = df_ratio_combined,
             df_euii_no_fut    = df_euii_no_fut),
        file = "Output/SEQ_vs_DUAL_euii.RDS")
cat("EUII ratio saved to Output/SEQ_vs_DUAL_euii.RDS\n")
