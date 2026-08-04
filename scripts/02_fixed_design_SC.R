# =============================================================================
# Fixed Design — Conditional and Average Operating Characteristics
# =============================================================================
# Part 1 – Conditional T1E and Power 
# Part 2 – Average T1E and Power (all analysis prior × design prior combinations)
## =============================================================================

library(RBesT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
source("scripts/00_functions_exact.R")


# ---------------------------------------------------------------------------
# Design parameters
# ---------------------------------------------------------------------------
sigma      <- 88
n_T        <- 40
n_C        <- 20
delta_MCID <- 60 #I chose to plot Power using this delta
mu_A       <- -50   # prior mean — reference vline


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
  "Dirac (-50)"    = mixnorm(c(1, -50, 1e-16), sigma = sigma, param = "ms"),
  "MAP"               = mixnorm(c(0.4848, -52.457, 21.154), c(0.4598, -47.465, 7.843), c(0.0554, -50.355, 48.164),
                                sigma = sigma, param = "ms"),
  "Skeptical (-90)" = mixnorm(c(1, -90, 17.6),  sigma = sigma, param = "ms"),
  "Misspecified (-10)"  = mixnorm(c(1, -10, 25),   sigma = sigma, param = "ms")
)
dprior_names <- names(design_priors)


# ---------------------------------------------------------------------------
# Decision criterion (Sign only)
# ---------------------------------------------------------------------------
sign.crit <- decision2S(pc = 0.95, qc = 0, lower.tail = TRUE)


# ---------------------------------------------------------------------------
# Shared theme
# ---------------------------------------------------------------------------
my_theme <- theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())


# =============================================================================
# PART 1 — Conditional T1E and Power
# =============================================================================
theta_c <- seq(-110, 10, by = 2)

oc_funs <- lapply(analysis_priors, function(prior_C) {
  oc2S(prior1 = p_vague, prior2 = prior_C,
       n1 = n_T, n2 = n_C, decision = sign.crit)
})

df_cond <- bind_rows(lapply(prior_names, function(pname) {
  f <- oc_funs[[pname]]
  data.frame(
    Prior   = pname,
    theta_C = theta_c,
    T1E     = f(theta_c, theta_c),
    Power   = f(theta_c - delta_MCID, theta_c)
  )
})) |> mutate(Prior = factor(Prior, levels = prior_names))

# --- Left panel: Conditional T1E ---
p_cond_T1E <- ggplot(df_cond,
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

# --- Right panel: Conditional Power at delta == 60 ---
p_cond_Power <- ggplot(df_cond,
                       aes(x = theta_C, y = Power,
                           color = Prior, linetype = Prior == "Vague")) +
  geom_line(linewidth = 1.2) +
  geom_vline(xintercept = mu_A, linetype = "dotted",
             color = "gray", linewidth = 0.8) +
  scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                        guide = "none") +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = expression("True Control Mean (" * theta[C] * ")"),
       y = bquote("Conditional Power (" * delta * " = 60)")) +
  my_theme

grid.arrange(p_cond_T1E, p_cond_Power, ncol = 2)

saveRDS(list(df_cond = df_cond), file = "Output/fixed_SC_cond.RDS")
cat("Conditional results saved to Output/fixed_SC_cond.RDS\n")


# =============================================================================
# PART 2 — Average T1E and Power (all analysis prior x design prior x delta)
# =============================================================================
# avgoc2S.normMix returns design_fun(delta_new, design_prior2_new).
# theta_T = theta_C + delta_new, so delta_new = -delta (lower is better).
# avgT1E is delta-independent (evaluated at delta_new = 0 only once per pair).

delta_values <- seq(0, 100, by = 10)

rows_avg <- lapply(prior_names, function(aname) {
  avg_fun <- avgoc2S.normMix(
    prior1        = p_vague,
    prior2        = analysis_priors[[aname]],
    n1            = n_T,
    n2            = n_C,
    decision      = sign.crit,
    delta         = 0,
    design_prior2 = design_priors[[1]]
  )
  lapply(dprior_names, function(dname) {
    dp         <- design_priors[[dname]]
    t1e_val    <- avg_fun(delta_new = 0, design_prior2_new = dp)
    bind_rows(lapply(delta_values, function(dv) {
      data.frame(
        Analysis_Prior = aname,
        Design_Prior   = dname,
        delta          = dv,
        avgT1E         = t1e_val, # a column with all t1e identical 
        avgPower       = avg_fun(delta_new = -dv, design_prior2_new = dp) # first of these have delta = 0, so it is avgT1E
      )
    }))
  })
})

# rows_avg is a list of list of dataframes. SO rows_avg[[1]][[2]] recovers the df with analysis prior 1 (MAP)
# and design priro 2 (MAP)
rows_avg[[1]][[2]]
#    Analysis_Prior Design_Prior delta     avgT1E   avgPower
# 1             MAP          MAP     0 0.04999885 0.04999885
# 2             MAP          MAP    10 0.04999885 0.12897861
# 3             MAP          MAP    20 0.04999885 0.27192287
# 4             MAP          MAP    30 0.04999885 0.46851064
# 5             MAP          MAP    40 0.04999885 0.67210689
# 6             MAP          MAP    50 0.04999885 0.83162425
# 7             MAP          MAP    60 0.04999885 0.92778874
# 8             MAP          MAP    70 0.04999885 0.97363421
# 9             MAP          MAP    80 0.04999885 0.99154804
# 10            MAP          MAP    90 0.04999885 0.99753338
# 11            MAP          MAP   100 0.04999885 0.99932471



df_avg <- bind_rows(unlist(rows_avg, recursive = FALSE)) |> # use recursive so it unlist only one layer
  mutate(
    Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
    Design_Prior   = factor(Design_Prior,   levels = dprior_names)
  )

# --- Left panel: Average T1E ---
p_avg_T1E <- ggplot(df_avg,
                    aes(x = Design_Prior, y = avgT1E,
                        color = Analysis_Prior,
                        group = Analysis_Prior,
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

# --- Right panel: Average Power ---
p_avg_Power <- ggplot(df_avg,
                      aes(x = Design_Prior, y = avgPower,
                          color = Analysis_Prior,
                          group = Analysis_Prior,
                          linetype = Analysis_Prior == "Vague")) +
  geom_point(size = 2.5) +
  geom_line(linewidth = 1.0) +
  scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"),
                        guide = "none") +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = "Design Prior",
       y = bquote("Average Power (" * delta * " = 60)")) +
  my_theme

grid.arrange(p_avg_T1E, p_avg_Power, ncol = 2)

saveRDS(list(df_avg = df_avg), file = "Output/fixed_SC_avg.RDS")
cat("Average results saved to Output/fixed_SC_avg.RDS\n")


# =============================================================================
# PART 3 — EUII 
# =============================================================================
# DOR  = (avgPower / (1-avgPower)) / (avgT1E / (1-avgT1E))
# EUII = DOR ^ (1 / n_trial)   where n_trial = n_T + n_C
# At delta = 0: avgPower = avgT1E  =>  DOR = 1  =>  EUII = 1 for all designs.

n_trial <- n_T + n_C

omega_map <- c(
  "MAP"            = 0.0,
  "Robust (w=0.2)" = 0.2,
  "Robust (w=0.4)" = 0.4,
  "Robust (w=0.6)" = 0.6,
  "Robust (w=0.8)" = 0.8,
  "Vague"          = 1.0
)

df_euii <- df_avg |>
  mutate(
    omega    = omega_map[as.character(Analysis_Prior)],
    # log space: qlogis(p) = logit(p) = log(p / (1 - p)), so
    # log_DOR = logit(avgPower) - logit(avgT1E)
    log_DOR  = qlogis(avgPower) - qlogis(avgT1E),
    log_EUII = log_DOR / n_trial,
    EUII     = exp(log_EUII)
  )

# --- Plot A: EUII vs omega at delta = 60 (one line per design prior) ---
df_euii_50 <- filter(df_euii, delta == 60)

euii_vague_60 <- df_euii_50 |>
  filter(Analysis_Prior == "Vague") |>
  pull(EUII) |>
  unique()

p_EUII_omega <- ggplot(df_euii_50,
                       aes(x = omega, y = EUII,
                           color = Design_Prior, group = Design_Prior)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = euii_vague_60, linetype = "dashed",
             color = "black", linewidth = 0.6) +
  scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = expression("Robustness weight (" * omega * ")"),
       y = expression("EUII (" * delta * " = 60)")) +
  my_theme

print(p_EUII_omega)

# --- Plot B: EUII vs delta (facet by design prior, color by analysis prior) ---
p_EUII_delta <- ggplot(df_euii,
                       aes(x = delta, y = EUII,
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
       y = "EUII") +
  my_theme

print(p_EUII_delta)

saveRDS(list(df_euii = df_euii), file = "Output/fixed_SC_euii.RDS")
cat("EUII results saved to Output/fixed_SC_euii.RDS\n")
