# =============================================================================
# Closed-Form Running Example 
# =============================================================================
# Reproduction of the 8 R chunks of manuscript.Rnw that use the 
# closed-form running example 
# =============================================================================

# Shared colour palette, reused across every figure in the manuscript so that
# SC/DC and Average/Conditional always match
col_sc   <- "#0072B2"
col_dc   <- "#D55E00"
col_avg  <- "#E69F00"
col_cond <- "#56B4E9"


# =============================================================================
# Chunk: plot_t1e
# =============================================================================
library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)

# --- Trial Parameters  ---
sigma <- 88; n_C <- 20; n_T <- 40; n_A <- 20; N_total <- n_C + n_T
v_C <- sigma^2 / n_C; v_T <- sigma^2 / n_T; sigma_A2 <- sigma^2 / n_A; mu_A <- -50
q <- 0; z_p <- qnorm(0.95)
W <- sigma_A2 / (sigma_A2 + v_C)
V_delta <- W * v_C + v_T
V_samp <- W^2 * v_C + v_T
k_SC <- q + z_p * sqrt(V_delta)   # SC (significance) threshold

# rejection probabilities for a generic threshold k
# (SC: k = k_SC; dual criterion: k = max(DV, k_SC))
condPower <- function(theta_C_star, delta, k)
  1 - pnorm((k - ((1 - W) * (mu_A - theta_C_star) + delta)) / sqrt(V_samp))

avgPower <- function(mu_D, sigma_D2, delta, k)
  1 - pnorm((k - ((1 - W) * (mu_A - mu_D) + delta)) / sqrt(V_samp + (1 - W)^2 * sigma_D2))

calc_EUII <- function(power, t1e, N) {
  eps <- 1e-10
  p <- pmax(pmin(power, 1 - eps), eps); t <- pmax(pmin(t1e, 1 - eps), eps)

  log_dor <- qlogis(p) - qlogis(t)
  return(exp(log_dor / N))
}

my_theme <- theme_minimal() + theme(legend.position = "bottom", legend.title = element_blank())
my_x_scale <- scale_x_continuous(trans = "log1p", breaks = c(0, 10, 100, 1000, 10000), labels = c("0\n(Cond.)", "10", "100", "1K", "10K"))

# --- T1E Data & Plots ---
df_t1e_loc <- data.frame(mu_D = seq(-90, -10, length.out = 100)) |>
  mutate(
    `conditional T1E` = condPower(mu_D, 0, k_SC),
    `average T1E` = avgPower(mu_D, sigma_A2, 0, k_SC)
  ) |>
  pivot_longer(cols = c(`conditional T1E`, `average T1E`), names_to = "Metric", values_to = "Value")

sigmas <- c(0, 10^seq(0, 5, length.out = 99))
df_t1e_scale <- data.frame(sigma_D = sigmas) |>
  mutate(Avg_T1E = avgPower(mu_A, sigma_D^2, 0, k_SC))

p1a <- ggplot(df_t1e_loc, aes(x = mu_D, y = Value, color = Metric)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_vline(xintercept = mu_A, linetype = "dotted", color = "gray", linewidth = 0.8) +
  scale_color_manual(values = c("average T1E" = col_avg, "conditional T1E" = col_cond)) +
  labs(x = expression("True control mean (" * mu["D,C"] * ")"),
       y = "T1E") +
  my_theme

p1b <- ggplot(df_t1e_scale, aes(x = sigma_D, y = Avg_T1E)) +
  geom_line(color = "#E69F00", size = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_vline(xintercept = sqrt(sigma_A2), linetype = "dotted", color = "gray", linewidth = 0.8) +
  my_x_scale +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(x = expression(sigma["D,C"]),
       y = "avgT1E") +
  my_theme

grid.arrange(p1a, p1b, ncol = 2)


# =============================================================================
# Chunk: plot_power
# =============================================================================
# --- Power Data & Plots ---
df_pow_loc <- data.frame(delta = seq(0, 100, length.out = 100)) |>
  mutate(
    `conditional power` = condPower(mu_A, delta, k_SC),
    `average power` = avgPower(mu_A, sigma_A2, delta, k_SC)
  ) |>
  pivot_longer(cols = c(`conditional power`, `average power`), names_to = "Metric", values_to = "Value")

df_pow_scale <- data.frame(sigma_D = sigmas) |>
  mutate(Avg_Power = avgPower(mu_A, sigma_D^2, 60, k_SC))

p2a <- ggplot(df_pow_loc, aes(x = delta, y = Value, color = Metric)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("average power" = col_avg, "conditional power" = col_cond)) +
  labs(x = expression("True effect (" * delta * ")"),
       y = "Power") +
  my_theme

p2b <- ggplot(df_pow_scale, aes(x = sigma_D, y = Avg_Power)) +
  geom_line(color = "#E69F00", size = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color="black", linewidth = 0.6) +
  geom_vline(xintercept = sqrt(sigma_A2), linetype = "dotted", color = "gray", linewidth = 0.8) +
  my_x_scale +
  coord_cartesian(ylim = c(0.4, 1.0)) +
  labs(x = expression(sigma["D,C"]),
       y = "avgPower") +
  my_theme

grid.arrange(p2a, p2b, ncol = 2)


# =============================================================================
# Chunk: plot_avgpower_nc
# =============================================================================
avgPower_ess <- function(delta, n_C, n_A, r_T, sigma, z_p, q = 0) {
  n_T <- r_T * n_C
  1 - pnorm(z_p - (delta - q) / (sigma * sqrt(1 / (n_C + n_A) + 1 / n_T)))
}

n_C_seq <- seq(4, 40, by = 0.5)
n_A_vals <- c(0, 10, 20, 30)
df_ss <- do.call(rbind, lapply(n_A_vals, function(na) {
  data.frame(
    n_C   = n_C_seq,
    power = avgPower_ess(delta = 60, n_C = n_C_seq, n_A = na,
                         r_T = 2, sigma = 88, z_p = qnorm(0.95)),
    n_A   = factor(paste0("n_A = ", na))
  )
}))

#
K_req   <- 88^2 * (qnorm(0.95) + qnorm(0.80))^2 / 60^2 # 13.30
b_req   <- K_req * (2 + 1) - 2 * 20
n_C_req <- (b_req + sqrt(b_req^2 + 4 * 2 * K_req * 20)) / (2 * 2)  # 11.5

legend_labels <- expression(n[A] == 0, n[A] == 10, n[A] == 20, n[A] == 30)

ggplot(df_ss, aes(x = n_C, y = power, colour = n_A)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = n_C_req, linetype = "dotted", colour = "grey40") +
  geom_vline(xintercept = 20,      linetype = "dotted", colour = "grey40") +
  annotate("text", x = n_C_req - 0.6, y = 0.15, label = paste0("n[C] == ", round(n_C_req, 1)),
           hjust = 1, size = 3, colour = "grey30", parse = TRUE) +
  annotate("text", x = 20 + 0.6, y = 0.15, label = "n[C] == 20",
           hjust = 0, size = 3, colour = "grey30", parse = TRUE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_colour_brewer(palette = "Set1", labels = legend_labels) +
  labs(x = expression(n[C]), y = "average power",
       colour = NULL, linetype = NULL) +
  my_theme +
  theme(legend.position = "right")
 

# =============================================================================
# Chunk: plot_t1e_ratio
# =============================================================================
# Dual-criterion threshold for the running example (k_SC defined in setup)
DV   <- 50
k_DC <- max(DV, k_SC)

df_rt1e_loc <- data.frame(mu_D = seq(-90, -10, length.out = 100)) |>
  mutate(
    `conditional` = condPower(mu_D, 0, k_DC) / condPower(mu_D, 0, k_SC),
    `average`     = avgPower(mu_D, sigma_A2, 0, k_DC) / avgPower(mu_D,  sigma_A2, 0, k_SC)
  ) |>
  pivot_longer(cols = c(`conditional`, `average`), names_to = "Metric", values_to = "Ratio")

df_t1e_abs <- data.frame(mu_D = seq(-90, -10, length.out = 100)) |>
  mutate(
     DC         = avgPower(mu_D, sigma_A2, 0, k_DC),
    SC = avgPower(mu_D, sigma_A2, 0, k_SC)
  ) |>
  pivot_longer(cols = c(DC, SC), names_to = "Criterion", values_to = "avgT1E")

pr1a <- ggplot(df_rt1e_loc, aes(mu_D, Ratio, color = Metric)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = mu_A, linetype = "dotted", color = "gray", linewidth = 0.8) +
  scale_color_manual(values = c("average" = col_avg, "conditional" = col_cond)) +
  labs(x = expression("True control mean (" * mu["D,C"] * ")"),
       y = "T1E ratio (DC / SC)") + my_theme

pr1b <- ggplot(df_t1e_abs, aes(mu_D, avgT1E, color = Criterion)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_vline(xintercept = mu_A, linetype = "dotted", color = "gray", linewidth = 0.8) +
  scale_color_manual(values = c("DC" = col_dc, "SC" = col_sc)) +
  labs(x = expression("True control mean (" * mu["D,C"] * ")"),
       y = "average T1E") + my_theme

grid.arrange(pr1a, pr1b, ncol = 2)


# =============================================================================
# Chunk: plot_power_ratio
# =============================================================================
df_pow_dc <- data.frame(delta = seq(0, 100, length.out = 100)) |>
  mutate(
    `conditional power` = condPower(mu_A, delta, k_DC),
    `average power`     = avgPower(mu_A, sigma_A2, delta, k_DC)
  ) |>
  pivot_longer(cols = c(`conditional power`, `average power`),
               names_to = "Metric", values_to = "Value")

df_pow_abs <- data.frame(delta = seq(0, 100, length.out = 100)) |>
  mutate(
    DC         = avgPower(mu_A, sigma_A2, delta, k_DC),
    SC = avgPower(mu_A, sigma_A2, delta, k_SC)
  ) |>
  pivot_longer(cols = c(DC, SC), names_to = "Criterion", values_to = "avgPower")

pr2a <- ggplot(df_pow_dc, aes(delta, Value, color = Metric)) +
  geom_line(size = 1) +
  geom_vline(xintercept = DV, linetype = "dotted", color = "gray", linewidth = 0.8) +
  scale_color_manual(values = c("average power" = col_avg, "conditional power" = col_cond)) +
  labs(x = expression("True effect (" * delta * ")"), y = "Power (DC)") + my_theme

pr2b <- ggplot(df_pow_abs, aes(delta, avgPower, color = Criterion)) +
  geom_line(size = 1) +
  geom_vline(xintercept = DV, linetype = "dotted", color = "gray", linewidth = 0.8) +
  scale_color_manual(values = c("DC" = col_dc, "SC" = col_sc)) +
  labs(x = expression("True effect (" * delta * ")"), y = "average power") + my_theme

grid.arrange(pr2a, pr2b, ncol = 2)


# =============================================================================
# Chunk: plot_euii
# =============================================================================
# --- EUII Data & Plots ---
df_euii_loc <- data.frame(delta = seq(0, 100, length.out = 100)) |>
  mutate(
    Cond_EUII = calc_EUII(condPower(mu_A, delta, k_SC), condPower(mu_A, 0, k_SC), N_total),
    Avg_EUII  = calc_EUII(avgPower(mu_A, sigma_A2, delta, k_SC), avgPower(mu_A, sigma_A2, 0, k_SC), N_total)
  ) |> pivot_longer(cols = c(Cond_EUII, Avg_EUII), names_to = "Metric", values_to = "Value")

df_euii_scale <- data.frame(sigma_D = sigmas) |>
  mutate(Avg_EUII = calc_EUII(avgPower(mu_A, sigma_D^2, 60, k_SC), avgPower(mu_A, sigma_D^2, 0, k_SC), N_total))

p3a <- ggplot(df_euii_loc, aes(x = delta, y = Value, color = Metric)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Avg_EUII" = col_avg, "Cond_EUII" = col_cond),
                     labels = c("Avg_EUII"  = expression(EUII[avg]),
                                "Cond_EUII" = expression(EUII[cond]))) +
  labs(x = expression("True effect (" * delta * ")"), y = "EUII") + my_theme

p3b <- ggplot(df_euii_scale, aes(x = sigma_D, y = Avg_EUII)) +
  geom_line(color = col_avg, size = 1.2) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color="black", linewidth = 0.6) +
  geom_vline(xintercept = sqrt(sigma_A2), linetype = "dotted", color = "gray", linewidth = 0.8) +
  my_x_scale +
  labs(x = expression(sigma["D,C"]), y = expression(EUII[avg])) + my_theme

grid.arrange(p3a, p3b, ncol = 2)


# =============================================================================
# Chunk: plot_euii_ratio
# =============================================================================
df_reuii_loc <- data.frame(delta = seq(50, 100, length.out = 100)) |>
  mutate(
    `conditional` = calc_EUII(condPower(mu_A, delta, k_DC), condPower(mu_A, 0, k_DC), N_total) /
                  calc_EUII(condPower(mu_A, delta, k_SC), condPower(mu_A, 0, k_SC), N_total),
    `average`     = calc_EUII(avgPower(mu_A, sigma_A2, delta, k_DC), avgPower(mu_A, sigma_A2, 0, k_DC), N_total) /
                  calc_EUII(avgPower(mu_A, sigma_A2, delta, k_SC), avgPower(mu_A, sigma_A2, 0, k_SC), N_total)
  ) |>
  pivot_longer(cols = c(`conditional`, `average`), names_to = "Metric", values_to = "Ratio")

df_euii_abs <- data.frame(delta = seq(50, 100, length.out = 100)) |>
  mutate(
    DC         = calc_EUII(avgPower(mu_A, sigma_A2, delta, k_DC), avgPower(mu_A, sigma_A2, 0, k_DC), N_total),
    SC = calc_EUII(avgPower(mu_A, sigma_A2, delta, k_SC), avgPower(mu_A, sigma_A2, 0, k_SC), N_total)
  ) |>
  pivot_longer(cols = c(DC, SC), names_to = "Criterion", values_to = "avgEUII")

pr3a <- ggplot(df_reuii_loc, aes(delta, Ratio, color = Metric)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("average" = col_avg, "conditional" = col_cond)) +
  labs(x = expression("True effect (" * delta * ")"), y = "EUII ratio (DC / SC)") + my_theme

pr3b <- ggplot(df_euii_abs, aes(delta, avgEUII, color = Criterion)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("DC" = col_dc, "SC" = col_sc)) +
  labs(x = expression("True effect (" * delta * ")"), y = "EUII") + my_theme

grid.arrange(pr3a, pr3b, ncol = 2)


# =============================================================================
# Chunk: plot_euii_asymptotic
# =============================================================================
library(gridExtra)
sigma <- 88; rT_a <- 2; DV_a <- 50; n_A_a <- 20; z_p_a <- qnorm(0.95)
cols_euii <- c("DC" = col_dc, "SC" = col_sc)

## asymptotic value with limiting threshold kappa (q for SC, DV for DC)
euii_asymptotic <- function(delta, kappa)
  exp(rT_a * ((delta - kappa)^2 + kappa^2) / (2 * sigma^2 * (1 + rT_a)^2))

## --- Left panel: asymptotic values as a function of delta ---
dd <- seq(0, 100, length.out = 200)
df_asym <- rbind(
  data.frame(delta = dd,             Criterion = "SC", EUII = euii_asymptotic(dd, 0)),
  data.frame(delta = dd[dd >= DV_a], Criterion = "DC", EUII = euii_asymptotic(dd[dd >= DV_a], DV_a))
)
p_left <- ggplot(df_asym, aes(delta, EUII, color = Criterion)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = DV_a, linetype = "dotted", color = "gray", linewidth = 0.8) +
  scale_color_manual(values = cols_euii) +
  labs(x = expression("True effect (" * delta * ")"),
       y = "Asymptotic EUII" ) + my_theme

## --- Right panel: finite-n EUII vs n_C at delta = 60, with asymptotes ---
Lstable <- function(x) pnorm(x, log.p = TRUE) - pnorm(x, lower.tail = FALSE, log.p = TRUE)  # logit Phi in log space

euii_finite <- function(n_C, criterion, delta = 60) {
  s <- sigma * sqrt(1 / (n_C + n_A_a) + 1 / (rT_a * n_C))   # see Formula (12). avgPower for aligned priors

  k  <- z_p_a * s                                          # SC threshold (q = 0)
  if (criterion == "DC") k <- max(k, DV_a)                 # DC threshold = max(k_SC, DV)
  exp((Lstable((delta - k) / s) - Lstable(-k / s)) / (n_C * (1 + rT_a)))
}

nc <- exp(seq(log(20), log(1e5), length.out = 200))

df_nc <- rbind(
  data.frame(n_C = nc, Criterion = "SC", EUII = sapply(nc, euii_finite, criterion = "SC")),
  data.frame(n_C = nc, Criterion = "DC", EUII = sapply(nc, euii_finite, criterion = "DC"))
)

asy <- data.frame(Criterion = c("SC", "DC"),
                  asy_value = c(euii_asymptotic(60, 0), euii_asymptotic(60, DV_a)))

p_right <- ggplot(df_nc, aes(n_C, EUII, color = Criterion)) +
  geom_hline(data = asy, aes(yintercept = asy_value, color = Criterion),
             linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE) +
  geom_line(linewidth = 1.1) +
  scale_x_log10() +
  scale_color_manual(values = cols_euii) +
  labs(x = expression("Control sample size " * n[C] * " (log scale)"),
       y = expression(EUII * " (" * delta * " = 60)")) + my_theme

grid.arrange(p_left, p_right, ncol = 2)

## I tried to get this asymptotic results using RBest, but the EUII converges to 1 even in the 1 Normal case prior
## (Analysis prior == Design prior ). Maybe it is numerical stabiity ?
