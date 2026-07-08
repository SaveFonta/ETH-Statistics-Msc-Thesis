# =============================================================================
# Closed-Form Running Example — Single-Component Normal Analysis Prior
# =============================================================================
# Reproduces the closed-form figures of the running example:
#   Part 1 – Conditional vs Average Type I Error (location and scale shift)
#   Part 2 – Conditional vs Average Power
#   Part 3 – Sample-size curve (average power vs n_C for several ESS n_A)
#   Part 4 – Conditional vs Average EUII
#   Part 5 – Dual criterion vs significance-only (T1E, Power, EUII)
# All metrics use the Normal closed forms; the analysis prior is the
# single-component Normal approximation of the MAP prior.
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(scales)

# ---------------------------------------------------------------------------
# Trial parameters
# ---------------------------------------------------------------------------
sigma <- 88; n_C <- 20; n_T <- 40; n_A <- 20; N_total <- n_C + n_T
v_C <- sigma^2 / n_C; v_T <- sigma^2 / n_T; sigma_A2 <- sigma^2 / n_A; mu_A <- -50
q <- 0; z_p <- qnorm(0.95)
W <- sigma_A2 / (sigma_A2 + v_C)
V_delta <- W * v_C + v_T # sigma_A2 * v_C / (sigma_A2 + v_C) + v_T
V_samp  <- W^2 * v_C + v_T # Variance of the Sampling distribution of delta_hat | delta, theta_C
k_SC <- q + z_p * sqrt(V_delta) # single-criterion (significance) threshold

# ---------------------------------------------------------------------------
# Closed-form operating characteristics for a generic threshold k
# (single criterion: k = k_SC; dual criterion: k = max(DV, k_SC))
# ---------------------------------------------------------------------------
condPower <- function(theta_C_star, delta, k)
  1 - pnorm((k - ((1 - W) * (mu_A - theta_C_star) + delta)) / sqrt(V_samp))
avgPower <- function(mu_D, sigma_D2, delta, k)
  1 - pnorm((k - ((1 - W) * (mu_A - mu_D) + delta)) / sqrt(V_samp + (1 - W)^2 * sigma_D2))
calc_EUII <- function(power, t1e, N) {
  eps <- 1e-10
  p <- pmax(pmin(power, 1 - eps), eps)
  t <- pmax(pmin(t1e, 1 - eps), eps)
  log_dor <- qlogis(p) - qlogis(t) #qlogis = logit = log(Odds)
  exp(log_dor / N)
}

# average power for the assumptions for the sample size computation
avgPower_ess <- function(delta, n_C, n_A, r_T, sigma, z_p, q = 0) {
  n_T <- r_T * n_C
  1 - pnorm(z_p - (delta - q) / (sigma * sqrt(1 / (n_C + n_A) + 1 / n_T)))
}

# ---------------------------------------------------------------------------
# Shared plotting theme / scales
# ---------------------------------------------------------------------------
my_theme <- theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())
my_x_scale <- scale_x_continuous(
  trans = "log1p",
  breaks = c(0, 10, 100, 1000, 10000),
  labels = c("0\n(Cond.)", "10", "100", "1K", "10K")
)

sigmas <- c(0, 10^seq(0, 5, length.out = 99))


# ===========================================================================
# Part 1 – Type I Error (location shift and scale shift)
# ===========================================================================
df_t1e_loc <- data.frame(mu_D = seq(-90, -10, length.out = 100))  |> 
  mutate(
    `Conditional T1E` = condPower(mu_D, 0, k_SC),
    `Average T1E`     = avgPower(mu_D, sigma_A2, 0, k_SC)
  )  |> 
  pivot_longer(cols = c(`Conditional T1E`, `Average T1E`),
               names_to = "Metric", values_to = "Value")

df_t1e_scale <- data.frame(sigma_D = sigmas)  |> 
  mutate(Avg_T1E = avgPower(mu_A, sigma_D^2, 0, k_SC))

p1a <- ggplot(df_t1e_loc, aes(x = mu_D, y = Value, color = Metric, linetype = Metric)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_vline(xintercept = mu_A, linetype = "dotted", color = "gray", linewidth = 0.8) +
  scale_color_manual(values = c("Average T1E" = "#E69F00", "Conditional T1E" = "#56B4E9")) +
  labs(x = expression("True Control Mean (" * mu["D,C"] * ")"), y = "Type I Error") +
  my_theme

p1b <- ggplot(df_t1e_scale, aes(x = sigma_D, y = Avg_T1E)) +
  geom_line(color = "#E69F00", size = 1.2) +
  geom_hline(yintercept = 0.5,  linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red",   linewidth = 0.4) +
  geom_vline(xintercept = sqrt(sigma_A2), linetype = "dotted", color = "gray", linewidth = 0.8) +
  my_x_scale +
  coord_cartesian(ylim = c(0, 0.6)) +
  labs(x = expression(sigma["D,C"]), y = "avgT1E") +
  my_theme

grid.arrange(p1a, p1b, ncol = 2)


# ===========================================================================
# Part 2 – Power (across deltas and with scale shift)
# ===========================================================================
df_pow_loc <- data.frame(delta = seq(0, 100, length.out = 100))  |> 
  mutate(
    `Conditional Power` = condPower(mu_A, delta, k_SC),
    `Average Power`     = avgPower(mu_A, sigma_A2, delta, k_SC)
  )  |> 
  pivot_longer(cols = c(`Conditional Power`, `Average Power`),
               names_to = "Metric", values_to = "Value")

df_pow_scale <- data.frame(sigma_D = sigmas)  |> 
  mutate(Avg_Power = avgPower(mu_A, sigma_D^2, 60, k_SC))

p2a <- ggplot(df_pow_loc, aes(x = delta, y = Value, color = Metric, linetype = Metric)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Average Power" = "#E69F00", "Conditional Power" = "#56B4E9")) +
  labs(x = expression("True Effect (" * delta * ")"), y = "Power") +
  my_theme

p2b <- ggplot(df_pow_scale, aes(x = sigma_D, y = Avg_Power)) +
  geom_line(color = "#E69F00", size = 1.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_vline(xintercept = sqrt(sigma_A2), linetype = "dotted", color = "gray", linewidth = 0.8) +
  my_x_scale +
  coord_cartesian(ylim = c(0.4, 1.0)) +
  labs(x = expression(sigma["D,C"]), y = "avgPower") +
  my_theme

grid.arrange(p2a, p2b, ncol = 2)


# ===========================================================================
# Part 3 – Sample size: average power vs n_C for several ESS n_A
# ===========================================================================
n_C_seq  <- seq(4, 40, by = 0.5)
n_A_vals <- c(0, 10, 20, 30)
df_ss <- do.call(rbind, lapply(n_A_vals, function(na) {
  data.frame(
    n_C   = n_C_seq,
    power = avgPower_ess(delta = 60, n_C = n_C_seq, n_A = na,
                         r_T = 2, sigma = 88, z_p = qnorm(0.95)),
    n_A   = factor(paste0("n_A = ", na))
  )
}))

# minimum n_C for 80% average power with n_A = 20 (closed form)
K_req   <- 88^2 * (qnorm(0.95) + qnorm(0.80))^2 / 60^2   # ~ 13.30
b_req   <- K_req * (2 + 1) - 2 * 20
n_C_req <- (b_req + sqrt(b_req^2 + 4 * 2 * K_req * 20)) / (2 * 2)  # ~ 11.5

legend_labels <- c("n_A = 0 (Z-test)", "n_A = 10", "n_A = 20", "n_A = 30")

p_ss <- ggplot(df_ss, aes(x = n_C, y = power, colour = n_A, linetype = n_A)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0.8, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = n_C_req, linetype = "dotted", colour = "grey40") +
  geom_vline(xintercept = 20,      linetype = "dotted", colour = "grey40") +
  annotate("text", x = n_C_req - 0.6, y = 0.15, label = paste0("n_C=", ceiling(n_C_req)),
           hjust = 1, size = 3, colour = "grey30") +
  annotate("text", x = 20 + 0.6, y = 0.15, label = "n_C=20",
           hjust = 0, size = 3, colour = "grey30") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  scale_colour_brewer(palette = "Set1", labels = legend_labels) +
  scale_linetype_discrete(labels = legend_labels) +
  labs(x = expression(n[C]), y = "Average Power", colour = NULL, linetype = NULL) +
  my_theme +
  theme(legend.position = "right")

print(p_ss)


# ===========================================================================
# Part 4 – EUII (deta varying and scale shift)
# ===========================================================================
df_euii_loc <- data.frame(delta = seq(0, 100, length.out = 100))  |> 
  mutate(
    Cond_EUII = calc_EUII(condPower(mu_A, delta, k_SC), condPower(mu_A, 0, k_SC), N_total),
    Avg_EUII  = calc_EUII(avgPower(mu_A, sigma_A2, delta, k_SC), avgPower(mu_A, sigma_A2, 0, k_SC), N_total)
  )  |> 
  pivot_longer(cols = c(Cond_EUII, Avg_EUII), names_to = "Metric", values_to = "Value")

df_euii_scale <- data.frame(sigma_D = sigmas)  |> 
  mutate(Avg_EUII = calc_EUII(avgPower(mu_A, sigma_D^2, 60, k_SC), avgPower(mu_A, sigma_D^2, 0, k_SC), N_total))

p3a <- ggplot(df_euii_loc, aes(x = delta, y = Value, color = Metric, linetype = Metric)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Avg_EUII" = "#009E73", "Cond_EUII" = "#CC79A7")) +
  labs(x = expression("True Effect (" * delta * ")"), y = "EUII") +
  my_theme

p3b <- ggplot(df_euii_scale, aes(x = sigma_D, y = Avg_EUII)) +
  geom_line(color = "#009E73", size = 1.2) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_vline(xintercept = sqrt(sigma_A2), linetype = "dotted", color = "gray", linewidth = 0.8) +
  my_x_scale +
  labs(x = expression(sigma["D,C"]), y = expression(EUII[avg])) +
  my_theme

grid.arrange(p3a, p3b, ncol = 2)


# ===========================================================================
# Part 5 – Dual criterion vs significance-only
# ===========================================================================
DV   <- 50
k_DC <- max(DV, k_SC)   # dual-criterion threshold ( in this case it is 50)

# --- 5a. Type I Error (location shift) -------------------------------------
df_rt1e_loc <- data.frame(mu_D = seq(-90, -10, length.out = 100))  |> 
  mutate(
    Conditional = condPower(mu_D, 0, k_DC) / condPower(mu_D, 0, k_SC),
    Average     = avgPower(mu_D, sigma_A2, 0, k_DC) / avgPower(mu_D, sigma_A2, 0, k_SC)
  )  |> 
  pivot_longer(cols = c(Conditional, Average), names_to = "Metric", values_to = "Ratio")

df_t1e_abs <- data.frame(mu_D = seq(-90, -10, length.out = 100))  |> 
  mutate(
    Dual         = avgPower(mu_D, sigma_A2, 0, k_DC),
    Significance = avgPower(mu_D, sigma_A2, 0, k_SC)
  )  |> 
  pivot_longer(cols = c(Dual, Significance), names_to = "Criterion", values_to = "avgT1E")

pr1a <- ggplot(df_rt1e_loc, aes(mu_D, Ratio, color = Metric, linetype = Metric)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  geom_vline(xintercept = mu_A, linetype = "dotted", color = "gray", linewidth = 0.8) +
  scale_color_manual(values = c("Average" = "#E69F00", "Conditional" = "#56B4E9")) +
  labs(x = expression("True Control Mean (" * mu["D,C"] * ")"),
       y = "Type I Error ratio (DC / SC)") +
  my_theme

pr1b <- ggplot(df_t1e_abs, aes(mu_D, avgT1E, color = Criterion)) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.4) +
  geom_vline(xintercept = mu_A, linetype = "dotted", color = "gray", linewidth = 0.8) +
  scale_color_manual(values = c("Dual" = "#D55E00", "Significance" = "#0072B2")) +
  labs(x = expression("True Control Mean (" * mu["D,C"] * ")"),
       y = "Average Type I Error") +
  my_theme

grid.arrange(pr1a, pr1b, ncol = 2)
df_t1e_abs  |> filter(Criterion == "Significance") #so only around mu_D = 88 the avgT1e is not controlled anymore!

# --- 5b. Power  ----------------------------------------------
df_pow_dc <- data.frame(delta = seq(0, 100, length.out = 100))  |>
  mutate(
    `Conditional Power` = condPower(mu_A, delta, k_DC),
    `Average Power`     = avgPower(mu_A, sigma_A2, delta, k_DC)
  )  |>
  pivot_longer(cols = c(`Conditional Power`, `Average Power`),
               names_to = "Metric", values_to = "Value")

df_pow_abs <- data.frame(delta = seq(0, 100, length.out = 100))  |>
  mutate(
    Dual         = avgPower(mu_A, sigma_A2, delta, k_DC),
    Significance = avgPower(mu_A, sigma_A2, delta, k_SC)
  )  |>
  pivot_longer(cols = c(Dual, Significance), names_to = "Criterion", values_to = "avgPower")

pr2a <- ggplot(df_pow_dc, aes(delta, Value, color = Metric, linetype = Metric)) +
  geom_line(size = 1) +
  geom_vline(xintercept = DV, linetype = "dotted", color = "gray", linewidth = 0.8) +
  scale_color_manual(values = c("Average Power" = "#E69F00", "Conditional Power" = "#56B4E9")) +
  labs(x = expression("True Effect (" * delta * ")"), y = "Power (Dual Criterion)") +
  my_theme

pr2b <- ggplot(df_pow_abs, aes(delta, avgPower, color = Criterion)) +
  geom_line(size = 1) +
  geom_vline(xintercept = DV, linetype = "dotted", color = "gray", linewidth = 0.8) +
  scale_color_manual(values = c("Dual" = "#D55E00", "Significance" = "#0072B2")) +
  labs(x = expression("True Effect (" * delta * ")"), y = "Average Power") +
  my_theme

grid.arrange(pr2a, pr2b, ncol = 2)

# --- 5c. EUII (varying delta, restricted to delta > DV = 50) ----------------
# Below DV the dual power is < 50% (powerless region), so comparing
# evidentiary efficiency there is not meaningful.
df_reuii_loc <- data.frame(delta = seq(50, 100, length.out = 100))  |>
  mutate(
    Conditional = calc_EUII(condPower(mu_A, delta, k_DC), condPower(mu_A, 0, k_DC), N_total) /
                  calc_EUII(condPower(mu_A, delta, k_SC), condPower(mu_A, 0, k_SC), N_total),
    Average     = calc_EUII(avgPower(mu_A, sigma_A2, delta, k_DC), avgPower(mu_A, sigma_A2, 0, k_DC), N_total) /
                  calc_EUII(avgPower(mu_A, sigma_A2, delta, k_SC), avgPower(mu_A, sigma_A2, 0, k_SC), N_total)
  )  |>
  pivot_longer(cols = c(Conditional, Average), names_to = "Metric", values_to = "Ratio")

df_euii_abs <- data.frame(delta = seq(50, 100, length.out = 100))  |>
  mutate(
    Dual         = calc_EUII(avgPower(mu_A, sigma_A2, delta, k_DC), avgPower(mu_A, sigma_A2, 0, k_DC), N_total),
    Significance = calc_EUII(avgPower(mu_A, sigma_A2, delta, k_SC), avgPower(mu_A, sigma_A2, 0, k_SC), N_total)
  )  |> 
  pivot_longer(cols = c(Dual, Significance), names_to = "Criterion", values_to = "avgEUII")

pr3a <- ggplot(df_reuii_loc, aes(delta, Ratio, color = Metric, linetype = Metric)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("Average" = "#009E73", "Conditional" = "#CC79A7")) +
  labs(x = expression("True Effect (" * delta * ")"), y = "EUII ratio (DC / SC)") +
  my_theme

pr3b <- ggplot(df_euii_abs, aes(delta, avgEUII, color = Criterion)) +
  geom_line(size = 1) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("Dual" = "#D55E00", "Significance" = "#0072B2")) +
  labs(x = expression("True Effect (" * delta * ")"), y = "Average EUII") +
  my_theme

grid.arrange(pr3a, pr3b, ncol = 2)
