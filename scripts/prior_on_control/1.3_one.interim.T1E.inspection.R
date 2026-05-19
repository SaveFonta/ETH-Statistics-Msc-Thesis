

# =============================================================================
# OC ANALYSIS UNDER DIFFERENT CONTROL PRIORS
# =============================================================================
# For each control prior :
#   Compute T1E and Power at fixed deltas (replicating the base design)
#   Compute how T1E changes across a (theta_1, theta_2) grid
# =============================================================================

library(RBesT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)   
library(stringr)

source("scripts/prior_on_control/00.functions.MC.R")

# =============================================================================
# SETUP DESIGN
# =============================================================================

sigma <- 88

dual.crit.95 <- decision2S(pc = c(0.95, 0.5), qc = c(0, -50), lower.tail = TRUE)
fut.40        <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)

decision_list <- list(
  list(success = dual.crit.95, futility = fut.40),
  list(success = dual.crit.95, futility = NULL)  # no futility stop at final stage
)

n1_seq <- c(20, 40)
n2_seq <- c(10, 20)

# Uninformative treatment prior (effectively flat around 0), could also use p_vague, which is even more vague, maybe I should
prior.t <- mixnorm(c(1, 0, 0.001), sigma = sigma, param = "mn")

# Fixed deltas to evaluate (same as base design)
# theta_2 (control) fixed at -49; theta_1 = theta_2 + delta
theta_c_fixed <- -49
deltas        <- c(0, 40, 50, 60, 70)

N_SIM <- 1e7
SEED  <- 123

# =============================================================================
#  PRIORS
# =============================================================================

p_MAP  <- mixnorm(c(0.51, -51, 19.9), c(0.44, -46.8, 7.6), c(0.05, -54.1, 51.7),
                  sigma = sigma, param = "ms")

p_rob0.2  <- robustify(p_MAP, 0.2, mean = -50)

p_rob_0.5 <- robustify(p_MAP, 0.5, mean = -50)


p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")

prior.t <- p_vague # Ok I decided to make prior.t super vague

p_skep  <- mixnorm(c(1, -90, 25), sigma = sigma, param = "ms") #posterior distribution of a stand-alone analysis of the historical study with the “most extreme” placebo effect

Normal <- mixnorm(c(1, -49, 20), sigma = sigma, param = "mn")



###################
# In case I wante the prior of a t distribution:

# ---------------------------------------------------------
# parameters of heavy-tailed prior
# ---------------------------------------------------------
mu_robust <- -50    
scale_robust <- 88 
nu <- 3            # Degrees of freedom (3 is the standard for robust priors)


set.seed(123)
n_samples <- 500000
t_samples <- mu_robust + scale_robust * rt(n_samples, df = nu)


approx_robust_prior <- automixfit(t_samples, Nc = 15)

print(approx_robust_prior)

# Combine with Informative MAP

 p_rob_t.dist0.2 <- mixcombine(
   informative = p_MAP, 
   robust = approx_robust_prior, 
   weight = c(0.8, 0.2)
 )

p_rob_t.dist0.5 <- mixcombine(
   informative = p_MAP, 
   robust = approx_robust_prior, 
   weight = c(0.5, 0.5)
 )



# 
priors <- list(MAP = p_MAP,
               Robust_0.20 = p_rob0.2,
               Vague = p_vague,
               Skeptical = p_skep,
               Normal = Normal,
               Robust_t_0.20 = p_rob_t.dist0.2,
               Robust_0.5 = p_rob_0.5,
               Robust_t_0.5 = p_rob_t.dist0.5)
prior_names <- names(priors)





# =============================================================================
# T1E AND POWER AT FIXED theta_c = -49, and varying delta 
# So here we use the same hp by Gsponer, just see how it behaves but instead of his prior
# we use tìall the mixture we can think of 
# =============================================================================
# For each prior and each delta, run the MC and extract Overall Power.

cat("\n--- T1E and Power at fixed deltas ---\n")

results_fixed <- bind_rows(lapply(prior_names, function(pname) {
  cat(paste0(" Computing Prior: ", pname, "\n"))
  prior_c <- priors[[pname]]
  
  # neew to find how many informative components we have -> MAP, Robust and Robust_t all use the 3-component p_MAP as their base 
  n_info <- if(pname %in% c("MAP", "Robust_0.20", "Robust_t_0.20", "Robust_0.5", "Robust_t_0.5")) 3 else 1
  weight.track <- if(n_info == 3) TRUE else FALSE



  bind_rows(lapply(deltas, function(d) {
    # Note: lower.tail = TRUE in the decision, so treatment must score lower to win.
    # delta > 0 means treatment is d units LOWER (better) than control.
    res <- oc2_seq_mc.normMix(
      theta_1        = theta_c_fixed - d,   # treatment arm: lower = better
      theta_2        = theta_c_fixed,
      prior_1        = prior.t,
      prior_2        = prior_c,
      n1_seq         = n1_seq,
      n2_seq         = n2_seq,
      decisions_list = decision_list,
      sigma_1        = sigma,
      sigma_2        = sigma,
      n_sim          = N_SIM,
      seed           = SEED,
      weight.track = weight.track, 
      n_info_comps   = n_info

    )
    data.frame(
      Prior      = pname,
      Delta      = d,
      Theta_t    = theta_c_fixed - d,
      Theta_c    = theta_c_fixed,
      Power      = res$Overall["Power"],
      P_Fut      = res$Overall["Prob_Fut_seq"],
      EN_total       = res$Overall["EN_t"] + res$Overall["EN_c"],
      Weight_S1  = res$Per_Stage$Exp_info_weight[1], 
      Weight_S2  = res$Per_Stage$Exp_info_weight[2]
    )
  }))
}))

# Label delta = 0 rows as T1E
results_fixed <- results_fixed |>
  mutate(
    Metric = ifelse(Delta == 0, "T1E", paste0("Power (delta=", Delta, ")")),
    Prior  = factor(Prior, levels = prior_names)
  )

cat("Done.\n")



# =============================================================================
#  GRID theta_1 and theta_2
# here we go a bit over Gsponer, to see of the design behaves in case Gsponer was 
#  wrong by assuming theta_c = 49. 
# So we have a whole lot of different theta values we want to try 
# =============================================================================

# T1E is Power when theta_1 = theta_2, but it varies based on theta_1 value.
# So the idea is to vary both on a grid 2d --> the diagonal is T1E


theta_t_grid <- seq(-100, 0, by = 10)
theta_c_grid <- seq(-100, 0, by = 10)

grid_df <- expand.grid(Theta_t = theta_t_grid, Theta_c = theta_c_grid)

# we include only the rows where theta_t < theta_c. Since otw it would mean the treatment is harmful
grid_df<- grid_df  |> 
filter(Theta_t <= Theta_c)

results_grid <- bind_rows(lapply(prior_names, function(pname) { # loop over priors
  cat(paste0("Prior: ", pname, "\n"))
  prior_c <- priors[[pname]]
  n_info <- if(pname %in% c("MAP", "Robust_0.20", "Robust_t_0.20", "Robust_0.5", "Robust_t_0.5")) 3 else 1
  weight.track <- if(n_info == 3) TRUE else FALSE

  bind_rows(lapply(seq_len(nrow(grid_df)), function(i) { # loop each combination of theta_c and theta_t
    tt <- grid_df$Theta_t[i]
    tc <- grid_df$Theta_c[i]

    res <- oc2_seq_mc.normMix(
      theta_1        = tt,
      theta_2        = tc,
      prior_1        = prior.t,
      prior_2        = prior_c,
      n1_seq         = n1_seq,
      n2_seq         = n2_seq,
      decisions_list = decision_list,
      sigma_1        = sigma,
      sigma_2        = sigma,
      n_sim          = N_SIM,
      seed           = SEED,
      weight.track = weight.track, 
      n_info_comps   = n_info
    )
    data.frame(
      Prior   = pname,
      Theta_t = tt,
      Theta_c = tc,
      Power   = res$Overall["Power"],
      P_Fut      = res$Overall["Prob_Fut_seq"],
      EN_total       = res$Overall["EN_t"] + res$Overall["EN_c"],
      Weight_S1  = res$Per_Stage$Exp_info_weight[1], 
      Weight_S2  = res$Per_Stage$Exp_info_weight[2]
    )
  }))
}))




saveRDS(list(results_fixed = results_fixed, results_grid = results_grid), file= "data/T1E.sequential.rds")

cat ("RESULTS SAVED CORRECTLY")




t1e_seq <- readRDS("data/T1E.sequential.rds")




results_fixed <- t1e_seq$results_fixed
results_grid <- t1e_seq$results_grid






# =============================================================================
#  PLOTS
# =============================================================================
prior_colours <- c(
  MAP            = "#2f00ff",

  Normal         = "#0099cc",

  Vague          = "#8d7f7f",
  Skeptical      = "#232323",

  Robust_0.20    = "#ff0000",
  Robust_0.5     = "#ffa600",

  Robust_t_0.20  = "darkgreen",
  Robust_t_0.5   = "#11ff00"
)

results_fixed$Prior <- factor(
  results_fixed$Prior,
  levels = c("MAP", "Robust_0.20", "Robust_0.5", "Robust_t_0.20", "Robust_t_0.5", "Vague", "Skeptical", "Normal")
)


# ----------------------------------------------------------------------------
# T1E and Power at fixed deltas (line + point per prior)
# ----------------------------------------------------------------------------

results_fixed |>
  ggplot(aes(x = Delta, y = Power, colour = Prior, group = Prior)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40", linewidth = 0.5) +
  scale_colour_manual(values = prior_colours) +
  scale_x_continuous(breaks = deltas) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1)) +
  labs(
    title    = "Operating Characteristics",
    subtitle = "theta_c fixed at -49; theta_t = theta_c - delta",
    x        = "Delta",
    y        = "Rejection probability",
    colour   = "Control prior"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")  

# ----------------------------------------------------------------------------
# T1E for theta_c = -49
# ----------------------------------------------------------------------------

t1e_df <- results_fixed |> filter(Delta == 0)

t1e_df |>
  ggplot(aes(x = Prior, y = Power, fill = Prior)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0.05, linetype = "dashed", colour = "grey40", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", 100 * Power)),
            vjust = -0.4, size = 3.5) +
  scale_fill_manual(values = prior_colours) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, max(t1e_df$Power) * 1.25)) +
  labs(
    title    = "Type I Error by control prior (delta = 0, theta_t = theta_c = -49)",
    x        = "Control prior",
    y        = "T1E",
    fill     = "Control prior"
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none")















#
#The boundary crit_y1(y2) answers: "given observed control mean y2, how low must treatment mean y1 be to trigger the decision?"
#The decision dual.crit.95 requires high posterior probability that treatment is below 0 AND below -50. The critical y1 is determined by how the prior on the control arm shapes the posterior comparison.
#Skeptical prior is centered at -90 with sd=25 — it strongly believes the control is around -90, which is much better than the true -49. When data says y2 = -50, the Skeptical posterior for control is pulled toward -90 — i.e. it thinks control is better than observed. This raises the bar for treatment enormously, requiring y1 ≈ -115. T1E is near zero because treatment at truth -49 essentially never clears -115.
#Vague prior contributes almost nothing — the boundary is nearly y1 ≈ y2 - 50 (a pure data comparison shifted by the fixed threshold of -50). The bar is set by data alone, and treatment at -49 clears -100 reasonably often.
#MAP prior is well-calibrated to the true control value of -49, so it neither inflates nor deflates the boundary much relative to what the data says. Robust adds 20% flat weight which slightly loosens the prior's grip, nudging the boundary marginally easier than MAP — hence slightly higher T1E, consistent with what you see.




# 
# But this assumes that the true placebo is -49. 
# Imagine now the true placebo is -60. Which results in an alpha larger tahn 0.05 as T1E. 
# So maybe a single number would be better in this case ? Can achieve this using the avgT1E







# ----------------------------------------------------------------------------
# Pointwise T1E
# ----------------------------------------------------------------------------

null_diagonal <- results_grid |> 
  filter(Theta_t == Theta_c)

ggplot(null_diagonal, aes(x = Theta_c, y = Power, color = Prior, group = Prior)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
  scale_y_continuous(labels = scales::percent, breaks = seq(0, 0.4, by = 0.05)) +
  labs(
    title = "Pointwise Type I Error Risk Across True Control Spaces",
    subtitle = "Extracted along the Null Diagonal (Theta_T = Theta_C)",
    x = "True Control Parameter Value (Theta_C)",
    y = "Actual Type I Error Rate"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")



#2D GRID for POwe? I dont like it
power_surface <- results_grid |> 
  filter(Theta_t < Theta_c)

ggplot(power_surface, aes(x = Theta_c, y = Theta_t, fill = Power)) +
  geom_tile() +
  scale_fill_viridis_c(option = "viridis", labels = scales::percent) +
  facet_wrap(~Prior, nrow = 2) +
  labs(
    title = "Global Sequential Power Topographies",
    subtitle = "Faceted by Control Arm Prior Specification Structure",
    x = "True Control Parameter (Theta_C)",
    y = "True Treatment Parameter (Theta_T)",
    fill = "Power"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hj = 1),
    panel.spacing = unit(1, "lines")
  )










#
# SOME PLOTS of the weights? Maybe add it in the end 
#

plot_df <- results_grid |> 
select(Prior, Theta_t, Theta_c, Power, Weight_S1, Weight_S2) %>%
  mutate(
    Delta = Theta_c - Theta_t
  )


  plot_df %>%
  filter(grepl("Robust", Prior)) %>%
  ggplot(aes(
    x = Delta,
    y = Weight_S2,
    colour = Prior
  )) +
  geom_line(aes(group = interaction(Prior, Theta_c))) +
  facet_wrap(~ Theta_c) +
  theme_bw() +
  labs(
    title = "Borrowing weight vs treatment effect",
    y = "Stage 2 weight"
  )


results_grid %>%
  filter(grepl("Robust", Prior)) %>%
  filter(!Theta_c %in% c(-40, -50, -60))  |> 
  ggplot(aes(
    x = Weight_S1,
    y = Weight_S2,
    colour = Prior
  )) +
  geom_point(alpha = 0.7) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  theme_bw() +
  coord_equal()









results_grid %>%
  filter(grepl("Robust", Prior)) %>%
  ggplot(aes(
    x = Theta_c,
    y = Weight_S2,
    colour = Prior
  )) +
  geom_point() +
  geom_line() +
  theme_bw() +
  labs(
    title = "Adaptive borrowing as function of control-prior conflict",
    x = expression(theta[C]),
    y = "Stage 2 informative weight"
  )


