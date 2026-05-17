

# =============================================================================
# OC ANALYSIS UNDER DIFFERENT CONTROL PRIORS
# =============================================================================
# For each control prior (p_MAP, p_rob, p_vague, p_skep):
#   1. Compute T1E and Power at fixed deltas (replicating the base design)
#   2. Compute how T1E changes across a (theta_1, theta_2) grid
#   3. Plot all results
# =============================================================================

library(RBesT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)   
library(stringr)

source("scripts/00.functions.MC.R")

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

p_rob  <- robustify(p_MAP, 0.2, mean = -50)

p_rob_0.5 <- robustify(p_MAP, 0.5, mean = -50)


p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma)

prior.t <- p_vague # Ok I decided to make prior.t super vague

p_skep  <- mixnorm(c(1, -90, 25), sigma = sigma) #posterior distribution of a stand-alone analysis of the historical study with the “most extreme” placebo effect


###################
# In case I wante the prior of a t distribution:
library(RBesT)

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

 p_rob_t.dist <- mixcombine(
   informative = p_MAP, 
   robust = approx_robust_prior, 
   weight = c(0.8, 0.2)
 )

p_rob_t.dist0.5 <- mixcombine(
   informative = p_MAP, 
   robust = approx_robust_prior, 
   weight = c(0.5, 0.5)
 )




priors <- list(MAP = p_MAP, Robust = p_rob, Vague = p_vague, Skeptical = p_skep, Robust_t = p_rob_t.dist, Robust_0.5 = p_rob_0.5, Robust_t_0.5 = p_rob_t.dist0.5)
prior_names <- names(priors)





# =============================================================================
# T1E AND POWER AT FIXED theta_c = -49, and varying delta 
# =============================================================================
# For each prior and each delta, run the MC and extract Overall Power.
# delta = 0  -> T1E (treatment = control, no true effect)
# delta > 0  → Power (treatment better by delta)

cat("\n--- T1E and Power at fixed deltas ---\n")

results_fixed <- bind_rows(lapply(prior_names, function(pname) {
  cat(paste0(" Computing Prior: ", pname, "\n"))
  prior_c <- priors[[pname]]
  
  # neew to find how many informative components we have -> MAP, Robust and Robust_t all use the 3-component p_MAP as their base 
  n_info <- if(pname %in% c("MAP", "Robust", "Robust_t", "Robust_0.5", "Robust_t_0.5")) 3 else 1
  weight.track <- if(pname %in% c("MAP", "Robust", "Robust_t", "Robust_0.5", "Robust_t_0.5")) TRUE else FALSE


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
print(results_fixed |> select(Prior, Delta, Power, P_Fut, EN_total))

# =============================================================================
#  GRID theta_1 and theta_2
# =============================================================================

# T1E is Power when theta_1 = theta_2, but it varies based on theta_1 value.
# So the idea is to vary both on a grid 2d --> the diagonal is T1E

cat("\n--- Rejection probability over (theta_1, theta_2) grid ---\n")

theta_t_grid <- seq(-120, -30, by = 10)
theta_c_grid <- seq(-120, -30, by = 10)

grid_df <- expand.grid(Theta_t = theta_t_grid, Theta_c = theta_c_grid)

results_grid <- bind_rows(lapply(prior_names, function(pname) {
  cat(paste0("Prior: ", pname, "\n"))
  prior_c <- priors[[pname]]

  bind_rows(lapply(seq_len(nrow(grid_df)), function(i) {
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
      seed           = SEED
    )
    data.frame(
      Prior   = pname,
      Theta_t = tt,
      Theta_c = tc,
      Power   = res$Overall["Power"],
      P_Fut      = res$Overall["Prob_Fut_seq"],
      EN_total       = res$Overall["EN_t"] + res$Overall["EN_c"]
    )
  }))
}))

cat("Done.\n")



saveRDS(list(results_fixed = results_fixed, results_grid = results_grid), file= "data/T1E.sequential.rds")

cat ("RESULTS SAVED CORRECTLY")




res <- readRDS("data/T1E.sequential.rds")




results_fixed <- res$results_fixed
results_grid <- res$results_grid





















