# =============================================================================
# AVG-OC ANALYSIS: CROSSING ANALYSIS PRIORS × DESIGN PRIORS
# =============================================================================
# avgOC integrates the classic OC over a design prior on theta_c, capturing
# our uncertainty about what the true control value will be.
#
# Two prior roles:
#   Analysis prior  (prior_2): used INSIDE the decision rule to update the
#                              posterior at each interim/final analysis.
#                              Determines where the decision boundary sits.
#
#   Design prior (design_prior_c): distribution from which theta_c is drawn.
#                              Represents our belief about the true control
#                              effect BEFORE the trial runs.
#
# We cross all 4 × 4 = 16 combinations at each delta.
# delta = 0  → avg T1E  (theta_t = theta_c, no true treatment effect)
# delta > 0  → avg Power (theta_t = theta_c - delta, treatment is better)
#
# NOTE on sign convention:
#   lower.tail = TRUE means lower values are better for treatment.
#   avgoc2_seq_mc does: theta_t = theta_c_draw + delta
#   So for power we use NEGATIVE delta: delta = -50 → theta_t = theta_c - 50
# =============================================================================

library(RBesT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(scales)

source("scripts/00.functions.MC.R")

# =============================================================================
# SETUP DESIGN
# =============================================================================

sigma <- 88

dual.crit.95 <- decision2S(pc = c(0.95, 0.5), qc = c(0, -50), lower.tail = TRUE)
fut.40        <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)

decision_list <- list(
  list(success = dual.crit.95, futility = fut.40),
  list(success = dual.crit.95, futility = NULL)
)

n1_seq <- c(20, 40)
n2_seq <- c(10, 20)

# Flat treatment prior
prior.t <-  mixnorm(c(1, 0, 0.001), sigma = sigma, param = "mn")

# Deltas for avgOC: negative because lower.tail = TRUE (lower = better for trt)
deltas        <- c(0, -40, -50, -60, -70)
delta_labels  <- c("0 (T1E)", "40", "50", "60", "70")   

N_SIM <- 1e7
SEED  <- 123

# =============================================================================
# DEFINE PRIORS
# =============================================================================

p_MAP  <- mixnorm(c(0.51, -51, 19.9), c(0.44, -46.8, 7.6), c(0.05, -54.1, 51.7),
                  sigma = sigma, param = "ms")
p_rob   <- robustify(p_MAP, 0.2, mean = -50)
p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma)

prior.t <- p_vague # also here,use as prior of treatment the more vague 
p_skep  <- mixnorm(c(1, -90, 25),   sigma = sigma)

priors      <- list(MAP = p_MAP, Robust = p_rob, Vague = p_vague, Skeptical = p_skep)
prior_names <- names(priors)

# =============================================================================
# Run all combinations
# =============================================================================
# For each (analysis prior × design prior × delta) triple, call avgoc2_seq_mc. (i.e. 4x4x5 calls)
# Result is a long data frame with one row per combination.

cat("\n--- Running avg-OC: 16 prior combinations × 5 deltas =", 16 * length(deltas), "MC runs ---\n")

results_avg <- bind_rows(lapply(prior_names, function(ap_name) {      # analysis prior
  bind_rows(lapply(prior_names, function(dp_name) {                   # design prior
    cat(paste0(" Analysis prior: ", ap_name, " | Design prior: ", dp_name  ))

    bind_rows(lapply(seq_along(deltas), function(di) {
      d      <- deltas[di]
      dlabel <- delta_labels[di]

      res <- avgoc2_seq_mc.normMix(
        prior_1        = prior.t,
        prior_2        = priors[[ap_name]],   # analysis prior: help to compute the decision boundary
        n1_seq         = n1_seq,
        n2_seq         = n2_seq,
        decisions_list = decision_list,
        delta          = d,                   # negative for power (lower = better)
        design_prior_c = priors[[dp_name]],   # design prior: distribution of true theta_c
        sigma_1        = sigma,
        sigma_2        = sigma,
        n_sim          = N_SIM,
        seed           = SEED
      )

      data.frame(
        Analysis_Prior = ap_name,
        Design_Prior   = dp_name,
        Delta          = abs(d),              # store as positive
        Delta_Label    = dlabel,
        Metric         = ifelse(d == 0, "avg T1E", paste0("avg Power (delta=", abs(d), ")")),
        Power          = res$Overall[["Power"]],
        P_Fut          = res$Overall[["Prob_Fut_seq"]],
        EN_total           = res$Overall[["EN_t"]] + res$Overall[["EN_c"]],
        EN_Succ      = res$Overall[["EN_t_Succ"]] + res$Overall[["EN_c_Succ"]],
        EN_Fail      = res$Overall[["EN_t_Fail"]] + res$Overall[["EN_c_Fail"]]
      )
    }))
  }))
}))

# Factor ordering for plots
results_avg <- results_avg |>
  mutate(
    Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
    Design_Prior   = factor(Design_Prior,   levels = prior_names),
    Delta_Label    = factor(Delta_Label,    levels = delta_labels)
  )

cat("Done.\n")



saveRDS(results_avg, file = "data/avgT1E.sequential.rds")
cat ("RESULTS SAVED CORRECTLY")
