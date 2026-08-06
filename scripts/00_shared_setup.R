# =============================================================================
# Shared design parameters, priors, decision criteria, and plotting constants
# =============================================================================
library(RBesT)
library(ggplot2)


# ---------------------------------------------------------------------------
# Design-level constants
# ---------------------------------------------------------------------------
sigma      <- 88
mu_A       <- -50   # prior mean, reference vline
delta_MCID <- 60     # treatment effect used to plot/evaluate Power
DV         <- 50     # dual criterion decision value


# ---------------------------------------------------------------------------
# Analysis priors
# ---------------------------------------------------------------------------
p_MAP <- mixnorm(
  c(0.4848, -52.457, 21.154), c(0.4598, -47.465, 7.843), c(0.0554, -50.355, 48.164),
  sigma = sigma, param = "ms"
)
p_vague <- mixnorm(c(1, mu_A, 8800), sigma = sigma, param = "ms")

analysis_priors <- list(
  "MAP"            = p_MAP,
  "Robust (w=0.2)" = robustify(p_MAP, 0.2, mean = mu_A),
  "Robust (w=0.4)" = robustify(p_MAP, 0.4, mean = mu_A),
  "Robust (w=0.6)" = robustify(p_MAP, 0.6, mean = mu_A),
  "Robust (w=0.8)" = robustify(p_MAP, 0.8, mean = mu_A),
  "Vague"          = p_vague
)
prior_names <- names(analysis_priors)

omega_map <- c(
  "MAP"            = 0.0,
  "Robust (w=0.2)" = 0.2,
  "Robust (w=0.4)" = 0.4,
  "Robust (w=0.6)" = 0.6,
  "Robust (w=0.8)" = 0.8,
  "Vague"          = 1.0
)


# ---------------------------------------------------------------------------
# Design priors
# ---------------------------------------------------------------------------
design_priors <- list(
  "Dirac (-50)"        = mixnorm(c(1, mu_A, 1e-16), sigma = sigma, param = "ms"),
  "MAP"                = p_MAP,
  "Skeptical (-90)"    = mixnorm(c(1, -90, 17.6), sigma = sigma, param = "ms"),
  "Misspecified (-10)" = mixnorm(c(1, -10, 25),   sigma = sigma, param = "ms")
)
dprior_names <- names(design_priors)

# Every (analysis prior x design prior) combination, used by every MC sweep.
combos <- expand.grid(Analysis_Prior = prior_names,
                      Design_Prior   = dprior_names,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)


# ---------------------------------------------------------------------------
# Decision criteria
# ---------------------------------------------------------------------------
# qc is on the scale of theta_T - theta_C = -delta, so delta-thresholds flip
# sign: delta > 0 -> qc = 0, delta > DV -> qc = -DV.
#
# Single criterion:            Pr(delta > 0  | data) > 0.95
sign.crit <- decision2S(pc = 0.95, qc = 0, lower.tail = TRUE)
# Dual criterion (Gsponer et al. 2014): Pr(delta > 0 | data) > 0.95 AND
#                                        Pr(delta > DV | data) > 0.50
dual.crit <- decision2S(pc = c(0.95, 0.50), qc = c(0, -DV), lower.tail = TRUE)
# Futility (interim looks only): Pr(delta < 40 | data) >= 0.90, i.e.
#   Pr(theta_T - theta_C > -40) >= 0.90, hence the upper tail at qc = -40.
fut.crit <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)


# ---------------------------------------------------------------------------
# Shared ggplot theme and scales
# ---------------------------------------------------------------------------
my_theme <- theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())

# One shared colour ramp for the 6 analysis priors, used on every OC/EUII plot.
scale_color_prior <- function(...) scale_color_viridis_d(option = "turbo", direction = -1, ...)

# The "Vague" (no borrowing) reference line is always dashed, every other
# analysis prior solid.
scale_linetype_vague <- function() scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"), guide = "none")


# ---------------------------------------------------------------------------
# Sequential-design comparison constants (fixed vs. one-interim vs. two-interim
# scripts all compare the same four criteria, with the same colours)
# ---------------------------------------------------------------------------
design_labels <- c("SC", "SC + Futility", "DC", "DC + Futility")

crit_cols <- c("SC" = "#0072B2", "SC + Futility" = "#56B4E9",
              "DC" = "#D55E00", "DC + Futility" = "#E69F00",
              "Fixed SC" = "#0072B2", "Fixed DC" = "#D55E00")

crit_ltys <- c("SC" = "solid", "SC + Futility" = "solid",
              "DC" = "solid", "DC + Futility" = "solid",
              "Fixed SC" = "twodash", "Fixed DC" = "twodash")



ph1_ref <- 0.5
