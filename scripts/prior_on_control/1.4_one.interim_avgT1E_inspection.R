# =============================================================================
# assurance
# =============================================================================

library(RBesT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(scales)
library(parallel)
library(pbmcapply)

source("scripts/prior_on_control/00.functions.MC.R")

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




N_SIM <- 1e7
SEED  <- 123

# =============================================================================
# DEFINE PRIORS
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
# mu_robust <- -50    
# scale_robust <- 88 
# nu <- 3            # Degrees of freedom (3 is the standard for robust priors)


#set.seed(123)
#n_samples <- 500000
#t_samples <- mu_robust + scale_robust * rt(n_samples, df = nu)


#approx_robust_prior <- automixfit(t_samples, Nc = 15)

approx_robust_prior <- readRDS("data/t.prior.rds")

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

cat("Prior defined \n")






deltas <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)

jobs <- expand.grid(
  analysis_prior = prior_names,
  design_prior   = prior_names,
  Delta = deltas,
  stringsAsFactors = FALSE
) 

cat("Total jobs:", nrow(jobs), "\n")

job_list <- split(jobs, seq_len(nrow(jobs)))

results_raw <- pbmclapply(job_list, function(job) {
  
  ap <- job$analysis_prior
  dp <- job$design_prior
  d  <- job$Delta
  
  oc <- tryCatch(
    avgoc2_seq_mc.normMix(
      prior_1        = prior.t,
      prior_2        = priors[[ap]],
      n1_seq         = n1_seq,
      n2_seq         = n2_seq,
      decisions_list = decision_list,
      delta          = d,
      design_prior_c = priors[[dp]],
      sigma_1        = sigma,
      sigma_2        = sigma,
      n_sim          = N_SIM,
      seed           = SEED
    ),
    error = function(e) {
      message("ERROR | analysis=", ap, " design=", dp, " delta=", d, " | ", e$message)
      NULL
    }
  )
  
  if (is.null(oc)) return(NULL)
  
  data.frame(
    Analysis_Prior = ap,
    Design_Prior   = dp,
    Delta          = d,
    Assurance      = oc$Overall["Power"],
    MCE            = oc$Overall["MCE_Power"]
  )
  
}, mc.cores = 15)

avgT1e <- do.call(rbind, Filter(Negate(is.null), results_raw))
rownames(avgT1e) <- NULL

saveRDS(avgT1e, file = "data/avgT1e.sequential.rds")
cat("Results saved \n")

# -----------------------------------------------------------------------------
# PLOT 1: heatmap per delta — one facet per delta value
# -----------------------------------------------------------------------------
#avgT1e <- readRDS("data/avgT1e")

p_heatmap <- avgT1e |>
  ggplot(aes(x = Design_Prior, y = Analysis_Prior, fill = Assurance)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", Assurance)), size = 2.8) +
  scale_fill_gradient2(
    low      = "#2E6DA4",
    mid      = "white",
    high     = "#E04F39",
    midpoint = 0.5,
    labels   = percent_format(accuracy = 0.1),
    name     = "Assurance"
  ) +
  facet_wrap(~ Delta, ncol = 3,
             labeller = labeller(Delta = function(x) paste0("\u03b4 = ", x))) +
  labs(
    title = "Assurance: analysis prior \u00d7 design prior",
    subtitle = "Facet = true delta | Colour = assurance",
    x = "Design Prior (theta_C draws)",
    y = "Analysis Prior (posterior)"
  ) +
  theme_bw(base_size = 11) +
  theme(
    axis.text.x      = element_text(angle = 35, hjust = 1),
    panel.grid       = element_blank(),
    strip.background = element_rect(fill = "grey92")
  )

# -----------------------------------------------------------------------------
# PLOT 2: assurance curves across delta, faceted by analysis prior,
#         coloured by design prior
# -----------------------------------------------------------------------------

p_curves <- avgT1e |>
  ggplot(aes(x = Delta, y = Assurance,
             colour = Design_Prior, group = Design_Prior)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.4) +
  geom_ribbon(aes(ymin = Assurance - 1.96 * MCE,
                  ymax = Assurance + 1.96 * MCE,
                  fill = Design_Prior),
              alpha = 0.1, colour = NA) +
  scale_x_continuous(breaks = deltas) +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_colour_brewer(palette = "Dark2", name = "Design Prior") +
  scale_fill_brewer(palette = "Dark2", guide = "none") +
  facet_wrap(~ Analysis_Prior, ncol = 3) +
  labs(
    title    = "Assurance across \u03b4 by prior combination",
    subtitle = "Facet = analysis prior | Colour = design prior",
    x        = expression(delta == theta[C] - theta[T]),
    y        = "Assurance"
  ) +
  theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    panel.grid.minor = element_blank(),
    legend.position  = "bottom"
  )

print(p_heatmap)
print(p_curves)
