# =============================================================================
# assurance, to compute EUII I need to save the whole list of results
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
sign.95 <- decision2S(pc = 0.95, qc = 0, lower.tail = TRUE)
relev.95 <- decision2S(pc = 0.5, qc = -50, lower.tail = TRUE)



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




# 
priors <- list(MAP = p_MAP,
               Robust_0.20 = p_rob0.2,
               Vague = p_vague,
               Skeptical = p_skep,
               Normal = Normal,
               Robust_0.5 = p_rob_0.5)
               
               
prior_names <- names(priors)
cat("Prior defined \n")



decision_list <- list(
  list(success = dual.crit.95, futility = fut.40),
  list(success = dual.crit.95, futility = NULL)
)

decision_list_nofut <- list(
  list(success = dual.crit.95, futility = NULL),
  list(success = dual.crit.95, futility = NULL)
)

decision_list_nofut.norel <- list(
  list(success = sign.95, futility = NULL),
  list(success = sign.95, futility = NULL)
)


decision_list_norel <- list(
  list(success = sign.95, futility = fut.40),
  list(success = sign.95, futility = NULL)
)


deltas <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)

jobs <- expand.grid(
  analysis_prior = prior_names,
  design_prior   = prior_names,
  Delta = deltas,
  stringsAsFactors = FALSE
) 

cat("Total jobs:", nrow(jobs), "\n")

job_list <- split(jobs, seq_len(nrow(jobs)))


RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

simulate <- function(decision_list) {
  results <- pbmclapply(job_list, function(job) {
  
  ap <- job$analysis_prior
  dp <- job$design_prior
  d  <- job$Delta
  
  oc <- avgoc2_seq_mc.normMix(
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
    )
}, mc.cores = 15)

  return(results)
}




res <- simulate(decision_list)
cat("One done \n")

res_nofut <- simulate(decision_list_nofut)
cat("2 done \n")

res_nofut.norel <- simulate(decision_list_nofut.norel)
cat("3 done \n")

res_norel <- simulate(decision_list_norel)
cat("4 done \n")


saveRDS(list( res = res, res_nofut = res_nofut, res_nofut.norel = res_nofut.norel, res_norel = res_norel), file = "data/1.6")
cat("results saved!!")


# Now with this list should feed it to compute_euii


