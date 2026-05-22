# =============================================================================
# assurance, to compute EUII I need to save the whole list of results
# =============================================================================

# this GONNA REPLACE 1.5 at some point since it will also include the classic


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

all_dec <- list (decision_list = decision_list, decision_list_nofut = decision_list_nofut,
                   decision_list_nofut.norel= decision_list_nofut.norel, decision_list_norel= decision_list_norel)



deltas <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90)

jobs <- expand.grid(
  analysis_prior = prior_names,
  design_prior   = prior_names,
  stringsAsFactors = FALSE
) 

cat("Total jobs:", nrow(jobs), "\n")

job_list <- split(jobs, seq_len(nrow(jobs)))

RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

# I want to have a super nested list: 

# $decisions[["ap.MAP_dp.MAP"]][["delta.0"]]
simulate <- function(dec_list) {
  results <- pbmclapply(job_list, function(job) {
    ap <- job$analysis_prior
    dp <- job$design_prior

    avgoc2_seq_mc.normMix(
      prior_1        = prior.t,
      prior_2        = priors[[ap]],
      n1_seq         = n1_seq,
      n2_seq         = n2_seq,
      decisions_list = dec_list,
      delta          = deltas,
      design_prior_c = priors[[dp]],
      sigma_1        = sigma,
      sigma_2        = sigma,
      n_sim          = N_SIM,
      seed           = NULL
    )
  }, mc.cores = 15
  )

  names(results) <- paste0("ap.", jobs$analysis_prior, "_dp.", jobs$design_prior)
  return(results)
}

all_results <- lapply(names(all_dec), function(dec_name) {
  simulate(all_dec[[dec_name]])
})
names(all_results) <- names(all_dec)

# Then compute EUII per job per decision list:
# euii_results <- lapply(all_results, function(dec_res) {
#   lapply(dec_res, compute_euii)
# })






saveRDS(list( all_results = all_results), file = "data/1.6")
cat("results saved!!")




decision_list <- all_results$decision_list
decision_list_nofut<- all_results$decision_list_nofut
decision_list_norel <- all_results$decision_list_norel
decision_list_nofut.norel <- all_results$decision_list_nofut.norel




# formatted 
euii <- lapply(all_results, function(dec.func) {
  lapply(dec.func, function(job) {
      compute_euii(job) #or format_results
 })
}
)



# formatted 
formatted <- lapply(all_results, function(dec.func) {
  lapply(dec.func, function(job) {
      format_results(job) 
 })
}
)