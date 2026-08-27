

source("scripts/00_functions.R") 
library(gsDesign)


# =============================================================================
# 1. Example
# =============================================================================


  library(RBesT)
  sigma   <- 88

  p_MAP   <- mixnorm(
    c(0.4848, -52.457, 21.154), c(0.4598, -47.465, 7.843), c(0.0554, -50.355, 48.164),
    sigma = sigma, param = "ms")

  p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")

  sign.crit <- decision2S(pc = c(0.95, 0.50), qc = c(0, -50), lower.tail = TRUE)
  fut.crit  <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)

  decisions_fut <- list(
    list(success = sign.crit, futility = fut.crit),  
    list(success = sign.crit, futility = NULL)        
  )

  n1_seq <- c(30, 60) 
  n2_seq <- c(20, 40)  

  # Conditional OC
  res_cond <- oc2_seq_mc.normMix(
    theta_1 = -50, theta_2 = -50,
    prior_1 = p_vague, prior_2 = p_MAP,
    n1_seq  = n1_seq,  n2_seq  = n2_seq,
    decisions_list = decisions_fut,
    n_sim = 1000, seed = 1
  )
  res_cond$Overall
  res_cond$Per_Stage

  # Average OC
  res_avg <- avgoc2_seq_mc.normMix(
    prior_1 = p_vague, prior_2 = p_MAP,
    n1_seq  = n1_seq,  n2_seq  = n2_seq,
    decisions_list = decisions_fut,
    delta          = c(0, 30, 60),
    design_prior_c = p_MAP,
    n_sim = 1e6, seed = 123
  )

  # avgoc2_seq_mc.normMix returns a list keyed by delta; each entry has
  # $Overall and $Per_Stage
  res_avg$delta.0

  # Sequential EUII. prior_H1 enters only here, not in the simulation
  res_euii <- compute_euii(res_avg, prior_H1 = c(0.01, 0.1, 0.5))
  res_euii[["0.01"]]





n1_seq <- c(11, 27, 40)   
n2_seq <- c(6, 14, 20)   

resT1E <- calibrate_threshold(0.05, n1_seq = n1_seq, n2_seq = n2_seq,
                                prior_cntrl = p_MAP, prior_treat = p_vague,
                                criterion = "SC")


seeds <- seq(123, 135)
res_seeds <- parallel::mclapply(seq_along(seeds), function(x){
  res <- calibrate_threshold(0.05, n1_seq = n1_seq, n2_seq = n2_seq,
                                prior_cntrl = p_MAP, prior_treat = p_vague,
                                criterion = "SC", seed = seeds[[x]])
}, mc.cores = length(seeds))

df <- data.frame(seed = seeds, 
            p = sapply(res_seeds, function(x) x$p))

df  |> summary()

ggplot(df, aes(x = "", y = p)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_jitter(width = 0.05, height = 0) +
  labs(x = NULL)
ggsave("scripts/plot.png", width = 7, height = 5)






# OBF style

builder_obf_2int <- function(p) {
  p_1 <- pnorm(qnorm(p) / sqrt(1 / 3))   
  p_2 <- pnorm(qnorm(p) / sqrt(2 / 3))   
  p_3 <- p                               
  list(
    list(success = RBesT::decision2S(pc = p_1, qc = 0, lower.tail = TRUE), futility = NULL),
    list(success = RBesT::decision2S(pc = p_2, qc = 0, lower.tail = TRUE), futility = NULL),
    list(success = RBesT::decision2S(pc = p_3, qc = 0, lower.tail = TRUE), futility = NULL)
  )
}
builder_obf_2int(0.95)

resT1E_obf <- calibrate_threshold(0.05, n1_seq = n1_seq, n2_seq = n2_seq,
                                   prior_cntrl = p_MAP, prior_treat = p_vague,
                                   criterion = "SC",
                                   decisions_builder = builder_obf_2int)




res_seeds_obf <- parallel::mclapply(seq_along(seeds), function(x){
  res <- calibrate_threshold(0.05, n1_seq = n1_seq, n2_seq = n2_seq,
                                prior_cntrl = p_MAP, prior_treat = p_vague,
                                criterion = "SC", seed = seeds[[x]],
                                   decisions_builder = builder_obf_2int)
}, mc.cores = length(seeds))
