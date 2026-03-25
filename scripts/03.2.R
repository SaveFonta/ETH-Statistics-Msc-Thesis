# -----------------------------------------------------------------
# Example 3.2: Simon Two-Stage Design  
# -----------------------------------------------------------------
# Stage 1: n = 16. Stage 2: n = 30 (Total = 46)
#early stopping for futility, THERE IS NO EARLY STOPPING for success 


source("00_compute_EUII_2.R")
library(lattice)
library(patchwork)

pi_c = 0.4; pi_E = 0.6; alpha = 0.05; beta = 0.20


logit <- function(p) log(p / (1 - p))
expit <- function(x) { exp(x)/(1+exp(x)) } 

sigma_E <- sqrt(1/pi_E + 1/(1-pi_E))
prior_ctrl <- c(logit(0.4), 1e20) #centerd at real logit value with huge weight

patients <- cbind( c( 0, 0), 
                   c(16, 30))


succ_stage2 <- logit(23.5/46) - logit(pi_c)
crit_succ_mat <- rbind(c(NA, NA), #NO EARLY STOPPING for success
                       c( succ_stage2, 0.5))

fut_stage1 <- logit(7.5/16) - logit(0.4) 
crit_fut_mat <- rbind(c(fut_stage1, 0.5),
                      c(succ_stage2, 0.5))



design_onco <- gsbDesign(nr.stages = 2, 
                        patients = patients, 
                        sigma = c(sigma_E, sigma_E), #they use sigma_E = 2 in the paper
                        criteria.success = crit_succ_mat, 
                        criteria.futility = crit_fut_mat, 
                        prior.control = prior_ctrl, 
                        prior.treatment = c(logit(pi_c), 0)) #non informative


true_treatment_rates <- seq(0.4, 0.6, by=0.05)
deltas <- logit(true_treatment_rates) - logit(0.4)

sim_onco <- gsbSimulation(
  truth = list(control = logit(0.4), delta = deltas), # Passed as a named list for sliced grid
  type.update = "per arm",
  method = "simulation",
  grid.type = "sliced", # sliced better!!
  nr.sim = 100000,
  seed = 1
)

res_onco <- gsb(design_onco, sim_onco)
table_onco <- tab(res_onco, what = "cumulative all")
tab_n_onco <- tab(res_onco, what = "sample size")

#change the delta scale to probabilities
table_onco[, "control"] <- expit(table_onco[, "control"])
table_onco[, "treatment"] <- expit(table_onco[, "treatment"])
table_onco[, "delta"] <- table_onco[, "treatment"] - table_onco[, "control"]


tab_n_onco[, "control"] <- expit(tab_n_onco[, "control"])
tab_n_onco[, "treatment"] <- expit(tab_n_onco[, "treatment"])
tab_n_onco[, "delta"] <- tab_n_onco[, "treatment"] - tab_n_onco[, "control"]


sample_size <- plot(res_onco, what="sample size", sliced = TRUE)

OC_onco <- plot(res_onco, what="cumulative all", sliced = TRUE)






#Note, in this specific design, we know that if a trial is a success, then it must be to the end
#This mean that E(N_+) = 46, while E(N_-) can be not


EUII_onco <- compute_euii(
    table_onco, N = c(16, 46), Pr_H1 = 0.1,
                     first_order = FALSE, exact = TRUE)

comp_clean_onco <- EUII_onco$result_df

# Bind with N data
comp_clean_onco <- cbind(comp_clean_onco, EUII_onco$N[, c("E_N_plus", "E_N_minus")])

colnames(comp_clean_onco) <- c("True Rate delta ($\\pi_E - \\pi_C$)",
                               "DOR",
                               "EUII",
                               "$E[N_+]$",
                               "$E[N_-]$")

comp_tab_onco <- xtable(comp_clean_onco,
                        caption = "EUII and other operating characteristics",
                        label = "tab:comp_euii_onco",
                        digits = c(0, 2, 2, 4, 1, 1))



plot_euii_onco <- plot_euii_comparison(table_onco,N = c(16, 46), exact = TRUE, first_order = FALSE, priors = c(0.01,0.10, 0.99))






