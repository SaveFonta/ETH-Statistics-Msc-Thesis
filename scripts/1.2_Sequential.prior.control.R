source("scripts/00.functions.R")
source("scripts/00.functions.MC.R")





# REPRODUCE gsponer results for the Crohn's disease
sigma <- 88

dual.crit.95 <- decision2S(pc = c(0.95, 0.5), qc = c(0, -50), lower.tail = TRUE)
fut.40 <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)

decision_list <- list ( 
  list (success = dual.crit.95, futility = fut.40),
  list (success = dual.crit.95, futility = NULL) # not adding the futility stop in the end 
)

prior.c <- mixnorm(c(1, -49, 20), sigma = sigma, param = "mn")
prior.t <- mixnorm(c(1, 0, 0.001), sigma = sigma, param = "mn")





res <- list(
delta.0 = oc2_seq_mc.normMix(-49 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7), #can also use p.vague instead of prior.t
delta.40 = oc2_seq_mc.normMix(-89 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7),
delta.50 = oc2_seq_mc.normMix(-99 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7),
delta.60 = oc2_seq_mc.normMix(-109 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7),
delta.70 = oc2_seq_mc.normMix(-119 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7)
)



final_tab <- format.results(res)  # same results as Gsponer 

euii <- compute_euii(res)





# The result coincide at the exact version
engine.exact <- oc2S_seq.dual.normMix(
  prior.t, prior.c,
  c(20, 40), c(10,20), decision_list)

res.exact <- engine.exact(-49, -49)

delta.0 = oc2_seq_mc.normMix(-49 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list)
delta.1 = oc2_seq_mc.normMix(0 ,0 , prior.t, prior.c, c(20, 40), c(10,20), decision_list) #T1E changes pointwise!! 



# ---------------------------------
# what if Gsponer used Bayesian OC ? 



run.avgoc <- function (analysis.prior = prior.c, design.prior = prior.c){
res <- list(
delta.0 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = 0, design_prior_c = design.prior, decision_list, n_sim = 1e7),
delta.40 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = -40, design_prior_c = design.prior,  decision_list, n_sim = 1e7),
delta.50 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = -50, design_prior_c = design.prior, decision_list, n_sim = 1e7),
delta.60 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = -60, design_prior_c = design.prior, decision_list, n_sim = 1e7),
delta.70 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = -70, design_prior_c = design.prior, decision_list, n_sim = 1e7)
)

return(res)
}

res.avg <- run.avgoc()




saveRDS(list(res = res, res.avg = res.avg), file= "data/reproduce.Gsponer.rds")





data <- readRDS("data/reproduce.Gsponer.rds")
res <- data$res
res.avg <- data$res.avg


final_tab <- format.results(res)
final_tab$per_stage
final_tab$overall


euii <- compute_euii(res)



final_tab.avg <- format.results(res.avg)
euii.avg <- compute_euii(res.avg)
# so the euii for the avg OC is lower, also the Power, while the T1E is increased. Does iit make sense ?
# I mean we know that assurance is always closer to 50 wrt to conditional Power. 

# 





## In the 1.3 script evaluates the classic conditional OC under different mixtures 
## In the 1.4 evaluates predictive OC under different mixture. 





















