source("scripts/prior_on_control/00.functions.R")
source("scripts/prior_on_control/00.functions.MC.R")





# REPRODUCE gsponer results for the Crohn's disease
sigma <- 88

dual.crit.95 <- decision2S(pc = c(0.95, 0.5), qc = c(0, -50), lower.tail = TRUE)
fut.40 <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)

decision_list <- list ( 
  list (success = dual.crit.95, futility = fut.40),
  list (success = dual.crit.95, futility = NULL) # not adding the futility stop in the end 
)

prior.c <- mixnorm(c(1, -49, 20), sigma = sigma, param = "mn")
prior.t <- mixnorm(c(1, 0, 8800), sigma = sigma)





res <- list(
delta.0 = oc2_seq_mc.normMix(-49 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7), #can also use p.vague instead of prior.t
delta.10 = oc2_seq_mc.normMix(-59 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7), 
delta.20 = oc2_seq_mc.normMix(-69 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7), 
delta.30 = oc2_seq_mc.normMix(-79 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7), #
delta.40 = oc2_seq_mc.normMix(-89 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7),
delta.50 = oc2_seq_mc.normMix(-99 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7),
delta.60 = oc2_seq_mc.normMix(-109 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7),
delta.70 = oc2_seq_mc.normMix(-119 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7),
delta.70 = oc2_seq_mc.normMix(-129 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7),
delta.70 = oc2_seq_mc.normMix(-139 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7)
)



final_tab <- format.results(res)  # same results as Gsponer 

euii <- compute_euii(res)





# The result coincide at the exact version
# engine.exact <- oc2S_seq.dual.normMix(
#   prior.t, prior.c,
#   c(20, 40), c(10,20), decision_list)

# res.exact <- engine.exact(-49, -49)

# delta.0 = oc2_seq_mc.normMix(-49 , -49 , prior.t, prior.c, c(20, 40), c(10,20), decision_list)
# delta.1 = oc2_seq_mc.normMix(0 ,0 , prior.t, prior.c, c(20, 40), c(10,20), decision_list) #T1E changes pointwise!! 



# ---------------------------------
# what if Gsponer used Bayesian OC ? 



run.avgoc <- function (analysis.prior = prior.c, design.prior = prior.c){
res <- list(
delta.0 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = 0, design_prior_c = design.prior, decision_list, n_sim = 1e7),
delta.10 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = 10, design_prior_c = design.prior, decision_list, n_sim = 1e7),
delta.20 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = 20, design_prior_c = design.prior, decision_list, n_sim = 1e7),
delta.30 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = 30, design_prior_c = design.prior, decision_list, n_sim = 1e7),
delta.40 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = 40, design_prior_c = design.prior,  decision_list, n_sim = 1e7),
delta.50 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = 50, design_prior_c = design.prior, decision_list, n_sim = 1e7),
delta.60 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = 60, design_prior_c = design.prior, decision_list, n_sim = 1e7),
delta.70 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = 70, design_prior_c = design.prior, decision_list, n_sim = 1e7),
delta.80 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = 70, design_prior_c = design.prior, decision_list, n_sim = 1e7),
delta.90 = avgoc2_seq_mc.normMix(prior.t, analysis.prior, c(20, 40), c(10,20), delta = 70, design_prior_c = design.prior, decision_list, n_sim = 1e7)
)

return(res)
}

res.avg <- run.avgoc()


# we can show the more robustness of avgt1e, since T1E varyies ppointwise and we are trusting a T1E that assumes the true value 
#of control and treatment is the same as the mean of the historical borrowing of control:




library(parallel)
true_c <- seq(-91, 26, 3)

pointwiseT1E <- mclapply(true_c, function(x) 
{oc2_seq_mc.normMix(x , x , prior.t, prior.c, c(20, 40), c(10,20), decision_list, n_sim = 1e7)
}, mc.cores = 5)

cat("Parallel computation worked")

saveRDS(list(res = res, res.avg = res.avg, pointwiseT1E = pointwiseT1E), file= "data/reproduce.Gsponer.rds")


cat("Results saved \n")

# inspection




#data <- readRDS("data/reproduce.Gsponer.rds")
#res <- data$res
#res.avg <- data$res.avg


#final_tab <- format.results(res)
#final_tab$per_stage
#final_tab$overall


#euii <- compute_euii(res)


#res.avg <- data$res.avg
#final_tab.avg <- format.results(res.avg)
#euii.avg <- compute_euii(res.avg)
# so the euii for the avg OC is lower, also the Power, while the T1E is increased. Does iit make sense ?
# I mean we know that assurance is always closer to 50 wrt to conditional Power. 

# 




## In the 1.3 script evaluates the classic conditional OC under different mixtures 
## In the 1.4 evaluates predictive OC under different mixture. 


#pointwiseT1E <- data$pointwiseT1E

#T1E_col <- sapply(pointwiseT1E, function(x) x$Overall[["Power"]])
#true_c <- seq(-91, 26, 3)

 
#df_T1E <- data.frame(theta_c = true_c, T1E = T1E_col)
