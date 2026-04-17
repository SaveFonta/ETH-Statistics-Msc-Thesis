library("RBesT")


# FOR STARTERS: let's implement the Crohn's Disease example!

# OF course we assume the endpoint is Normally distributed. 
# RECALL --> negative values means an improvement

dat <- crohn
crohn_sigma <- 88 #standard deviation from previous studies 
dat$y.se <- crohn_sigma / sqrt(dat$n)


# OBTAIN THE MAP PRIOR
# they use an interesting tau prior as half normal. Read documentation if neeed, maybe for now skip it

library(RBesT)
set.seed(689654)
map_mcmc <- gMAP(cbind(y, y.se) ~ 1 | study,
  weights = n, data = dat,
  family = gaussian,
  beta.prior = cbind(0, crohn_sigma),
  tau.dist = "HalfNormal", tau.prior = cbind(0, crohn_sigma / 2)
)
print(map_mcmc)



# Approximation of MAP Prior using a Mixture Distribution
# Now we need to convert the MAP prior to a parametric representation 
# we basically fit a parametric mixture representation using expectation-maximization (EM). The number of mixture components is chosen automatically using AIC. 
# One can also specify the number of components for the mixture via mixfit function and compare with the automixfit outcome.


map <- automixfit(map_mcmc)
print(map)
plot(map)$mix


# Now lets fit only one Normal
one <- mixfit(map_mcmc, Nc = 1) 
round(ess(one)) #almost 20 patients, not bad, similar to the 20 from Gsponer .

# Compute effective sample size (ESS)
round(ess(map)) ## default elir method

# there are many methods to compute ESS but I dont really care abt them 





# ROBUSTIFY
# to protects against T1E inflation in presence of prior-data conflict. For the 
# normal case we strongly recommend explicitly choosing the mean of the robust 
# component. 

# We use -50 consistent with the mean of the MAP prior. Furthermore, 20% probability 
# is used for the additional robust (unit-information) mixture component. T

# this 20 is the the possibility of a non-exchangable control group to be enrolled per inclusion/exclusion criteria in the current trial as compared to the historical control group population.
# Note that robustification decreases the ESS.

## add a 20% non-informative mixture component
map_robust <- robustify(map, weight = 0.2, mean = -50)




# setup the DECISION FUNCTION

# Note that the way we specify the the criteria success is the same as GsbDesign 

# Here we have: DOUBLE CRITERION 
# P (delta < 0| data ) > 0.95
# P (delta < -50| data) > 0.5

# where delta = active - prior

# TWO NOTES: 1) negative values are beneficial 
# 2) in the https://opensource.nibr.com/RBesT/articles/introduction_normal.html, they dont contiditon the criterion on data
# i think because they forgot?  


decision <- decision2S(pc = c(0.95, 0.5), qc = c(0, -50), lower.tail = TRUE) #two sample decision funtion 
print(decision)

# NOTE THAT WE HAVE lower.tail !! Foundamental, we can have futility criterion











# LETS GO WITH THE DESIGN:

## set up prior for active group (non informative)
weak_prior <- mixnorm(c(1, -50, 1), sigma = crohn_sigma, param = "mn")

# so now we can evaluate different priors on control:
# weak_prior
# map 
# map_robust

n_act <- 40
n_pbo <- 20





## four designs

# oc2S = Operating characteristincs for 2 sample design
# creates a FUNCTION that given theta1 and theta2 as input gives the operating char

## "b" means a balanced design, 1:1
## "ub" means 40 in active and 20 in placebo

design_noprior_b <- oc2S(weak_prior, weak_prior, n_act, n_act, decision,
  sigma1 = crohn_sigma, sigma2 = crohn_sigma
)
design_noprior_ub <- oc2S(weak_prior, weak_prior, n_act, n_pbo, decision,
  sigma1 = crohn_sigma, sigma2 = crohn_sigma
)
design_nonrob_ub <- oc2S(weak_prior, map, n_act, n_pbo, decision,
  sigma1 = crohn_sigma, sigma2 = crohn_sigma
)
design_rob_ub <- oc2S(weak_prior, map_robust, n_act, n_pbo, decision,
  sigma1 = crohn_sigma, sigma2 = crohn_sigma
)

# NOTE THAT: 
# 1 = active
# 2 = placebo

# there is also the function pos2S, but instead assume parameters are distributed following a distribution 
# not as fixed values

# basically pos2S integrates also over the uncertainty of the parameters 



# this creates a function that accepts mixture distributios that we will then feed
pos_fun  <- pos2S(prior1 = weak_prior, prior2 = map_robust, n1 = n_act, n2 = n_act, decision =decision,
  sigma1 = crohn_sigma, sigma2 = crohn_sigma
)

# let's integrate! 
# assume true is: 
# active is weak 
# placebo is map robust

res <- pos_fun(map_robust, map_robust)
# this is the avergae T1E? --> NO, this is not the average T1E like the paper, 
# since here the two distributions are assumed to be independent, this means that 
# we are not evaluating the case that theta_c = theta_t with theta sampled from p(theta_c)





# T1E 
# the range for true values
cfb_truth <- seq(-120, -40, by = 1)

typeI1 <- design_noprior_b(cfb_truth, cfb_truth)
typeI2 <- design_noprior_ub(cfb_truth, cfb_truth)
typeI3 <- design_nonrob_ub(cfb_truth, cfb_truth)
typeI4 <- design_rob_ub(cfb_truth, cfb_truth)

ocI <- rbind(
  data.frame(
    cfb_truth = cfb_truth, typeI = typeI1,
    design = "40:40 with non-informative priors"
  ),
  data.frame(
    cfb_truth = cfb_truth, typeI = typeI2,
    design = "40:20 with non-informative priors"
  ),
  data.frame(
    cfb_truth = cfb_truth, typeI = typeI3,
    design = "40:20 with non-robust prior for placebo"
  ),
  data.frame(
    cfb_truth = cfb_truth, typeI = typeI4,
    design = "40:20 with robust prior for placebo"
  )
)
ggplot(ocI, aes(cfb_truth, typeI, colour = design)) +
  geom_line() +
  ggtitle("Type I Error") +
  xlab(expression(paste("True value of change from baseline ", mu[act] == mu[pbo]))) +
  ylab("Type I error") +
  coord_cartesian(ylim = c(0, 0.2)) +
  theme(legend.justification = c(1, 1), legend.position = c(0.95, 0.85))








# POWER 

delta <- seq(-80, 0, by = 1)
m <- summary(map)["mean"]
cfb_truth1 <- m + delta # active for 1
cfb_truth2 <- m + 0 * delta # pbo for 2

power1 <- design_noprior_b(cfb_truth1, cfb_truth2)
power2 <- design_noprior_ub(cfb_truth1, cfb_truth2)
power3 <- design_nonrob_ub(cfb_truth1, cfb_truth2)
power4 <- design_rob_ub(cfb_truth1, cfb_truth2)

ocP <- rbind(
  data.frame(
    cfb_truth1 = cfb_truth1, cfb_truth2 = cfb_truth2,
    delta = delta, power = power1,
    design = "40:40 with non-informative priors"
  ),
  data.frame(
    cfb_truth1 = cfb_truth1, cfb_truth2 = cfb_truth2,
    delta = delta, power = power2,
    design = "40:20 with non-informative priors"
  ),
  data.frame(
    cfb_truth1 = cfb_truth1, cfb_truth2 = cfb_truth2,
    delta = delta, power = power3,
    design = "40:20 with non-robust prior for placebo"
  ),
  data.frame(
    cfb_truth1 = cfb_truth1, cfb_truth2 = cfb_truth2,
    delta = delta, power = power4,
    design = "40:20 with robust prior for placebo"
  )
)

ggplot(ocP, aes(delta, power, colour = design)) +
  geom_line() +
  ggtitle("Power") +
  xlab("True value of difference (act - pbo)") +
  ylab("Power") +
  scale_y_continuous(breaks = c(seq(0, 1, 0.2), 0.9)) +
  scale_x_continuous(breaks = c(seq(-80, 0, 20), -70)) +
  geom_hline(yintercept = 0.9, linetype = 2) +
  geom_vline(xintercept = -70, linetype = 2) +
  theme(legend.justification = c(1, 1), legend.position = c(0.95, 0.85))





























# REAL WORLD

# Now imagine in a real world scenarion, we observe the real values, 

## one can either use summary data or individual data. See ?postmix.
y.act <- -29.2
y.act.se <- 14.0
n.act <- 39

y.pbo <- -63.1
y.pbo.se <- 13.9
n.pbo <- 20

## first obtain posterior distributions
post_act <- postmix(weak_prior, m = y.act, se = y.act.se)
post_pbo <- postmix(map_robust, m = y.pbo, se = y.pbo.se)

## then calculate probability for the dual criteria
## and compare to the predefined threshold values
p1 <- pmixdiff(post_act, post_pbo, 0)
print(p1)


# the funtion postmix can really save my ass, we can develop interim analysis 
# from here --> use the posterior at interim 1 as the prior and use interim 2! 

# Now we can evaluate the criteria and see that the real deisgn is rejected
decision(post_act, post_pbo)


























# SEQUENTIAL DESIGN
ia <- data.frame(
  n = c(12, 14),
  median_count = c(20.5, 21),
  mean_count = c(23.3, 27),
  mean_log = c(2.96, 3.03),
  sd_log = c(0.67, 0.774),
  row.names = c("active", "placebo")
)  |> 
  transform(se_log = round(sd_log / sqrt(n), 3))
sd_log_pooled <- with(ia, sqrt(sum(sd_log^2 * (n - 1)) / (sum(n) - 2)))


n <- 21 # planned total n per arm
rules <- decision2S(c(0.9, 0.5), c(0, -0.357), lower.tail = TRUE)
print(rules)




# first, define non informative prior 
priorP <- priorT <- mixnorm(c(1, log(20), 1), sigma = 0.47, param = "mn")


# now let's use the data we have from the interim analysis (ia) and create posterior from there. 
postT_interim <- postmix(priorT, m = ia["active", "mean_log"], se = ia["active", "se_log"])
postP_interim <- postmix(priorP, m = ia["placebo", "mean_log"], se = ia["placebo", "se_log"])


# NOW WE HAVE A FUNDAMENTAL FUNCTION. pos2S
# 

pos_final <- pos2S(
  postT_interim, 
  postP_interim, 
  n - ia["active", "n"], #patients left to recruit
  n - ia["placebo", "n"],  # patients left to recruit
  rules, 
  sigma1 = sd_log_pooled, 
  sigma2 = sd_log_pooled
)


# If the patients we haven't seen yet are like the ones we just saw 
# (represented by the interim posteriors), what is the probability we will 
# meet our rules at the final analysis?
pos_final(postT_interim, postP_interim)

#VEERY low





# We could even compute conditional power assuming some true value
ia_oc <- oc2S(
  postT_interim,
  postP_interim,
  n - ia["active", "n"],
  n - ia["placebo", "n"],
  rules,
  sigma1 = sd_log_pooled,
  sigma2 = sd_log_pooled
)

delta <- seq(0, 0.9, 0.01) # pct diff from pbo
pbomean <- ia["placebo", "mean_log"]
y1 <- log(exp(pbomean) * (1 - delta)) # active
y2 <- log(exp(pbomean) * (1 - 0 * delta)) # placebo

out <-
  data.frame(
    diff_pct = delta,
    diff = round(y1 - y2, 2),
    y_act = y1,
    y_pbo = y2,
    cp = ia_oc(y1, y2)
  )

ggplot(data = out, aes(x = diff_pct, y = cp)) +
  geom_line() +
  scale_x_continuous(labels = scales::percent) +
  labs(
    y = "Conditional power",
    x = "True percentage difference from placebo in lesion count",
    title = "Conditional power at interim for success at final analysis"
  )
































  # LLets play
###############################
try <- mixnorm(rob = c(0.2, 0, 2), inf = c(0.8, 2, 2), sigma = 5)

# you can choose param = "ms" (default) or "mn"
# with "ms", you specify the two components as c(weight, mean , sd)

# with param = "mn"
# you specify as (weight, mean, n)

# note that it doesnt matter much, since you need to specify the sigma
# ofc se = sigma / sqrt(n) where infact sigma = standard deviation
try
summary(try)  
plot(try)  
ess(try)

# usually this mixture come from historical data, so we want to robustiify it to be 
# covered against T1E iincrease
try_rob <- robustify(try, 0.2, -50, sigma = 5)
#methods(class = "normMix")
plot(try_rob)


# once you have this amaxing mixture that will represent your prior, we need 
# to mix it with your data!!!

# postmix is the solution 
# imagine trial has 20 patients, and the observed mean is -45
post_mix_summary <- postmix(try_rob, m = -45, n = 20)

# Note you can even pass raw data and it will compute m and n based on sigma 

# this post_mix_summary is just another mixture of normals! That can be used for an interim analysis,
# we can just add new data with postmix and then bum 






# I TRY TO CODE a interim analysis

sigma_trial <- 5
prior_ctrl <- mixnorm(inf = c(0.8, -50, 2), rob = c(0.2, 0, 2), sigma = sigma_trial)
prior_trt <- mixnorm(vague = c(1, 0, 100), sigma = sigma_trial)


# Interim Futility: Stop if we are 90% sure the difference (Trt - Ctrl) is LESS than 2
futility_interim <- decision2S(pc = 0.90, qc = 2, lower.tail = TRUE)

# Interim Success: Stop early ONLY if we are 99% sure the difference is GREATER than 0
success_interim <- decision2S(pc = 0.99, qc = 0, lower.tail = FALSE)

success_final <- decision2S(
  pc = c(0.95, 0.50), 
  qc = c(0, 5), 
  lower.tail = FALSE
)




sequential_simulation <- function(n_sim = 10000, true_mean = -50, sigma_trial = 5) {
  early_success.count <- 0
  early_fut.count <- 0
  late_success.count <- 0

  for (i in seq_len(n_sim)) {
    data_ctrl <- rnorm(50, mean = true_mean, sd = sigma_trial)
    data_trt <- rnorm(50, mean = true_mean, sd = sigma_trial)

    post_ctrl_int <- postmix(prior_ctrl, data = data_ctrl[1:20])
    post_trt_int <- postmix(prior_trt, data = data_trt[1:20])

    if (success_interim(post_trt_int, post_ctrl_int) == 1) {
      early_success.count <- early_success.count + 1
      next
    }

    if (futility_interim(post_trt_int, post_ctrl_int) == 1) {
      early_fut.count <- early_fut.count + 1
      next
    }

    post_ctrl_fin <- postmix(prior_ctrl, data = data_ctrl)
    post_trt_fin <- postmix(prior_trt, data = data_trt)

    if (success_final(post_trt_fin, post_ctrl_fin) == 1) {
      late_success.count <- late_success.count + 1
    }
  }

  list(
    early_success.count = early_success.count,
    early_fut.count = early_fut.count,
    late_success.count = late_success.count,
    type_1_error = (early_success.count + late_success.count) / n_sim
  )
}

sequential_results <- sequential_simulation(1000)
print(sequential_results)



# need to compute power 

# note that witht this we could even put futility and success in the same criteria. 
# BUT this si not helpful for us cause in that way we couldnt be able to 
# distinguish if we stopped for futility or for success .
intermediateCrit2 <- decision2S(
  c(0.95, 1 - 0.5, 1 - 0.90),
  c(0, 50, 40),
  c(FALSE, TRUE, TRUE)
)




























# EXPERIMENT
sequential_simulation <- function(n_sim = 10000, design_prior_cntrl = NULL, true_delta = 0,  fixed_mean_ctrl = -50, sigma_trial = 5) {

  early_success.count <- 0
  early_fut.count <- 0
  late_success.count <- 0


  # Generate the parameters
  if (!is.null(design_prior_ctrl)) {
    # If a design prior is provided, sample true control means from it (Average Type I error)
    true_means_ctrl <- rmix(design_prior_ctrl, n_sim)
  } else {
    # Fallback to a fixed point null (Classis Type I error)
    true_means_ctrl <- rep(fixed_mean_ctrl, n_sim)
  }

  # Deterministically tied
  true_means_trt <- true_means_ctrl + true_delta

  for (i in seq_len(n_sim)) {
    data_ctrl <- rnorm(50, mean = true_mean, sd = sigma_trial)
    data_trt <- rnorm(50, mean = true_mean, sd = sigma_trial)

    post_ctrl_int <- postmix(prior_ctrl, data = data_ctrl[1:20])
    post_trt_int <- postmix(prior_trt, data = data_trt[1:20])

    if (success_interim(post_trt_int, post_ctrl_int) == 1) {
      early_success.count <- early_success.count + 1
      next
    }

    if (futility_interim(post_trt_int, post_ctrl_int) == 1) {
      early_fut.count <- early_fut.count + 1
      next
    }

    post_ctrl_fin <- postmix(prior_ctrl, data = data_ctrl)
    post_trt_fin <- postmix(prior_trt, data = data_trt)

    if (success_final(post_trt_fin, post_ctrl_fin) == 1) {
      late_success.count <- late_success.count + 1
    }
  }

  list(
    early_success.count = early_success.count,
    early_fut.count = early_fut.count,
    late_success.count = late_success.count,
    type_1_error = (early_success.count + late_success.count) / n_sim
  )
}

sequential_results <- sequential_simulation(1000)
print(sequential_results)


























