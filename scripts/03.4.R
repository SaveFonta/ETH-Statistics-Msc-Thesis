
# =================================================================
# Example 3.4: Time-to-Event Endpoint 
# =================================================================
source("00_compute_EUII_2.R") 

margin_HR <- 1.075 
theta_NI <- log(margin_HR) 
delta_NI  <- -theta_NI #= -0.0723

#the margin is log(1.075), if the drug goes above 1.075, it mean it is harmful, so we should stop  

#we use delta = - log(HR) since the gsb package is designed for endpoints where smaller values are better, and a hazard ratio greater than 1 indicates worse outcomes for the experimental arm.

# 1800 total events / 3 stages / 2 arms = 300 events per arm per stage
events_per_stage <- 300


# As proven algebraically, standard deviation is 1 for the survival approx
sigma_survival <- c(1, 1)

#Note that success is the same for both designs
crit_succ <- rbind(c(0, 0.99995),  # Stage 1: Practically impossible to stop for success
                   c(0, 0.9995),   # Stage 2: Extremely hard to stop for success
                   c(delta_NI, 0.975))   # Stage 3: Standard 1-sided 2.5% alpha


#the futility looks P(delta < delta_NI | D) < 0.95, 

fut_bayes <- rbind(c(delta_NI, 0.95),
                   c(delta_NI, 0.90),
                   c(delta_NI, 0.80))

#original stopping criterian translated from modified Haybittle–Peto
fut_orig <- rbind(c(0, 0.99995),
                  c(0, 0.9995),
                  NA) # Forced to continue until end unless massive harm





design_orig <- gsbDesign(
  nr.stages = 3, patients = events_per_stage, sigma = sigma_survival,
  criteria.success = crit_succ, criteria.futility = fut_orig, prior.difference = "non-informative"
)

design_bayes <- gsbDesign(
  nr.stages = 3, patients = events_per_stage, sigma = sigma_survival,
  criteria.success = crit_succ, criteria.futility = fut_bayes, prior.difference = "non-informative"
)

#now set the grid of true delta values we want to evaluate where delta = - log(HR), so delta = 0 means HR = 1, delta > 0 means HR < 1 (beneficial), delta < 0 means HR > 1 (harmful)
# in particular, for a truth value of 0.1, we have delta = 0.1 = - log (HR) => HR = exp(-0.1) = 0.905, which means a 9.5% reduction in the hazard, which is a clinically meaningful benefit. 
base_grid <- seq(-0.15, 0.15, by=0.005) 

#add the null delta_NI inside
truth_grid <- sort(unique(c(base_grid, delta_NI)))

sim_profess <- gsbSimulation(
  truth = truth_grid, grid.type = "manually",
  type.update = "treatment effect", method = "numerical integration"
)



res_orig  <- gsb(design = design_orig, simulation = sim_profess)
res_bayes <- gsb(design = design_bayes, simulation = sim_profess)

#replace the delta values with the corresponding true hazard ratios for better interpretability in the output tables and plots
res_orig$OC$delta <- exp(-res_orig$OC$delta)
res_bayes$OC$delta <- exp(-res_bayes$OC$delta)


#Now compute EUII
table_orig <- tab(res_orig, what = "cumulative all" )
table_bayes <- tab(res_bayes, what = "cumulative all")

cum_events <- c(600, 1200, 1800)

EUII_orig <- compute_euii(
  table_total = table_orig, 
  N = cum_events, 
  exact = TRUE, null_value = 1.075
)

EUII_bayes <- compute_euii(
  table_total = table_bayes, 
  N = cum_events, 
  exact = TRUE, null_value = 1.075
)

#what if we used the patients, not the events?
cum_pat <- c(600, 1200, 1800) / 0.09
EUII_orig_pat <- compute_euii(
  table_total = table_orig, 
  N = cum_pat, 
  exact = TRUE, null_value = 1.075
)

EUII_bayes_pat <- compute_euii(
  table_total = table_bayes, 
  N = cum_pat, 
  exact = TRUE, null_value = 1.075
)



orig_df  <- cbind(EUII_orig$result_df,  EUII_orig$N[, c("E_N_plus", "E_N_minus")])
bayes_df <- cbind(EUII_bayes$result_df, EUII_bayes$N[, c("E_N_plus", "E_N_minus")])

# Merge them together by Delta
comparison_df <- merge(orig_df, bayes_df, by = "Delta", suffixes = c("_Orig", "_Bayes"))



# Only shows some rows for table  (which ones????????????????????)
target_HRs <- c(0.86, 0.90, 0.95, 1.00, 1.05, 1.075, 1.10, 1.15)

# Find the rows in our granular grid that closest match these target HRs
keep_idx <- sapply(target_HRs, function(hr) which.min(abs(comparison_df$Delta - hr)))
table_subset <- comparison_df[unique(keep_idx), ]

table_clean <- table_subset[, c("Delta", 
                                "EUII_Exact_Orig", "EUII_Exact_Bayes", 
                                "E_N_plus_Orig", "E_N_plus_Bayes", 
                                "E_N_minus_Orig", "E_N_minus_Bayes")]

colnames(table_clean) <- c("True HR", 
                           "EUII (Orig)", "EUII (Bayes)", 
                           "$E[Events_+]$ (Orig)", "$E[Events_+]$ (Bayes)", 
                           "$E[Events_-]$ (Orig)", "$E[Events_-]$ (Bayes)")

# Print nicely to the console to verify
print("Filtered Comparison Table:")
print(table_clean, row.names = FALSE)



#EUII smaller than 1 when HR smaller than 1. INspection: 
EUII_bayes$result_df #WHEN HR bigger than 1, DOR is smaller 1, and this mean also EUII must be smaller 1. 

EUII_bayes$LR_df #when hr BIGGER ONE, LR_+ smaller than one, this mean Power < T1E. --> as established by Lehmann (1959), a test where the power to reject the null is smaller than the TIE rate is formally called a "biased test".
#that is why the DOR is smaller than 1 

# WHY?
# think at DOR = O(H1|+) / O (H1|-). If it is smaller than 1 menan that a significant result actually leaves you with lower odds of H1  being true thana  non significant result.
# this happens since H1: elta > 0-
# The consequence here is that every additional patient you enroll actively degrades the evidentiary value of the trial. You are burning resources to collect data that actively confuses the test's ability to answer the primary question
#adding the futilitu we can save a


#The EUII < 1.0 reflects that each additional unit of data is actively reducing our posterior probability that the treatment is superior to the control.
# infact it is smaller than one when the true HR is bigger than 1, means the drug is harmful












# some thing to remember, from the paper: First recurrence of stroke occurred in 916 of 10,181
#patients (9%) in the test group and 898 of 10,151 patients (8.8%)
# in the control group.

#this mean that the event rate is around 0.9
# so to see 1800 strokes, you need to enroll 1800 / 0.09 = 20k patients






library(lattice)

#  Extract the clean data
dep_orig <- subset(res_orig$OC, type %in% c("cumulative futility", "cumulative success"))
dep_orig$stage <- paste("Original Design:", dep_orig$stage)

dep_bayes <- subset(res_bayes$OC, type %in% c("cumulative futility", "cumulative success"))
dep_bayes$stage <- paste("Bayesian Design:", dep_bayes$stage)

dep_combined <- rbind(dep_orig, dep_bayes)

# Plot 1
p1_hazard <- xyplot(value ~ delta| stage, 
             groups = type, 
             data = dep_combined, 
             col = c("red", "blue"), 
             lty = c(2, 1), 
             layout = c(3, 2), 
             as.table = TRUE, 
             scales = list(alternating = 1),
             ylab = "Cumulative Probability of Stopping",
             ylim = c(-0.05, 1.05), 
             xlab = "True Hazard Ratio\n",
             panel = function(...){
                 panel.grid(h = -1, v = -1, col = "grey", lwd = 1, lty = 2)
                 panel.xyplot(..., type = "l", lwd = 3)
             },
             key = list(columns = 2, space = "bottom", 
                        lines = list(lty = c(2,1), col=c("red", "blue"), lwd = 2.5),
                        text = list(c("Futility", "Success"), col = 1), border = FALSE))

print(p1_hazard)

# Plot 2
xx_orig <- subset(res_orig$OC, type == "sample size")
xx_bayes <- subset(res_bayes$OC, type == "sample size")

xx_combined <- rbind(cbind(xx_orig, design = "Original Design"),
                     cbind(xx_bayes, design = "Bayesian Futility Design"))

xx_final <- subset(xx_combined, stage == "stage 3")
xx_final$value_patients <- xx_final$value / 0.09 # Convert events to patients

p2_hazard <- xyplot(value_patients ~ delta| design, 
             data = xx_final,
             scales = list(alternating = 1),
             panel = function(...){
                 panel.grid(h = -1, v = -1, col = "grey", lwd = 1, lty = 2)
                 panel.xyplot(..., type = "l", lwd = 3, col = "black")
             }, 
             ylim = c(1300 / 0.09, 1850 / 0.09), 
             ylab = "Expected Number of Enrolled Patients", 
             xlab = "True Hazard Ratio")

print(p2_hazard)

 


