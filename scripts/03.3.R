
# -----------------------------------------------------------------
# Example 3.3: Multiple Sclerosis
# -----------------------------------------------------------------
source("00_compute_EUII_2.R")
library(xtable)

#Some params to vary
stages <- 2
patients_uninfo <- 25 
patients_info <- 20
prior_patients <- 10 #worth 10 each group 


# -------- NON informative prior
mu_c <- 14.8       # True mean lesion count for control
rate_c <- log(mu_c) #log-rate for control
kappa <- 2        
redscen <- seq(0, 0.9, 0.05) #reduction scenarios from 0% to 90%

# Note that log (r) = N (log(mu), 1/(mu*n) + k/n)  where r is the estimate of mu

sigma_c <- sqrt(kappa + (1 / mu_c))  #control is always the same 
OCout_noninf=NULL

#loop since sigma_E is tied to mu_E --> we need to recompute it for each mu_E scenario! 
for(i in seq(along=redscen)){
  red <- redscen[i]
  
  # Calculate true mean for experimental based on this reduction
  mu_E <- mu_c * (1 - red) 
  rate_E <- log(mu_E)
  
  # Delta is the log relative risk (positive value = benefit)
  delta <- rate_c - rate_E
  
  # CORRECTED VARIANCE: Use the raw mean (mu), in the paper they were using r instead of mu...
  sigma_E <- sqrt(kappa + (1 / mu_E))
  
  # Create the design for this specific variance scenario
  design1 <- gsbDesign(
    nr.stages = stages,
    patients = patients_uninfo, #100 patients are used w a 1:1 ratio 
    sigma = c(sigma_c, sigma_E), # from documentation: c(sigma.control, sigma.treatment). --> in the paper they used in the wrong order
    criteria.success = c(0, 0.9, log(2.5), 0.5), 
    criteria.futility = c(log(1.43), 0.5, log(2.5), 0.9),
    prior.difference = "non-informative" 
  )
  
  simulation1 <- gsbSimulation(
    truth = c(-1, delta, 2), # evaluate this delta range, truth = c(min, max, n), which is the same as evaluating only -1 and delta
    type.update = "treatment effect", #with type.update "trt effect" we can use numerical integratin
    method = "numerical integration"
  )
  
  x1 <- gsb(design = design1, simulation = simulation1)
  
  sub <- x1$OC[,"delta"] != -1 #we dont care about the -1
  xout <- x1$OC[sub,]
  xout$delta <- red * 100 # Store as % reduction for plotting
  OCout_noninf <- rbind(OCout_noninf, xout)
}





# -------- Informative prior


# A 75% reduction means the prior mean is 25% of the placebo mean
mu_prior <- mu_c * (1 - 0.75) 
sigma_prior_E <- sqrt(kappa + (1 / mu_prior))

V_target = (sigma_c^2 + sigma_prior_E^2) / prior_patients

# Calculate the log relative risk of the prior
delta_prior <- rate_c - log(mu_prior) #log(100/25) = log(4) 

OCout_inf <- NULL

for(i in seq(along=redscen)){
  red <- redscen[i]
  
  # Calculate true mean for experimental arm based on reduction
  mu_E <- mu_c * (1 - red) 
  rate_E <- log(mu_E)
  delta <- rate_c - rate_E
  
  sigma_E <- sqrt(kappa + (1 / mu_E))
  
  n_adj <- (sigma_c^2 + sigma_E^2) / V_target

  # Design incorporating the prior: c(Effect, n_control, n_treatment)
  design_inf <- gsbDesign(
    nr.stages = 2,
    patients = patients_info,
    sigma = c(sigma_c, sigma_E), 
    criteria.success = c(0, 0.9, log(2.5), 0.5),
    criteria.futility = c(log(1.43), 0.5, log(2.5), 0.9),
    prior.difference = c(delta_prior, n_adj, n_adj) 
  )
  
  sim_inf <- gsbSimulation(
    truth = c(-1, delta, 2), 
    type.update = "treatment effect",
    method = "numerical integration"
  )
  
  # Compute operating characteristics
  x_inf <- gsb(design = design_inf, simulation = sim_inf)
  
  sub <- x_inf$OC[,"delta"] != -1
  xout <- x_inf$OC[sub,]
  xout$delta <- red * 100 # Convert back to % for the x-axis
  OCout_inf <- rbind(OCout_inf, xout)
}



# Let's report

# Create 'gsbMainOut' objects to use tab()
x_noninf_final <- x1;    x_noninf_final$OC <- OCout_noninf
x_inf_final    <- x_inf; x_inf_final$OC    <- OCout_inf



table_noninf <- tab(x_noninf_final, what = "cumulative all")
table_inf    <- tab(x_inf_final, what = "cumulative all")

T1E_noninf <- table_noninf[1, "stage2.cum.suc"]
Power_noninf <-  table_noninf[-1, "stage2.cum.suc"]

T1E_inf <- table_inf[1, "stage2.cum.suc"]
Power_inf <-  table_inf[-1, "stage2.cum.suc"]


N_physical_uninfo <- c(patients_uninfo*stages,patients_uninfo*stages*2)
N_physical_info <- c(patients_info*stages,patients_info*stages*2)

res_noninf <- compute_euii(table_total = table_noninf, N = 30, Pr_H1 = 0.1)
res_inf    <- compute_euii(table_total = table_inf,    N = N_physical_info, Pr_H1 = 0.1)

n_noninf <- tab(x_noninf_final, what = "sample size")
n_inf    <- tab(x_inf_final,    what = "sample size")

colnames(n_noninf)[1] <- "reduction" ; colnames(n_inf)[1] <- "reduction" 

# TABLE for Sample Sizes Comparison
N_table <- merge(
  n_noninf, 
  n_inf, 
  by = "reduction", 
  suffixes = c("_NonInf", "_Inf")
)

N_table_clean <- N_table %>%
  select(reduction, `stage 2_NonInf`, `stage 2_Inf`) %>%

  rename(
    `True Reduction (%)` = reduction,
    `E[N] Non-Informative` = `stage 2_NonInf`,
    `E[N] Informative` = `stage 2_Inf`
  )

rownames(N_table_clean) <- NULL




# Construct Comparison Dataframe
comparison_df <- data.frame(
  Delta  = res_noninf$result_df$Delta,
  EUII_NonInf      = res_noninf$result_df$EUII_Exact,
  EUII_Inf         = res_inf$result_df$EUII_Exact,
  ESS_Plus_NonInf  = res_noninf$effective_N$effective_N_plus,
  ESS_Plus_Inf     = res_inf$effective_N$effective_N_plus,
  ESS_Minus_NonInf = res_noninf$effective_N$effective_N_minus,
  ESS_Minus_Inf    = res_inf$effective_N$effective_N_minus
)
comparison_df
additional_df <- data.frame( 
    Delta             = res_noninf$result_df$Delta,
    LR_plus_inf       = res_inf$LR_df$LR_pos,
    LR_plus_noninf    = res_noninf$LR_df$LR_pos,
    LR_minus_inf      = res_inf$LR_df$LR_neg,
    LR_minus_noninf   = res_noninf$LR_df$LR_neg,
    Power_noninf = Power_noninf, 
    Power_inf = Power_inf
)


# ------------------------------------------------ 
# PLOTS VERY INTERESTING
# -------------------------------------------------

# Define  prior probabilities of H1
priors_to_test <- c(0.01, 0.1, 0.9) 
all_plot_data <- list()

# Loop through each P(H1) and compute metrics for both designs
for(p in priors_to_test) {
  
  # Compute for this specific P(H1)
  res_noninf <- compute_euii(table_total = table_noninf, N = N_physical_uninfo, Pr_H1 = p, exact = TRUE)
  res_inf    <- compute_euii(table_total = table_inf,    N = N_physical_info, Pr_H1 = p, exact = TRUE)
  

  res_noninf$result_df
  # Build a clean dataframe for Non-Informative
  df_noninf <- data.frame(
    Delta    = res_noninf$result_df$Delta,
    Design   = "Non-Informative",
    Prior_H1 = paste0("P(H1) = ", p),
    EUII     = res_noninf$result_df$EUII_Exact, 
    effective_N_Plus  = res_noninf$effective_N$effective_N_plus,
    effective_N_Minus = res_noninf$effective_N$effective_N_minus
  )
  
  # Build a clean dataframe for Informative
  df_inf <- data.frame(
    Delta    = res_inf$result_df$Delta,
    Design   = "Informative",
    Prior_H1 = paste0("P(H1) = ", p),
    EUII     = res_inf$result_df$EUII_Exact,
    effective_N_Plus  = res_inf$effective_N$effective_N_plus,
    effective_N_Minus = res_inf$effective_N$effective_N_minus
  )
  
  # Append to our list
  all_plot_data[[length(all_plot_data) + 1]] <- rbind(df_noninf, df_inf)
}

#  Combine everything into one master dataframe
master_df <- bind_rows(all_plot_data)

# Factor the Prior_H1 column so it orders nicely in the legend
master_df$Prior_H1 <- factor(master_df$Prior_H1, levels = paste0("P(H1) = ", priors_to_test))

# Prepare LR plotting dataframe (one curve per design for LR+ and LR-)
lr_plot_df <- additional_df %>%
  select(Delta, LR_plus_inf, LR_plus_noninf, LR_minus_inf, LR_minus_noninf) %>%
  pivot_longer(
    cols = -Delta,
    names_to = "metric",
    values_to = "LR"
  ) %>%
  mutate(
    Design = case_when(
      metric %in% c("LR_plus_inf", "LR_minus_inf") ~ "Informative",
      metric %in% c("LR_plus_noninf", "LR_minus_noninf") ~ "Non-Informative"
    ),
    LR_Type = case_when(
      metric %in% c("LR_plus_inf", "LR_plus_noninf") ~ "LR+",
      metric %in% c("LR_minus_inf", "LR_minus_noninf") ~ "LR-"
    )
  )


# -----------------------------------------------------------------
# PLOT
# -----------------------------------------------------------------

# Base theme to keep things looking clean
my_theme <- theme_minimal() +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.title = element_text(face = "bold")
  )

# EUII Plot
p_euii <- ggplot(master_df, aes(x = Delta, y = EUII, color = Design, linetype = Prior_H1, shape = Prior_H1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Informative" = "royalblue", "Non-Informative" = "firebrick")) +
  labs(title = "EUII", x = "% Reduction", y = "EUII", 
       color = "Design", linetype = "Prior Belief", shape = "Prior Belief") +
  my_theme+ theme(legend.position = "none") # Hide legend to avoid repeating it

# ESS+ Plot
p_enp <- ggplot(master_df, aes(x = Delta, y = effective_N_Plus, color = Design, linetype = Prior_H1, shape = Prior_H1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Informative" = "royalblue", "Non-Informative" = "firebrick")) +
  labs(title = "Effective Sample Size for Success (ESS+)", x = "% Reduction", y = "ESS") +
  my_theme 

# ESS- Plot
p_enm <- ggplot(master_df, aes(x = Delta, y = effective_N_Minus, color = Design, linetype = Prior_H1, shape = Prior_H1)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Informative" = "royalblue", "Non-Informative" = "firebrick")) +
  labs(title = "Effective Sample Size for Non-Significant (ESS-)", x = "% Reduction", y = "ESS") +
  my_theme + theme(legend.position = "none")

# LR+ Plot
p_lr_plus <- ggplot(
  lr_plot_df %>% filter(LR_Type == "LR+"),
  aes(x = Delta, y = LR, color = Design)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Informative" = "royalblue", "Non-Informative" = "firebrick")) +
  labs(title = "Positive Likelihood Ratio (LR+)", x = "% Reduction", y = "LR+") +
  my_theme + theme(legend.position = "none")

# LR- Plot
p_lr_minus <- ggplot(
  lr_plot_df %>% filter(LR_Type == "LR-"),
  aes(x = Delta, y = LR, color = Design)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("Informative" = "royalblue", "Non-Informative" = "firebrick")) +
  labs(title = "Negative Likelihood Ratio (LR-)", x = "% Reduction", y = "LR-") +
  my_theme

# Combine them vertically
p_final <- patchwork::wrap_plots(p_euii, p_enp, p_enm, ncol = 1)
 # p_euii / p_enp / p_enm horiz terrible to see

# Combined LR panel
p_lr_final <- patchwork::wrap_plots(p_lr_plus, p_lr_minus, ncol = 1)









# -----------------------------------------------------------------
#  PLOTS FROM PAPER
# -----------------------------------------------------------------

#  Create the Non-Informative Plot 
p1_lat <- plot(x_noninf_final, what = "cumulative all")
p1_lat$main <- "Non-Informative Prior"
p1_lat$ylab <- "Cumulative prob. of stopping"
p1_lat$xlab <- "" # Hide x-axis for the top row
p1_lat$panel.args.common$col <- 1 
p1_lat$panel.args.common$lty <- 1:3
p1_lat$legend$bottom$args$key$lines$col <- 1 
p1_lat$legend$bottom$args$key$lines$lty <- 1:3

#  Create the Informative Plot
p2_lat <- plot(x_inf_final, what = "cumulative all")
p2_lat$main <- "Informative Prior"
p2_lat$ylab <- "Cumulative prob. of stopping"
p2_lat$xlab <- "% relative reduction"
p2_lat$panel.args.common$col <- 1 
p2_lat$panel.args.common$lty <- 1:3
p2_lat$legend$bottom$args$key$lines$col <- 1 
p2_lat$legend$bottom$args$key$lines$lty <- 1:3

# Convert Lattice to Patchwork-compatible objects
# we use as.ggplot to wrap the lattice objects
p1_wrap <- ggplotify::as.ggplot(p1_lat)
p2_wrap <- ggplotify::as.ggplot(p2_lat)

# Use patchwork explicitly so composition works even if package is not attached
final_oc_plot <- patchwork::wrap_plots(p1_wrap, p2_wrap, ncol = 1)

# He does it smartly by using same legend .... but whatever 
