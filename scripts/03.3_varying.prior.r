
# -----------------------------------------------------------------
# Variant of Example 3.3
# -----------------------------------------------------------------

# Not anymore a sequential bayesian DC design, but now a  fIXED DC bayesian DESIGN

source("00_compute_EUII_2.R")
library(xtable)
library(dplyr)
library(tidyr)
library(ggplot2)

ex.3.function <- function(
  patients_uninfo = NULL, 
  patients_info = NULL, 
  prior_patients = 10, # 10 +10 = 20
  prior_reduction = 0.75, 
  mu_c = 14.8,
  kappa = 2,
  stages = 1,
  pr_h1 = 0.1
) {

if (is.null(patients_uninfo)) {
  patients_uninfo <- if (stages == 1) 30 else 15
}
if (is.null(patients_info)) {
  patients_info <- if (stages == 1) 30 else 15
}

# patient vector based on how many patients we have
patients_vec_uninfo <- (1:stages) * patients_uninfo
patients_vec_info <- (1:stages) * patients_info



rate_c <- log(mu_c) #log-rate for control
sigma_c <- sqrt(kappa + (1 / mu_c))  #control is always the same 
redscen <- c(0, seq(0.5, 0.9, 0.05)) #reduction scenarios from 50% to 90% --> 

# QUESTION: do we even care about reduction scenarios lower than 50%? Why should we inspect them 

# Nrormal approx for log(negative binomial) is ->  log (r) = N (log(mu), 1/(mu*n) + k/n)  where r is the estimate of mu


  OCout_noninf <- NULL
  OCout_inf <- NULL


  # --- Prior parameters (kept constant across truth scenarios) ---
  mu_prior <- mu_c * (1 - prior_reduction)
  sigma_prior_E <- sqrt(kappa + (1 / mu_prior))
  V_target <- (sigma_c^2 + sigma_prior_E^2) / prior_patients
  delta_prior <- rate_c - log(mu_prior)



#loop since sigma_E is tied to mu_E --> we need to recompute it for each mu_E scenario! 
 for(i in seq_along(redscen)){
    red <- redscen[i]
    mu_E <- mu_c * (1 - red)
    rate_E <- log(mu_E)
    delta <- rate_c - rate_E
    sigma_E <- sqrt(kappa + (1 / mu_E))
    
    # ---------------------------------------------------------
    #  Non-Informative Design
    # ---------------------------------------------------------
    design1 <- gsbDesign(
      nr.stages = stages,
      patients = patients_vec_uninfo, 
      sigma = c(sigma_c, sigma_E), 
      criteria.success = c(0, 0.9, log(2.5), 0.5), 
      criteria.futility = c(log(1.43), 0.5, log(2.5), 0.9),
      prior.difference = "non-informative" 
    )
    
    simulation1 <- gsbSimulation(
      truth = c(-1, delta, 2), 
      type.update = "treatment effect", 
      method = "numerical integration"
    )
    
    x1 <- gsb(design = design1, simulation = simulation1)
    
    # Extract only the targeted delta (ignore -1)
    sub1 <- x1$OC[,"delta"] != -1
    xout1 <- x1$OC[sub1, , drop=FALSE]
    xout1$delta <- red * 100 
    OCout_noninf <- rbind(OCout_noninf, xout1)
    
    # ---------------------------------------------------------
    #  Informative Design
    # ---------------------------------------------------------

    n_adj <- (sigma_c^2 + sigma_E^2) / V_target
    
    design_inf <- gsbDesign(
      nr.stages = stages,
      patients = patients_vec_info,
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
    
    x_inf <- gsb(design = design_inf, simulation = sim_inf)
    
    sub_inf <- x_inf$OC[,"delta"] != -1
    xout_inf <- x_inf$OC[sub_inf, , drop=FALSE]
    xout_inf$delta <- red * 100 
    OCout_inf <- rbind(OCout_inf, xout_inf)
  }
  
  # ---------------------------------------------------------
  # Report
  # ---------------------------------------------------------
  x_noninf_final <- x1
  x_noninf_final$OC <- OCout_noninf
  
  x_inf_final <- x_inf
  x_inf_final$OC <- OCout_inf
  
  table_noninf <- tab(x_noninf_final, what = "cumulative all")
  table_inf    <- tab(x_inf_final, what = "cumulative all")
  
  if (stages == 1 ){
  col_name <- paste0("stage", stages, ".suc")
  } else if (stages == 2){
  col_name <- paste0("stage", stages, ".cum.suc")
  }

  T1E_noninf <- table_noninf[1, col_name]
  Power_noninf <- table_noninf[-1, col_name]
  
  T1E_inf <- table_inf[1, col_name]
  Power_inf <- table_inf[-1, col_name]
  
  if (stages == 1 ){
  N_physical_uninfo <- patients_uninfo * 2
  N_physical_info <- patients_info * 2
  } else if (stages == 2 ){
    N_physical_uninfo <- c(patients_uninfo * 2, patients_uninfo * 4)
    N_physical_info <- c(patients_info * 2, patients_info * 4)
  }


  res_noninf <- compute_euii(table_total = table_noninf, N = N_physical_uninfo, Pr_H1 = pr_h1)
  res_inf    <- compute_euii(table_total = table_inf,    N = N_physical_info, Pr_H1 = pr_h1)
  
  # Construct final return list
  comparison_df <- data.frame(
    Reduction_Pct = res_noninf$result_df$Delta,
    EUII_NonInf      = res_noninf$result_df$EUII_Exact,
    EUII_Inf         = res_inf$result_df$EUII_Exact
  )
  
  diagnostics_df <- data.frame(
    Reduction_Pct = res_noninf$result_df$Delta,
    Power_NonInf = Power_noninf,
    Power_Inf = Power_inf,
    T1E_NonInf = rep(T1E_noninf, length(Power_noninf)),
    T1E_Inf = rep(T1E_inf, length(Power_inf)),
    DOR_NonInf = res_noninf$DOR_df$DOR,
    DOR_Inf = res_inf$DOR_df$DOR
  )

  sample_size_df <- data.frame(
    Reduction_Pct = res_noninf$result_df$Delta,
    Eff_N_Plus_NonInf = res_noninf$effective_N$effective_N_plus,
    Eff_N_Minus_NonInf = res_noninf$effective_N$effective_N_minus,
    Eff_N_Plus_Inf = res_inf$effective_N$effective_N_plus,
    Eff_N_Minus_Inf = res_inf$effective_N$effective_N_minus
  )
  
  return(list(
    EUII_Comparison = comparison_df,
    Diagnostics = diagnostics_df,
    SampleSizes = sample_size_df
  ))
}


results_entusiastic <- ex.3.function(
  prior_reduction = 0.75, 
  prior_patients = 10, stages = 1
)
print(results_entusiastic$Diagnostics)



#more skeptical prior
results_skeptical <- ex.3.function(
  prior_reduction = 0.3, 
  prior_patients = 10
)
results_skeptical$EUII_Comparison



























# -----------------------------------------------------------------
# Let's evaluate a grid of prior
# -----------------------------------------------------------------


ex.3.function_grid.prior <- function(stages = 1, patients_uninfo = NULL, patients_info = NULL) {

  if (is.null(patients_uninfo)) {
    patients_uninfo <- if (stages == 1) 30 else 15
  }
  if (is.null(patients_info)) {
    patients_info <- if (stages == 1) 30 else 15
  }
  
  prior_beliefs <- c(0, 0.20, 0.40, 0.60, 0.70, 0.80)
  names(prior_beliefs) <- c("0%", "20%", "40%", "60%", "70%", "80%")
  
  target_effects <- c(60, 70, 80, 90)


  results_list <- list()
  
  for (i in seq_along(prior_beliefs)) {
    prior_val <- prior_beliefs[i]
    prior_name <- names(prior_beliefs)[i]
    
    cat("Running simulation for Prior:", prior_name, "...\n")
    
    suppressWarnings({
      sim_res <- ex.3.function(
        stages = stages,
        prior_reduction = prior_val,
        prior_patients = 10,
        patients_uninfo = patients_uninfo, 
        patients_info = patients_info
      )
    })
    
full_res <- cbind(sim_res$Diagnostics, 
                      EUII_NonInf = sim_res$EUII_Comparison$EUII_NonInf, 
                      EUII_Inf = sim_res$EUII_Comparison$EUII_Inf,
                      Eff_N_Plus_NonInf = sim_res$SampleSizes$Eff_N_Plus_NonInf,
                      Eff_N_Minus_NonInf = sim_res$SampleSizes$Eff_N_Minus_NonInf,
                      Eff_N_Plus_Inf = sim_res$SampleSizes$Eff_N_Plus_Inf,
                      Eff_N_Minus_Inf = sim_res$SampleSizes$Eff_N_Minus_Inf)
    
    # Extract T1E (Row 1 where True Reduction is 0)
    t1e_inf <- full_res$T1E_Inf[1]
    t1e_noninf <- full_res$T1E_NonInf[1]
    
    # Extract Power and EUII for our specified Target Effects
    wanted_rows <- full_res  |>  
                    filter(round(Reduction_Pct) %in% target_effects)  |> 
                    mutate (
                      Prior_Reduction_Assumed = prior_val * 100,
                      True_Effect = paste0("True Reduction: ", round(Reduction_Pct), "%"),
                      T1E_Inf = t1e_inf,
                      T1E_NonInf = t1e_noninf
                    )  |> 
                    select(
                      Prior_Reduction_Assumed, True_Effect, 
                                              T1E_Inf, T1E_NonInf, 
                                              Power_Inf, Power_NonInf, 
                                              EUII_Inf, EUII_NonInf,
                                              Eff_N_Plus_Inf, Eff_N_Plus_NonInf,
                                              Eff_N_Minus_Inf, Eff_N_Minus_NonInf)

        # add the created df to the list
        results_list[[i]] <- wanted_rows
  }
  return(bind_rows(results_list))
}

# Run the simulations
grid_results <- ex.3.function_grid.prior(patients_uninfo = 30, patients_info = 30)

# Reshape for the 3x4 Grid Plot
plot_data <- grid_results %>%
  pivot_longer(
    cols = c(T1E_Inf, T1E_NonInf, Power_Inf, Power_NonInf, EUII_Inf, EUII_NonInf),
    names_to = c("Metric", "Design"),
    names_sep = "_"
  ) %>%
  mutate(
    Design = ifelse(Design == "Inf", "Informative Prior", "Non-Informative"),
    Metric = factor(Metric, levels = c("T1E", "Power", "EUII")),
    True_Effect = factor(True_Effect, levels = c("True Reduction: 60%", 
                                                 "True Reduction: 70%", 
                                                 "True Reduction: 80%", 
                                                 "True Reduction: 90%"))
  )

# 3x4 Grid Plot
grid_plot <- ggplot(plot_data, aes(x = Prior_Reduction_Assumed, y = value, 
                                                      color = Design, linetype = Design)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  facet_grid(Metric ~ True_Effect, scales = "free_y") +
  scale_color_manual(values = c("Informative Prior" = "#D55E00", 
                                "Non-Informative" = "#0072B2")) +
  scale_linetype_manual(values = c("Informative Prior" = "solid", 
                                   "Non-Informative" = "dashed")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Operating Characteristics & EUII for Fixed Design",
    x = "Assumed Prior Reduction (%)",
    y = "Value"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "#e5e5e5"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )







# We are in a fixed design. So N is not r.v., but N = 50. 
#So both informative and non informative are being elevated to 1/N both. And the difference in EUII is made
#only from the DOR.

# Maybe there is some literature regarding that. But what happens is that the DOR is higher for the noninformative desig if the prior is too entusiastic. 
# evem when prior mean = true mean!!!

# This happenns cause the T1e with a prior so entusiastic is inflated, leading in turn to a reduced LR_+ and then to a reduced DOR

# A way to dedinflate the T1E would be to increase the threshold of the bayesian criteria for success

# increasing real sample size of the prior (e.g setting patients_info = 15)  inflates even more T1E
 # while decreasing it, patients_info = 5, reduces T1E.

# Maybe we could compare the designs at the same T1E, but this would mean using some kind of uniroot on the gsb function to find the best treshold 
# or sample size that lead to T1E = 0.05 and then the comparison is fair.... 


# Another good idea to decrease T1E --> use a mixture prior? But how to implement?








# Not relevant check to see if compute_euii works fine for fixed design 
grid_results |> 
  mutate(
    # Informative Check
    EUII_check = ( (Power_Inf / (1 - Power_Inf)) / (T1E_Inf / (1 - T1E_Inf)) )^(1/50),
    
    # Non-Informative Check
    EUII_check_noninf = ( (Power_NonInf / (1 - Power_NonInf)) / (T1E_NonInf / (1 - T1E_NonInf)) )^(1/50)
  ) |>
  # Optional: Select just the relevant columns to easily compare them side-by-side
  select(True_Effect, Prior_Reduction_Assumed, 
         EUII_Inf, EUII_check, 
         EUII_NonInf, EUII_check_noninf)



























#################################################


# -------------------------------------------
#              Sequential Design
# -------------------------------------------
grid_results_seq <- ex.3.function_grid.prior(stages =2)


# Reshape the Sequential Data
plot_data_seq <- grid_results_seq %>%
  pivot_longer(
    cols = c(T1E_Inf, T1E_NonInf, Power_Inf, Power_NonInf, EUII_Inf, EUII_NonInf),
    names_to = c("Metric", "Design"),
    names_sep = "_"
  ) %>%
  mutate(
    Design = ifelse(Design == "Inf", "Informative Prior", "Non-Informative"),
    Metric = factor(Metric, levels = c("T1E", "Power", "EUII")),
    True_Effect = factor(True_Effect, levels = c("True Reduction: 60%", 
                                                 "True Reduction: 70%", 
                                                 "True Reduction: 80%", 
                                                 "True Reduction: 90%"))
  )

# We bind both datasets together to find the global max and min for each Metric
combined_data <- bind_rows(plot_data, plot_data_seq)

y_limits <- combined_data %>%
  group_by(Metric) %>%
  summarize(
    min_val = min(value, na.rm = TRUE),
    max_val = max(value, na.rm = TRUE)
  ) %>%
  pivot_longer(cols = c(min_val, max_val), values_to = "value") %>%
  mutate(
    Prior_Reduction_Assumed = 50, 
    True_Effect = factor("True Reduction: 60%", levels = levels(plot_data_seq$True_Effect)),
    Design = "Informative Prior"
  )

#  Generate the 3x4 Grid Plot for the Sequential Design
grid_plot_seq <- ggplot(plot_data_seq, aes(x = Prior_Reduction_Assumed, y = value, 
                                                  color = Design, linetype = Design)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  # This invisible layer forces the axes to match the global min/max!
  geom_blank(data = y_limits, aes(x = Prior_Reduction_Assumed, y = value)) +
  facet_grid(Metric ~ True_Effect, scales = "free_y") +
  scale_color_manual(values = c("Informative Prior" = "#D55E00", 
                                "Non-Informative" = "#0072B2")) +
  scale_linetype_manual(values = c("Informative Prior" = "solid", 
                                   "Non-Informative" = "dashed")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Operating Characteristics & EUII for Sequential Design",
    x = "Assumed Prior Reduction (%)",
    y = "Value"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "#e5e5e5"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

# Rebuild fixed-design plot using shared y-limits so it is directly comparable to sequential
grid_plot <- ggplot(plot_data, aes(x = Prior_Reduction_Assumed, y = value,
                                   color = Design, linetype = Design)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 1) +
  geom_blank(data = y_limits, aes(x = Prior_Reduction_Assumed, y = value)) +
  facet_grid(Metric ~ True_Effect, scales = "free_y") +
  scale_color_manual(values = c("Informative Prior" = "#D55E00",
                                "Non-Informative" = "#0072B2")) +
  scale_linetype_manual(values = c("Informative Prior" = "solid",
                                   "Non-Informative" = "dashed")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Operating Characteristics & EUII for Fixed Design",
    x = "Assumed Prior Reduction (%)",
    y = "Value"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "#e5e5e5"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

print(grid_plot)
print(grid_plot_seq)



























#################################################


# -------------------------------------------
#              Fixed vs Sequential Design
# -------------------------------------------

#  Merge  Fixed (grid_results) and Sequential (grid_results_seq)
combined_wide <- grid_results %>%
  rename_with(~paste0(., "_fixed"), -c(Prior_Reduction_Assumed, True_Effect)) %>% #renam all cols exctept those two with suffix
  left_join(
    grid_results_seq %>% rename_with(~paste0(., "_seq"), -c(Prior_Reduction_Assumed, True_Effect)),
    by = c("Prior_Reduction_Assumed", "True_Effect") #join by those
  )

#  compute ratios (seq/fix)  
ratio_data <- combined_wide %>%
  mutate(
    # Type I Error Ratios
    T1E_Inf_Ratio = T1E_Inf_seq / T1E_Inf_fixed,
    T1E_NonInf_Ratio = T1E_NonInf_seq / T1E_NonInf_fixed,
    
    # Power Ratios
    Power_Inf_Ratio = Power_Inf_seq / Power_Inf_fixed,
    Power_NonInf_Ratio = Power_NonInf_seq / Power_NonInf_fixed,
    
    # EUII Ratios
    EUII_Inf_Ratio = EUII_Inf_seq / EUII_Inf_fixed,
    EUII_NonInf_Ratio = EUII_NonInf_seq / EUII_NonInf_fixed
  ) %>%
  select(Prior_Reduction_Assumed, True_Effect, 
         T1E_Inf_Ratio, T1E_NonInf_Ratio, 
         Power_Inf_Ratio, Power_NonInf_Ratio, 
         EUII_Inf_Ratio, EUII_NonInf_Ratio)

# longer
plot_ratio <- ratio_data %>%
  pivot_longer(
    cols = ends_with("_Ratio"),
    names_to = c("Metric", "Design"),
    names_pattern = "(.*)_(Inf|NonInf)_Ratio", #  split Metric and Design
    values_to = "Ratio"
  ) %>%
  mutate(
    Design = ifelse(Design == "Inf", "Informative Prior", "Non-Informative"),
    # Rename metrics
    Metric = factor(Metric, levels = c("T1E", "Power", "EUII"), 
                    labels = c("T1E Ratio", "Power Ratio", "EUII Ratio")),
    # Clean any undefined 
    Ratio = ifelse(is.nan(Ratio) | is.infinite(Ratio), NA, Ratio)
  ) %>%
  # Drop NAs 
  filter(!is.na(Ratio))

# Generate the 3x4 grid Plot
ratio_plot <- ggplot(plot_ratio, aes(x = Prior_Reduction_Assumed, y = Ratio, 
                                     color = Design, linetype = Design)) +
  #  reference line at 1.0 (Where Seq = Fixed)
  geom_hline(yintercept = 1, color = "black", linewidth = 0.8) + 
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_grid(Metric ~ True_Effect, scales = "free_y") +
  scale_color_manual(values = c("Informative Prior" = "#D55E00", 
                                "Non-Informative" = "#0072B2")) +
  scale_linetype_manual(values = c("Informative Prior" = "solid", 
                                   "Non-Informative" = "dashed")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Ratio Sequential vs. Fixed Design",
    x = "Assumed Prior Reduction (%)",
    y = "Ratio (Sequential / Fixed)"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "#e5e5e5"),
    strip.text = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )

print(ratio_plot)























# ---------------------------------------------------
# Summary of sample size used
######## --------------------------------------------

plot_data_n <- grid_results_seq %>%
  select(Prior_Reduction_Assumed, True_Effect, 
         Eff_N_Plus_Inf, Eff_N_Plus_NonInf, 
         Eff_N_Minus_Inf, Eff_N_Minus_NonInf) %>%
  #: break the names into "Direction" (Plus/Minus) and "Design" (Inf/NonInf)
  pivot_longer(
    cols = starts_with("Eff_N"),
    names_to = c("Direction", "Design"),
    names_pattern = "Eff_N_(Plus|Minus)_(.*)", 
    values_to = "Effective_N"
  ) %>%
  mutate(
    Effective_N = ifelse(is.nan(Effective_N) | is.infinite(Effective_N), NA, Effective_N),
    
    # Clean up labels for the legend
    Design = ifelse(Design == "Inf", "Informative Prior", "Non-Informative"),
    Direction = ifelse(Direction == "Plus", "Plus", "Minus")
  ) %>%
  # Drop the NAs so ggplot doesn't even try to draw them
  filter(!is.na(Effective_N))


n_plot <- ggplot(plot_data_n, aes(x = Prior_Reduction_Assumed, y = Effective_N, 
                                  color = Design, linetype = Design, shape = Direction)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 3) +
  # Facet by true underlying effect
  facet_wrap(~ True_Effect, scales = "free_y", ncol = 4) +
  scale_shape_manual(
    values = c("Plus" = 16, "Minus" = 17),
    breaks = c("Plus", "Minus"),
    labels = c("E(1/N|+)^-1", "E(1/N|-)^-1")
  ) +
  scale_color_manual(values = c("Informative Prior" = "#D55E00", 
                                "Non-Informative" = "#0072B2")) +
  scale_linetype_manual(values = c("Informative Prior" = "solid", 
                                   "Non-Informative" = "dashed")) +
  theme_bw(base_size = 14) +
  labs(
    title = "Effective Sample Sizes for Sequential Design",
    x = "Assumed Prior Reduction (%)",
    y = "Effective Sample Size"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical", 
    legend.title = element_blank(),
    strip.background = element_rect(fill = "#e5e5e5"),
    strip.text = element_text(face = "bold", size = 11)
  )

print(n_plot)






































