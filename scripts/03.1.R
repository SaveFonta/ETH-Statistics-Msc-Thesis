#
# Script for creating the First Example 
#

#Load the compute_euii function

source("00_compute_EUII_2.R")

#define variables: 
mysigma <- 88
mycrit.success <- c(0,0.95,50,0.5)
mycrit.fut <- c(40,0.9)
alternative_deltas <- c(0,40,50,60,70) #those are the possible delta values, where delta = mean_control - mean_treatment, so a positive delta means the treatment is better than control.



# we dedfine bias = Real control value - Prior control value 


change_truth <- function(bias) { 
  
  # ---------------------------------------------------------------------
  # Now I write it using grid.type = "manually": 
  # Create an exact N x 2 matrix for the simulations. 
  # Column 1: True Control Mean. Column 2: True Treatment Mean.
  # ---------------------------------------------------------------------
  true_control <- 49 + bias
  true_treatment <- true_control + alternative_deltas
  truth_matrix <- cbind(true_control, true_treatment)
  
  # -------------------
  # Design 1 -> Informative Sequential
  design_informative <- gsbDesign(nr.stages = 2, patients = c(10,20), sigma=mysigma,
                                  criteria.success = mycrit.success, criteria.futility = mycrit.fut,
                                  prior.control = c(49,20), 
                                  prior.treatment = c(49,0)) 
  
  simulation_informative <- gsbSimulation(truth = truth_matrix,
                                          grid.type = "manually", 
                                          type.update = "per arm", 
                                          nr.sim = 100000, warnings.sensitivity = 500, seed=1)
  
  result_informative <- gsb(design=design_informative, simulation=simulation_informative)
  table_seq_informative <- tab(result_informative, what = "cumulative all")
  
  # -------------------
  # Design 2 -> Uninformative Sequential
  design_uninformative <- gsbDesign(nr.stages = 2, patients = c(20,20), sigma=mysigma,
                                    criteria.success = mycrit.success, criteria.futility = mycrit.fut,
                                    prior.control = c(49,0), prior.treatment = c(49,0))
  
  simulation_uninformative <- gsbSimulation(truth = truth_matrix,
                                            grid.type = "manually", 
                                            type.update = "per arm",
                                            nr.sim = 100000, warnings.sensitivity = 500, seed=1)
  
  result_uninformative <- gsb(design=design_uninformative, simulation=simulation_uninformative)
  table_seq_uninformative <- tab(result_uninformative, what = "cumulative all")
  
  # Calculate EUII for sequential
  res_seq_info <- compute_euii(table_seq_informative, N= c(30, 60), prior_N = 0, Pr_H1 = 0.1, exact = TRUE)
  res_seq_uninfo <- compute_euii(table_seq_uninformative, N = c(40,80), prior_N = 0, Pr_H1 = 0.1, exact = TRUE)
  
  # -------------------
  # Design 3 -> Uninformative Fixed
  design_fixed_uninformative <- gsbDesign(nr.stages = 1, patients = c(40,40), sigma=mysigma,
                                          criteria.success = mycrit.success, criteria.futility = mycrit.fut,
                                          prior.control = c(49,0), prior.treatment = c(49,0))
  
  simulation_fixed_uninformative <- gsbSimulation(truth = truth_matrix, 
                                                  grid.type = "manually", 
                                                  type.update = "per arm",
                                                  nr.sim=100000, warnings.sensitivity = 500, seed=1)
  
  result_fixed_uninformative <- gsb(design=design_fixed_uninformative, simulation=simulation_fixed_uninformative)
  table_fixed_uninformative <- tab(result_fixed_uninformative, what = "cumulative all")
  
  # -------------------
  # Design 4 -> Informative Fixed
  design_fixed_informative <- gsbDesign(nr.stages = 1, patients = c(20,40), sigma=mysigma,
                                        criteria.success = mycrit.success, criteria.futility = mycrit.fut,
                                        prior.control = c(49,20), 
                                        prior.treatment = c(49,0)) 
  
  simulation_fixed_informative <- gsbSimulation(truth = truth_matrix, 
                                                grid.type = "manually", 
                                                type.update = "per arm",
                                                nr.sim=100000, warnings.sensitivity = 500, seed=1)
  
  result_fixed_informative <- gsb(design=design_fixed_informative, simulation=simulation_fixed_informative)
  table_fixed_informative <- tab(result_fixed_informative, what = "cumulative all")
  
# --- Extraction ---
  delta_vals <- alternative_deltas[-1]
  
  T1E_fixed_uninfo <- table_fixed_uninformative[1, "stage1.suc"]
  Power_fixed_uninfo <- table_fixed_uninformative[-1, "stage1.suc"]
  EUII_fixed_uninfo <- ((Power_fixed_uninfo / T1E_fixed_uninfo) / 
                        ((1 - Power_fixed_uninfo) / (1 - T1E_fixed_uninfo))) ^ (1/80)
  
  T1E_fixed_info <- table_fixed_informative[1, "stage1.suc"]
  Power_fixed_info <- table_fixed_informative[-1, "stage1.suc"]
  EUII_fixed_info <- ((Power_fixed_info / T1E_fixed_info) / 
                      ((1 - Power_fixed_info) / (1 - T1E_fixed_info))) ^ (1/60)
  
  # --- Build Dataframes ---
  df_fixed_uninfo <- data.frame(Delta = delta_vals, EUII = EUII_fixed_uninfo, Design = "Fixed Uninformative (N=80)")
  df_fixed_info <- data.frame(Delta = delta_vals, EUII = EUII_fixed_info, Design = "Fixed Informative (N=60)")
  df_seq_info <- data.frame(Delta = res_seq_info$result_df$Delta, EUII = res_seq_info$result_df$EUII_Exact, Design = "Sequential Informative (Max N=60)")
  df_seq_uninfo <- data.frame(Delta = res_seq_uninfo$result_df$Delta, EUII = res_seq_uninfo$result_df$EUII_Exact, Design = "Sequential Uninformative (Max N=80)")
  
  df_merged <- rbind(df_fixed_uninfo, df_fixed_info, df_seq_info, df_seq_uninfo)
  df_merged$Bias <- bias
  

  # ---- Build df with additional values, such as power and T1E for each design
  T1E_S.I <- table_seq_informative [1, "stage2.suc"]
  T1E_S.U <- table_seq_uninformative [1, "stage2.suc"]
  Power_S.I <- table_seq_informative [-1, "stage2.suc"]
  Power_S.U <- table_seq_uninformative [-1, "stage2.suc"]

  T1E <- c(T1E_sequential.info = T1E_S.I, T1E_sequential.noninfo = T1E_S.U,  T1E_fixed.info = T1E_fixed_info, T1E_fixed.uninfo = T1E_fixed_uninfo)
  Power <- cbind(
    Delta = delta_vals,
    Power_sequential.info = Power_S.I,
    Power_sequential.noninfo = Power_S.U,
    Power_fixed.info = Power_fixed_info,
    Power_fixed.uninfo = Power_fixed_uninfo
  )

  res <- list(df_merged = df_merged, T1E = T1E, Power = Power, res_seq_info = res_seq_info, res_seq_uninfo = res_seq_uninfo )
  return(res)
}


# classic approach with bias = 0
df_merged <- change_truth(0)


final_plot_CRON <- ggplot(df_merged$df_merged, aes(x = Delta, y = EUII, 
                                      color = Design, shape = Design, linetype = Design)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  # Using a colorblind-friendly palette
  scale_color_manual(values = c("Fixed Uninformative (N=80)" = "#E69F00", 
                                "Fixed Informative (N=60)" = "#D55E00",
                                "Sequential Uninformative (Max N=80)" = "#56B4E9", 
                                "Sequential Informative (Max N=60)" = "#0072B2")) +
  theme_bw(base_size = 14) +
  labs(
    x = "True Treatment Effect",
    y = "EUII"
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    panel.grid.minor = element_blank()
  )


## Sample size df building

info <- df_merged$res_seq_info$N[, 2:3] 
uninfo <- df_merged$res_seq_uninfo$N[,2:3]
names(info) <- c("E_N_+_info", "E_N_-_info")
names(uninfo) <- c("E_N_+_uninfo", "E_N_-_uninfo")
delta <- df_merged$res_seq_uninfo$result_df$Delta
N_data_CRON <- cbind(delta, info, uninfo)

library(xtable)
colnames(N_data_CRON) <- c("$\\Delta$", "$E[N_+]^{\\text{info}}$", "$E[N_-]^{\\text{info}}$", "$E[N_+]^{\\text{uninfo}}$", "$E[N_-]^{\\text{uninfo}}$")
N_data_CRON <- xtable(N_data_CRON, digits = c(0, 0,2,2,2,2))





#####################################################################à
# What if we change the bias value??????????

bias_values <- seq(-20, 20, by = 5)
raw_results <- lapply(bias_values, function(b) {
  res <- change_truth(b)
  
  # Inject the bias value into the Power and T1E tables before returning
  # so we don't lose track of which bias generated them!
  res$Power <- as.data.frame(res$Power)
  res$Power$Bias <- b
  
  res$T1E <- as.data.frame(t(res$T1E)) # Transpose named vector to 1-row dataframe
  res$T1E$Bias <- b
  
  return(res)
})

df_euii_total <- bind_rows(lapply(raw_results, function(x) x$df_merged))
df_power_total <- bind_rows(lapply(raw_results, function(x) x$Power))
df_t1e_total <- bind_rows(lapply(raw_results, function(x) x$T1E))


# Filter for a target Alternative Hypothesis (e.g., Delta = 50)
df_plot_euii <- df_euii_total  |>  
                  filter(Delta == 50)  |> 
                  filter(is.finite(EUII))

# Plot the EUII vs. Prior Bias
EUII_CRON <- ggplot(df_plot_euii, aes(x = Bias, y = EUII, color = Design, linetype = Design, shape = Design)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("Fixed Uninformative (N=80)" = "#E69F00", 
                                "Fixed Informative (N=60)" = "#D55E00",
                                "Sequential Uninformative (Max N=80)" = "#56B4E9", 
                                "Sequential Informative (Max N=60)" = "#0072B2")) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "darkgray", linewidth = 1) +
  theme_bw(base_size = 14) +
  labs(
    title = "EUII Sensitivity to Prior-Data Conflict",
    subtitle = expression("True Effect Size " * Delta * " = 50"),
    x = "Prior Bias (True Control Mean - Prior Mean)",
    y = "EUII"
  ) +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank())





library(tidyr)

# Reshape T1E for ggplot
df_t1e_long <- df_t1e_total %>%
  pivot_longer(cols = starts_with("T1E"), names_to = "Design", values_to = "Alpha")

T1E_CRON <- ggplot(df_t1e_long, aes(x = Bias, y = Alpha, color = Design)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3.5) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + # Standard 5% line
  theme_bw(base_size = 14) +
  labs(  x = "Prior Bias", y = "Type-I Error Rate")




# Filter Power for Delta = 50 and reshape
df_power_long <- df_power_total %>% 
  filter(Delta == 50) %>%
  pivot_longer(cols = starts_with("Power"), names_to = "Design", values_to = "Power_Val")

Power_plot_CRON <- ggplot(df_power_long, aes(x = Bias, y = Power_Val, color = Design)) +
  geom_line(linewidth = 1.2) + geom_point(size = 3.5) +
  theme_bw(base_size = 14) +
  labs(x = "Prior Bias", y = "Statistical Power")







library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Extract the components directly from your raw_results list
df_components <- do.call(rbind, lapply(raw_results, function(res) {
  
  b <- res$T1E$Bias
  
  # Extract Type I Error and Power for the Sequential Informative design
  alpha <- res$T1E$T1E_sequential.info
  power <- res$Power$Power_sequential.info[res$Power$Delta == 50]
  
  # Calculate Likelihood Ratios (adding a tiny epsilon to prevent division by zero)
  epsilon <- 1e-6
  LR_plus <- power / max(alpha, epsilon)
  LR_minus <- (1 - power) / max(1 - alpha, epsilon)
  
  # Extract Expected Sample Sizes (E[N+] and E[N-]) for Delta = 50
  # We find the row index where Delta == 50 in the sequential info object
  idx <- which(res$res_seq_info$result_df$Delta == 50)
  
  # Based on your setup, column 2 is E[N+] and column 3 is E[N-]
  EN_plus <- res$res_seq_info$N[idx, 2]
  EN_minus <- res$res_seq_info$N[idx, 3]
  
  data.frame(
    Bias = b,
    LR_plus = LR_plus,
    LR_minus = LR_minus,
    EN_plus = as.numeric(EN_plus),
    EN_minus = as.numeric(EN_minus)
  )
}))

# 2. Reshape the data for a multi-panel ggplot
df_plot_components <- df_components %>%
  pivot_longer(
    cols = c("LR_plus", "LR_minus", "EN_plus", "EN_minus"), 
    names_to = "Metric", 
    values_to = "Value"
  ) %>%
  mutate(
    # Group metrics into two categories for faceting
    Category = ifelse(grepl("LR", Metric), "Evidentiary Value (Likelihood Ratios)", "Cost (Expected Sample Size)"),
    # Rename for cleaner plot legends
    Metric = case_when(
      Metric == "LR_plus"  ~ "LR+ (Evidence for Alternative)",
      Metric == "LR_minus" ~ "LR- (Evidence for Null)",
      Metric == "EN_plus"  ~ "E[N+] (Sample Size if Significant)",
      Metric == "EN_minus" ~ "E[N-] (Sample Size if Non-Sig.)"
    )
  )

# 3. Create the Decoupled Plot
decoupled_plot <- ggplot(df_plot_components, aes(x = Bias, y = Value, color = Metric, linetype = Metric)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  # Free Y-scales because LRs can be in the hundreds, while N is bounded at 60
  facet_wrap(~ Category, scales = "free_y", ncol = 1) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "darkgray", linewidth = 1) +
  scale_color_manual(values = c(
    "LR+ (Evidence for Alternative)" = "#0072B2", 
    "LR- (Evidence for Null)" = "#56B4E9",
    "E[N+] (Sample Size if Significant)" = "#D55E00", 
    "E[N-] (Sample Size if Non-Sig.)" = "#E69F00"
  )) +
  theme_bw(base_size = 14) +
  labs(
    title = "Mechanics of EUII Drop: Likelihood Ratios vs. Sample Size",
    subtitle = "Sequential Informative Design (True Treatment Effect = 50)",
    x = "Prior Bias (True Control Mean - Prior Mean)",
    y = "Metric Value"
  ) +
  theme(
    legend.position = "bottom", 
    legend.direction = "vertical", 
    legend.title = element_blank(),
    strip.text = element_text(face = "bold", size = 12),
    panel.grid.minor = element_blank()
  )

# Display the plot
print(decoupled_plot)