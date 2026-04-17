#
# Script for creating the First Example 
#

#Load the compute_euii function

source("00_compute_EUII_2.R")

#define variables: 
mysigma <- 88
mycrit.success <- c(0,0.95,50,0.5)
mycrit.fut <- c(40,0.9)
alternative_deltas <- c(0,10,20,30,40,50,60,70) #those are the possible delta values, where delta = mean_control - mean_treatment, so a positive delta means the treatment is better than control.
control_value <- 49


# we dedfine bias = Real control value - Prior control value 

# --------------------------------
# Function to create the values
# ---------------------

change_truth <- function(bias) { 
  
  # ---------------------------------------------------------------------
  # Now I write it using grid.type = "manually": 
  # Create an exact N x 2 matrix for the simulations. 
  # Column 1: True Control Mean. Column 2: True Treatment Mean.
  # ---------------------------------------------------------------------
  true_control <- control_value + bias
  true_treatment <- true_control + alternative_deltas
  truth_matrix <- cbind(true_control, true_treatment)
  
  # -------------------
  # Design 1 -> Informative Sequential
  design_informative <- gsbDesign(nr.stages = 2, patients = c(10,20), sigma=mysigma,
                                  criteria.success = mycrit.success, criteria.futility = mycrit.fut,
                                  prior.control = c(control_value,20), 
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

T1E_example1 <- df_merged$T1E
Power_example1 <- df_merged$Power

final_plot_CRON <- ggplot(df_merged$df_merged, aes(x = Delta, y = EUII, 
                                      color = Design, shape = Design, linetype = Design)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.2) +
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

info <- df_merged$res_seq_info$effective_N[, 2:5]
uninfo <- df_merged$res_seq_uninfo$effective_N[, 2:5]
names(info) <- c("N_eff_+_info", "N_eff_-_info", "CV_+_info", "CV_-_info")
names(uninfo) <- c("N_eff_+_uninfo", "N_eff_-_uninfo", "CV_+_uninfo", "CV_-_uninfo")
delta <- df_merged$res_seq_uninfo$result_df$Delta
N_data_CRON <- cbind(delta, info, uninfo)

library(xtable)
colnames(N_data_CRON) <- c(
  "$\\Delta$",
  "$N_{\\text{eff},+}^{\\text{info}}$", "$N_{\\text{eff},-}^{\\text{info}}$",
  "$CV_{+}^{\\text{info}}$", "$CV_{-}^{\\text{info}}$",
  "$N_{\\text{eff},+}^{\\text{uninfo}}$", "$N_{\\text{eff},-}^{\\text{uninfo}}$",
  "$CV_{+}^{\\text{uninfo}}$", "$CV_{-}^{\\text{uninfo}}$"
)
N_data_CRON <- xtable(N_data_CRON, digits = c(0, 0, 2, 2, 3, 3, 2, 2, 3, 3))










# --------------------------------------------------------
# What if we change the bias value??????????
# --------------------------------------------------------------
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

#######
# Plots to put in the paper:
EUII_CRON <- ggplot(df_plot_euii, aes(x = Bias, y = EUII, color = Design, linetype = Design, shape = Design)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.2) +
  scale_color_manual(values = c("Fixed Uninformative (N=80)" = "#E69F00", 
                                "Fixed Informative (N=60)" = "#D55E00",
                                "Sequential Uninformative (Max N=80)" = "#56B4E9", 
                                "Sequential Informative (Max N=60)" = "#0072B2")) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "darkgray", linewidth = 1) +
  theme_bw(base_size = 14) +
  labs(
    title = "EUII Sensitivity to Drift",
    subtitle = expression("True Effect Size " * Delta * " = 50"),
    x = "Drift (True Control Mean - Prior Mean)",
    y = "EUII"
  ) +
  theme(legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank())






# Reshape T1E for ggplot
df_t1e_long <- df_t1e_total %>%
  pivot_longer(cols = starts_with("T1E"), names_to = "Design", values_to = "Alpha")

T1E_CRON <- ggplot(df_t1e_long, aes(x = Bias, y = Alpha, color = Design)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.2) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red") + # Standard 5% line
  theme_bw(base_size = 14) +
  labs(  x = "Drift", y = "Type-I Error Rate")




# Filter Power for Delta = 50 and reshape
df_power_long <- df_power_total %>% 
  filter(Delta == 50) %>%
  pivot_longer(cols = starts_with("Power"), names_to = "Design", values_to = "Power_Val")

Power_plot_CRON <- ggplot(df_power_long, aes(x = Bias, y = Power_Val, color = Design)) +
  geom_line(linewidth = 1.2) + geom_point(size = 2.2) +
  theme_bw(base_size = 14) +
  labs(x = "Drift", y = "Statistical Power")







