# -------------------
# Design 1 -> informative sequential

# Creating the informative design
design_informative <- gsbDesign(nr.stages = 2,
                    patients = c(10,20),
                    sigma=mysigma,
                    criteria.success = mycrit.success,
                    criteria.futility = mycrit.fut,
                    prior.control = c(49,20),
                    prior.treatment = c(49,0))


simulation_informative <- gsbSimulation(
  truth=list(49,alternative_deltas),
  grid.type="sliced",
  method = "both",
  type.update="per arm",
  nr.sim=100000, #note this is exaclty the code of the paper, 1e6, not 1e5 simulations.
  warnings.sensitivity = 500,
  seed=1)


result_informative<- gsb(design=design_informative,simulation=simulation_informative)

table_seq_informative <- tab(result_informative, what = "cumulative all")




# -------------------
# Design 2 -> uninformative sequential
design_uninformative <- gsbDesign(nr.stages = 2,
                                patients = c(20,20),
                                sigma=mysigma,
                                criteria.success = mycrit.success,
                                criteria.futility = mycrit.fut,
                                prior.control = c(49,0),
                                prior.treatment = c(49,0))


simulation_uninformative <- gsbSimulation(
  truth=list(49,alternative_deltas),
  grid.type="sliced",
  method = "both",
  type.update="per arm",
  nr.sim=100000, #note this is exaclty the code of the paper, 1e6, not 1e5 simulations.
  warnings.sensitivity = 500,
  seed=1)


result_uninformative<- gsb(design=design_uninformative,simulation=simulation_uninformative)

table_seq_uninformative <- tab(result_uninformative, what = "cumulative all")



#res_seq_info <- compute_euii(table_seq_informative, N1 = 30, N2 = 60, prior_N = 0, Pr_H1 = 0.1, exact = TRUE)
#res_seq_uninfo <- compute_euii(table_seq_uninformative, N1 = 40, N2 =  80, prior_N = 0, Pr_H1 = 0.1, exact = TRUE)
#res_info_added_prior <- compute_euii(table_seq_informative, N1= 30, N2 =60, prior_N = 20, Pr_H1 = 0.01, exact = TRUE)
# As stated in the preliminary paper, if we add the 20 patients, the EUII drops, but 
#it is still larger than with an uninformative prior (without any virtual patients). 


#new compute_euii function:
res_seq_info <- compute_euii(table_seq_informative, N= c(30, 60), prior_N = 0, Pr_H1 = 0.1, exact = TRUE)
res_seq_uninfo <- compute_euii(table_seq_uninformative, N = c(40,80), prior_N = 0, Pr_H1 = 0.1, exact = TRUE)



# -------------------
# Design 3 -> uninformative fixed

design_fixed_uninformative <- gsbDesign(nr.stages = 1, 
                    patients = c(40,40), 
                    sigma=mysigma,
                    criteria.success = mycrit.success,
                    criteria.futility = mycrit.fut,
                    prior.control = c(49,0),
                    prior.treatment = c(49,0))
## simulation scenarios
simulation_fixed_uninformative <- gsbSimulation(truth=list(49,alternative_deltas), 
                            grid.type="sliced",
                            type.update="per arm", 
                            nr.sim=100000,
                            warnings.sensitivity = 500, 
                            method="both")
## operating characteristics
result_fixed_uninformative <- gsb(design=design_fixed_uninformative,simulation=simulation_fixed_uninformative)
table_fixed_uninformative <- tab(result_fixed_uninformative, what = "cumulative all")




# -------------------
# Design 4 -> informative  fixed

design_fixed_informative <- gsbDesign(nr.stages = 1, 
                                        patients = c(20,40), 
                                        sigma=mysigma,
                                        criteria.success = mycrit.success,
                                        criteria.futility = mycrit.fut,
                                        prior.control = c(49,20),
                                        prior.treatment = c(49,0))
## simulation scenarios
simulation_fixed_informative <- gsbSimulation(truth=list(49,alternative_deltas), 
                                                grid.type="sliced",
                                                type.update="per arm", 
                                                nr.sim=100000,
                                                warnings.sensitivity = 500, 
                                                method="both")
## operating characteristics
result_fixed_informative <- gsb(design=design_fixed_informative,simulation=simulation_fixed_informative)

table_fixed_informative <- tab(result_fixed_informative, what = "cumulative all")





#Extract the values for the fixed, one day you will replace with compute_euii finished

delta_vals <- alternative_deltas[-1]

# --- Fixed Uninformative (Physical N = 80) ---
T1E_fixed_uninfo <- table_fixed_uninformative[1, "stage1.suc"]
Power_fixed_uninfo <- table_fixed_uninformative[-1, "stage1.suc"]
LR_pos_uninfo <- Power_fixed_uninfo / T1E_fixed_uninfo
LR_neg_uninfo <- (1 - Power_fixed_uninfo) / (1 - T1E_fixed_uninfo)
DOR_fixed_uninfo <- LR_pos_uninfo / LR_neg_uninfo

EUII_fixed_uninfo <- DOR_fixed_uninfo ^ (1/80) 

# --- Fixed Informative (Physical N = 60) ---
T1E_fixed_info <- table_fixed_informative[1, "stage1.suc"]
Power_fixed_info <- table_fixed_informative[-1, "stage1.suc"]
LR_pos_info <- Power_fixed_info / T1E_fixed_info
LR_neg_info <- (1 - Power_fixed_info) / (1 - T1E_fixed_info)
DOR_fixed_info <- LR_pos_info / LR_neg_info

EUII_fixed_info <- DOR_fixed_info ^ (1/60) 



#NOw create the df to plot!
df_fixed_uninfo <- data.frame(
  Delta = delta_vals,
  EUII = EUII_fixed_uninfo,
  Design = "Fixed Uninformative (N=80)"
)

df_fixed_info <- data.frame(
  Delta = delta_vals,
  EUII = EUII_fixed_info,
  Design = "Fixed Informative (N=60)"
)


df_seq_info <- data.frame(
  Delta = res_seq_info$result_df$Delta,
  EUII = res_seq_info$result_df$EUII_Exact,
  Design = "Sequential Informative (Max N=60)"
)

df_seq_uninfo <- data.frame(
  Delta = res_seq_uninfo$result_df$Delta,
  EUII = res_seq_uninfo$result_df$EUII_Exact,
  Design = "Sequential Uninformative (Max N=80)"
)

df_merged <- rbind(df_fixed_uninfo, df_fixed_info, df_seq_info, df_seq_uninfo)

final_plot_CRON <- ggplot(df_merged, aes(x = Delta, y = EUII, 
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

info <- res_seq_info$N[, 2:3] 
uninfo <- res_seq_uninfo$N[,2:3]
names(info) <- c("E_N_+_info", "E_N_-_info")
names(uninfo) <- c("E_N_+_uninfo", "E_N_-_uninfo")
delta <- res_seq_uninfo$result_df$Delta
N_data_CRON <- cbind(delta, info, uninfo)

library(xtable)
colnames(N_data_CRON) <- c("$\\Delta$", "$E[N_+]^{\\text{info}}$", "$E[N_-]^{\\text{info}}$", "$E[N_+]^{\\text{uninfo}}$", "$E[N_-]^{\\text{uninfo}}$")
N_data_CRON <- xtable(N_data_CRON, digits = c(0, 0,2,2,2,2))

