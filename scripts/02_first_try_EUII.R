library(dplyr)
library(gsbDesign)

# First application fo EUII to bayesian criterion 
design <- gsbDesign(nr.stages = 2,
                    patients = c(10,20),
                    sigma=88,
                    criteria.success = c(0,0.95,50,0.5),
                    criteria.futility = c(40,0.9),
                    prior.control = c(49,20),
                    prior.treatment = c(49,0))


simulation <- gsbSimulation(
  truth=list(49,c(0,40,50,60,70)),
  grid.type="sliced",
  method = "both",
  type.update="per arm",
  nr.sim=100000, #note this is exaclty the code of the paper, 1e6, not 1e5 simulations.
  warnings.sensitivity = 500,
  seed=1)

#JUST a small NOTE --> Ming shows examples of 10k simulations, then all of a sudden uses the results from the paper that have 100k
#any way now I will be only using 100k

result_new <- gsb(design=design,simulation=simulation)

table_total <- tab(result_new, what = "cumulative all") #use always cumulative all, then also it contains I guess Power at T1E
# cause it basically gives you P (reject | null) and P(reject | H1) that is stage2.suc in this case. 

table_n <- tab(result_new, what = "sample size")
table_n

n_matrix <- table_n[, "stage2", drop = FALSE]
T1E <- table_total[1, "stage2.suc"] # so we take line one where delta = 0 ( i.e ground truth is H0) and then take cumulative success. So P(reject | H0)
T1E

powers <- table_total[2:5, c("delta", "stage2.suc")] # here Power = P(reject |H1)



#compute  LR+ and LR-
LR_pos <- cbind(powers[,"delta"], powers[, "stage2.suc"] / T1E)
LR_neg <- cbind(powers[,"delta"], (1 - powers[, "stage2.suc"]) / (1- T1E) ) 



#Plot them as Ming did
xlab = "Delta", ylab = "Positive Likelihood Ratio (LR+)")
axis(1, at=c(40, 50, 60, 70))
text(LR_pos[, 1] - 1.3, LR_pos[, 2] + 3, 
     labels = round(LR_pos[, 2], 2), pos = 4, cex = 0.8, col = "red")



plot(LR_neg[, 1], LR_neg[, 2], type="o", col="blue", lwd=2, xaxt="n", ylim = c(0,1),
     xlab = "Delta", ylab = "Negative Likelihood Ratio (LR-)")
axis(1, at=c(40, 50, 60, 70))
text(LR_neg[, 1] - 1, LR_neg[, 2] + 0.05, 
     labels = round(LR_neg[, 2], 2), pos = 4, cex = 0.8, col = "red")


#Now compute DOR
DOR <- cbind(delta = LR_pos[, 1], DOR = LR_pos[, 2] / LR_neg[, 2])

# Now we need to compute the EUII, for that we need E(N+) and E(N-)

#cannot use n_matrix, since it contains E(N_{H0}) in the first line and E(N_{H1}) for the others line
N1 <- 30  #first sample is (20 + 10 people) 
N2 <- 60  # second sample is (20 + 10 people)

#foundamental question to inspect: SHOULD WE USE THE SAMPLE SIZE OF THE PRIOR???? Maybe yes

#for now we avoid this , but we will inspect

# Extract exact probabilities from the 'cumulative all' table
#here, success is "+"

P_S1 <- table_total[, "stage1.suc"] # P(+ at stage 1)
P_S_cum <- table_total[, "stage2.suc"] # P(+)
P_S2 <- P_S_cum - P_S1  #  P(+ at stage 2)



#Now we put on the same bucket all those that are not +.

P_F1 <- table_total[, "stage1.fut"] # P (fut at 1)
P_notS_cum <- 1 - P_S_cum # P(that we don't get + in general) = P (fut at 1 )+ P(fut at 2) + P(indef at 2)
P_notS2 <- P_notS_cum - P_F1 # Probability of stopping at stage 2 without success = P(futility at 2) + P(indefined at 2)





# Compute conditional expected sample sizes
# E (N | +, H0) = n(stop at 1) * P (stop at 1| +, H0) + n(stop at 2) * P(stop at 2 | +, H0)

# P ( stop at i| +, H0) = P (success at i, H0) / P (success at 1 or 2, H0) and same for stop at 2
P_S1_sign <- P_S1 / P_S_cum
P_S2_sign <- P_S2 / P_S_cum

E_N_cond_plus <- N1 * P_S1_sign + N2 * P_S2_sign 
E_N_cond_minus <- (N1 * P_F1 + N2 * P_notS2) / P_notS_cum


# NOTE, for now, all those that are not significant findings (so are different from +) go in the -. 
#We should do a sensitivity analysis to understand wheter to use only the rejcection for the minus or also the ind 



# Split into H0 (row 1) and H1 (rows 2:5)
E_N0_plus <- E_N_cond_plus[1]   # E[N|+, H0] = E[N_{+, H0}]
E_N0_minus <- E_N_cond_minus[1] # E[N|+, H1]

E_N1_plus <- E_N_cond_plus[2:5] # E[N|-, H0]
E_N1_minus <- E_N_cond_minus[2:5] #  E[N|-, H1]

# Define prior
Pr_H1 <- 0.1
Pr_H0 <- 1 - Pr_H1

#Define Odds
O_H1 <- Pr_H1 / (1 - Pr_H1)

O_H1_given_sig    <- LR_pos[, 2] * O_H1
O_H1_given_nonsig <- LR_neg[, 2] * O_H1

#do back trasform
Pr_H1_given_sig <- O_H1_given_sig / (O_H1_given_sig+ 1)
Pr_H0_given_sig <- 1- Pr_H1_given_sig

Pr_H1_given_nonsig <- O_H1_given_nonsig / (O_H1_given_nonsig+ 1)
Pr_H0_given_nonsig <- 1- Pr_H1_given_nonsig


# Compute final E(N+) and E(N-) (Equations 16 and 17 from the paper)
# E (N | +) = E [N | +,0] * P (H0|+) + E [N | + , H1] * P (H1 | +)
E_N_plus <- E_N0_plus * Pr_H0_given_sig + E_N1_plus * Pr_H1_given_sig
E_N_minus <- E_N0_minus * Pr_H0_given_nonsig + E_N1_minus * Pr_H1_given_nonsig

# 8. Compute the First-Order EUII
EUII <- (LR_pos[, 2] ^ (1 / E_N_plus)) / (LR_neg[, 2] ^ (1 / E_N_minus))

# View everything together
final_table <- cbind(delta = LR_pos[, 1], 
                     DOR = DOR[, 2], 
                     E_N_plus = E_N_plus, 
                     E_N_minus = E_N_minus, 
                     EUII = EUII)
print(round(final_table, 3))




# Compute conditional second moments: E(N^2 | +) and E(N^2 | -)
E_N2_cond_plus <- (N1^2 * P_S1 + N2^2 * P_S2) / P_S_cum
E_N2_cond_minus <- (N1^2 * P_F1 + N2^2 * P_notS2) / P_notS_cum

# Variance (N | + , H_0) = E[N^2 | H_0, +] - (E[N | +, H_0])^2
Var_N_cond_plus <- E_N2_cond_plus - E_N_cond_plus^2
Var_N_cond_minus <- E_N2_cond_minus - E_N_cond_minus^2

# Split the variances into H0 (row 1) and H1 (rows 2:5)
Var_N0_plus <- Var_N_cond_plus[1] #Var (N | + , H_0)
Var_N0_minus <- Var_N_cond_minus[1] #Variance (N | - , H_0)

Var_N1_plus <- Var_N_cond_plus[2:5] 
Var_N1_minus <- Var_N_cond_minus[2:5]

#  Compute overall variance using the Law of Total Variance (Equation 18)
# Var(N) = E[Var(N|H)] + Var(E[N|H]) and then expand
Var_N_plus <- (Var_N0_plus * Pr_H0_given_sig + Var_N1_plus * Pr_H1_given_sig) +
  ((E_N0_plus - E_N_plus)^2 * Pr_H0_given_sig + 
     (E_N1_plus - E_N_plus)^2 * Pr_H1_given_sig)

Var_N_minus <- (Var_N0_minus * Pr_H0_given_nonsig + Var_N1_minus * Pr_H1_given_nonsig) +
  ((E_N0_minus - E_N_minus)^2 * Pr_H0_given_nonsig + 
     (E_N1_minus - E_N_minus)^2 * Pr_H1_given_nonsig)

# Compute the Coefficient of Variation (CV)
CV_N_plus <- sqrt(Var_N_plus) / E_N_plus
CV_N_minus <- sqrt(Var_N_minus) / E_N_minus

#  Calculate the adjusted exponents for the second-order formula
adj_exp_plus <- (1 + CV_N_plus^2) / E_N_plus
adj_exp_minus <- (1 + CV_N_minus^2) / E_N_minus

# Compute the Second-Order EUII (Equation 15)
EUII_2nd_order <- (LR_pos[, 2] ^ adj_exp_plus) / (LR_neg[, 2] ^ adj_exp_minus)

# Compare First-Order and Second-Order EUII
comparison_table <- cbind(delta = LR_pos[, 1],
                          EUII_1st = EUII,
                          EUII_2nd = EUII_2nd_order)
print(round(comparison_table, 4))





















#FUNCTION
compute_euii <- function(table_total, N1 = 30, N2 = 60, Pr_H1 = 0.1, strict_futility = FALSE) {
  
  delta <- table_total[, "delta"]
  
  # Extract exact probabilities from the 'cumulative all' table
  #here, success is "+"
  P_S1 <- table_total[, "stage1.suc"] 
  P_S_cum <- table_total[, "stage2.suc"]
  P_S2 <- P_S_cum - P_S1
  
  T1E <- P_S_cum[1]
  Power <- P_S_cum[-1]
  LR_pos <- Power / T1E
  
  P_F1 <- table_total[, "stage1.fut"]
  
  # IDEA: do we want to include in the - cases only futility cases or both futility and indifferent at last step? 
  if (strict_futility) { #include only futility cases 
    P_neg_cum <- table_total[, "stage2.fut"]
    P_neg2    <- P_neg_cum - P_F1
    
    Fut0 <- P_neg_cum[1]
    Fut1 <- P_neg_cum[-1]
    LR_neg <- Fut1 / Fut0 # P(- | H1) / P (-|H0)
    
  } else {
    #Now we put on the same bucket all those that are not +.
    # I.e. we get P (not significant result) = P (futility at 1) + P(futility at 2) + P(indefined at 2)
    # Note that  P (not significant result) = 1- P(success (i.e. significante result))
    P_neg_cum <- 1 - P_S_cum
    P_neg2    <- P_neg_cum - P_F1
    
    LR_neg <- (1 - Power) / (1 - T1E)
  }
  
  DOR <- LR_pos / LR_neg
  
  # Compute conditional expected sample sizes
  # E (N | +, H0) = n(stop at 1) * P (stop at 1| +, H0) + n(stop at 2) * P(stop at 2 | +, H0)
  # where P ( stop at i| +, H0) = P (success at i, H0) / P (success at 1 or 2, H0) and same for stop at 2
  
  E_N_cond_plus  <- (N1 * P_S1 + N2 * P_S2) / P_S_cum
  E_N_cond_minus <- (N1 * P_F1 + N2 * P_neg2) / P_neg_cum
  
  # Split into H0 (row 1) and H1 (rows 2:5)
  E_N0_plus  <- E_N_cond_plus[1]   # E[N|+, H0] = E[N_{+, H0}]
  E_N0_minus <- E_N_cond_minus[1]  # E[N|-, H0]
  
  E_N1_plus  <- E_N_cond_plus[-1]  # E[N|+, H1]
  E_N1_minus <- E_N_cond_minus[-1] # E[N|-, H1]
  
  # Define prior
  Pr_H0 <- 1 - Pr_H1
  
  #Define Odds
  O_H1 <- Pr_H1 / (1 - Pr_H1)
  
  O_H1_given_sig    <- LR_pos * O_H1
  O_H1_given_nonsig <- LR_neg * O_H1
  
  #do back trasform
  Pr_H1_given_sig <- O_H1_given_sig / (O_H1_given_sig + 1)
  Pr_H0_given_sig <- 1 - Pr_H1_given_sig
  
  Pr_H1_given_nonsig <- O_H1_given_nonsig / (O_H1_given_nonsig + 1)
  Pr_H0_given_nonsig <- 1 - Pr_H1_given_nonsig
  
  # Compute final E(N+) and E(N-) (Equations 16 and 17 from the paper)
  # E (N | +) = E [N | +,0] * P (H0|+) + E [N | + , H1] * P (H1 | +)
  E_N_plus  <- E_N0_plus * Pr_H0_given_sig + E_N1_plus * Pr_H1_given_sig
  E_N_minus <- E_N0_minus * Pr_H0_given_nonsig + E_N1_minus * Pr_H1_given_nonsig
  
  #  Compute EUII
  EUII <- (LR_pos ^ (1 / E_N_plus)) / (LR_neg ^ (1 / E_N_minus))
  
  result_df <- data.frame(
    Delta = delta[-1],
    Strict_Futility = strict_futility,
    DOR = round(DOR, 2),
    E_N_plus = round(E_N_plus, 2),
    E_N_minus = round(E_N_minus, 2),
    EUII = round(EUII, 4)
  )
  
  return(result_df)
}


# Run standard
res_standard <- compute_euii(table_total, N1 = 30, N2 = 60, Pr_H1 = 0.1, strict_futility = FALSE)

# Run strict
res_strict <- compute_euii(table_total, N1 = 30, N2 = 60, Pr_H1 = 0.1, strict_futility = TRUE)

# Combine and look at the differences
comparison <- left_join(res_standard, res_strict, by = "Delta")

print(comparison)









#####################################


library(dplyr)
library(tidyr)
library(ggplot2)

# 1. Define the priors we want to test (matching the paper's choices)
priors <- c(0.01, 0.1, 0.5)

# 2. Loop through each prior and calculate the EUII
# lapply runs our function for each value in 'priors', and bind_rows stitches the dataframes together
all_results <- bind_rows(lapply(priors, function(p) {
  
  # Calculate using the function we built in the last step
  res <- calculate_euii(table_total = table_total, N1 = 30, N2 = 60, Pr_H1 = p)
  
  # Add a column to track which prior generated these specific rows
  res$Pr_H1 <- p
  
  return(res)
}))

# 3. Create a clean factor label so the facet headers look nice
all_results <- all_results %>%
  mutate(Prior_Label = factor(Pr_H1, 
                              levels = c(0.01, 0.1, 0.5),
                              labels = c("Pr(H1) = 1%", "Pr(H1) = 10%", "Pr(H1) = 50%")))

# 4. Reshape the data for ggplot2 (combining 1st and 2nd order into one column)
plot_data_multi <- pivot_longer(all_results, 
                                cols = c("EUII_1st_Order", "EUII_2nd_Order"),
                                names_to = "Approximation",
                                values_to = "EUII")

# Clean up the legend labels
plot_data_multi$Approximation <- ifelse(plot_data_multi$Approximation == "EUII_1st_Order", 
                                        "1st Order", 
                                        "2nd Order")

# Assuming 'plot_data_multi' is already generated from the previous steps

# Generate the single combined plot
single_combined_plot <- ggplot(plot_data_multi, aes(x = Delta, y = EUII, 
                                                    color = Prior_Label, 
                                                    linetype = Approximation)) +
  geom_line(linewidth = 1) +
  # Adding shapes mapped to the prior for extra clarity alongside the colors
  geom_point(aes(shape = Prior_Label), size = 3) + 
  theme_bw() +
  labs(title = "Experimental Unit Information Index (EUII)",
       subtitle = "1st vs. 2nd Order Approximations across Prior Probabilities",
       x = expression("Treatment Effect (" * Delta * ")"),
       y = "EUII",
       color = "Prior Probability",
       shape = "Prior Probability",
       linetype = "Approximation") +
  # Map specific linetypes: solid for 1st order, dashed for 2nd order
  scale_linetype_manual(values = c("1st Order" = "solid", "2nd Order" = "dashed")) +
  # Use a nice, distinct color palette for the three priors
  scale_color_manual(values = c("Pr(H1) = 1%" = "#E69F00", 
                                "Pr(H1) = 10%" = "#56B4E9", 
                                "Pr(H1) = 50%" = "#009E73")) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.key.width = unit(2, "cm") # Make the lines wider in the legend so the dashes are clear
  )

# Display the final plot
print(single_combined_plot)

