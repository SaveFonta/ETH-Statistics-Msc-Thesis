###########################################
# Function to compute EUII for 2 stages (need to extend it)
##########################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(gsbDesign)


compute_euii <- function(table_total, N1 = 30, N2 = 60,  prior_N = 0, Pr_H1 = 0.1, first_order= TRUE, strict_futility = FALSE, exact = FALSE, separate = FALSE, second_order= FALSE) {
  
  #not really sure about this. We just increase our sample size 
  N1 <- N1 + prior_N
  N2 <- N2 + prior_N
  
  delta <- table_total[, "delta"]
  
  # Extract exact probabilities from the 'cumulative all' table
  # here, success is "+"
  P_S1 <- table_total[, "stage1.suc"] 
  P_S_cum <- table_total[, "stage2.suc"]
  P_S2 <- P_S_cum - P_S1
  
  T1E <- P_S_cum[1]
  Power <- P_S_cum[-1]
  LR_pos <- Power / T1E
  
  P_F1 <- table_total[, "stage1.fut"]
  
  if (strict_futility) {
    P_neg_cum <- table_total[, "stage2.fut"]
    P_neg2    <- P_neg_cum - P_F1
    
    Fut0 <- P_neg_cum[1]
    Fut1 <- P_neg_cum[-1]
    LR_neg <- Fut1 / Fut0 # P(- | H1) / P (-|H0)
    
  } else {
    # Now we put on the same bucket all those that are not + (both ind at last stage and fut).
    P_neg_cum <- 1 - P_S_cum
    P_neg2    <- P_neg_cum - P_F1
    
    LR_neg <- (1 - Power) / (1 - T1E)
  }
  
  DOR <- LR_pos / LR_neg
  
  
  E_N_cond_plus  <- (N1 * P_S1 + N2 * P_S2) / P_S_cum
  E_N_cond_minus <- (N1 * P_F1 + N2 * P_neg2) / P_neg_cum
  
  
  
  # Split into H0 (row 1) and H1 (rows 2:5)
  E_N0_plus  <- E_N_cond_plus[1]   # E[N|+, H0] 
  E_N0_minus <- E_N_cond_minus[1]  # E[N|-, H0]
  E_N1_plus  <- E_N_cond_plus[-1]  # E[N|+, H1]
  E_N1_minus <- E_N_cond_minus[-1] # E[N|-, H1]
  
  
  # Define prior and Odds
  Pr_H0 <- 1 - Pr_H1
  O_H1 <- Pr_H1 / (1 - Pr_H1)
  
  O_H1_given_sig    <- LR_pos * O_H1
  O_H1_given_nonsig <- LR_neg * O_H1
  
  # Do back transform to get posteriors
  Pr_H1_given_sig <- O_H1_given_sig / (O_H1_given_sig + 1)
  Pr_H0_given_sig <- 1 - Pr_H1_given_sig
  
  Pr_H1_given_nonsig <- O_H1_given_nonsig / (O_H1_given_nonsig + 1)
  Pr_H0_given_nonsig <- 1 - Pr_H1_given_nonsig
  
  
  # E(N+) and E(N-) (Equations 16 and 17 from the paper [cite: 322])
  E_N_plus  <- E_N0_plus * Pr_H0_given_sig + E_N1_plus * Pr_H1_given_sig
  E_N_minus <- E_N0_minus * Pr_H0_given_nonsig + E_N1_minus * Pr_H1_given_nonsig
  
  
  # First-Order EUII 
  EUII_1st_Order <- (LR_pos ^ (1 / E_N_plus)) / (LR_neg ^ (1 / E_N_minus))
  
  
  # -------------------------------------------------------------------
  # E[1/N | +, H0] = (1/n1) * P(stop at 1 | +, H0) + (1/n2) * P(stop at 2 | +, H0)
  # -------------------------------------------------------------------
  
  E_invN_cond_plus  <- ((1/N1) * P_S1 + (1/N2) * P_S2) / P_S_cum
  E_invN_cond_minus <- ((1/N1) * P_F1 + (1/N2) * P_neg2) / P_neg_cum
  
  E_invN0_plus  <- E_invN_cond_plus[1]   # E[1/N|+, H0]
  E_invN0_minus <- E_invN_cond_minus[1]  # E[1/N|-, H0]
  E_invN1_plus  <- E_invN_cond_plus[-1]  # E[1/N|+, H1]
  E_invN1_minus <- E_invN_cond_minus[-1] # E[1/N|-, H1]
  
  E_invN_plus  <- E_invN0_plus * Pr_H0_given_sig + E_invN1_plus * Pr_H1_given_sig
  E_invN_minus <- E_invN0_minus * Pr_H0_given_nonsig + E_invN1_minus * Pr_H1_given_nonsig
  
  
  # Exact EUII 
  EUII_Exact <- (LR_pos ^ E_invN_plus) / (LR_neg ^ E_invN_minus)
  
  #Another adding: we can remove the choice of how much I believe in P(H_0), just need to 
  # create two different EUII
  
  EUII_H0 <- (LR_pos ^ (1 / E_N0_plus)) / (LR_neg ^ (1 / E_N0_minus))
  EUII_H1 <- (LR_pos ^ (1 / E_N1_plus)) / (LR_neg ^ (1 / E_N1_minus))
  
  # Second order computations:
  E_N_plus_squared <- ( N1^2 * P_S1 + N2^2 * P_S2 )  /  P_S_cum
  E_N_minus_squared <- ( N1^2 * P_F1 + N2^2 * P_neg2 )  /  P_neg_cum
  
  E_N0_plus_squared <- E_N_plus_squared[1]
  E_N1_plus_squared <-  E_N_plus_squared[-1]
  
  E_N0_minus_squared <- E_N_minus_squared[1]
  E_N1_minus_squared <-  E_N_minus_squared[-1]
  
  Var_N0_plus <- E_N0_plus_squared - E_N0_plus^2
  Var_N1_plus <- E_N1_plus_squared - E_N1_plus^2
  
  Var_N0_minus <- E_N0_minus_squared - E_N0_minus^2
  Var_N1_minus <- E_N1_minus_squared - E_N1_minus^2
  
  # Law of total probability for variances (Equation 18)
  Var_N_plus <- Var_N0_plus * Pr_H0_given_sig + Var_N1_plus * Pr_H1_given_sig +
    (E_N0_plus - E_N_plus)^2 * Pr_H0_given_sig +
    (E_N1_plus - E_N_plus)^2 * Pr_H1_given_sig
  
  Var_N_minus <- Var_N0_minus * Pr_H0_given_nonsig + Var_N1_minus * Pr_H1_given_nonsig +
    (E_N0_minus - E_N_minus)^2 * Pr_H0_given_nonsig +
    (E_N1_minus - E_N_minus)^2 * Pr_H1_given_nonsig 
  
  # Coefficient of variation
  CV_N_plus <- sqrt(Var_N_plus) / E_N_plus
  CV_N_minus <- sqrt(Var_N_minus) / E_N_minus
  
  # Second-order EUII approximation (Equation 15)
  EUII_2nd_Order <- (LR_pos ^ ((1 + CV_N_plus^2) / E_N_plus)) / 
    (LR_neg ^ ((1 + CV_N_minus^2) / E_N_minus))
  
  result_df <- data.frame(
    Delta = delta[-1],
    DOR = round(DOR, 2)
  )
  
  if (first_order){
    result_df$EUII_1st_Order = round(EUII_1st_Order, 4)
  }
  
  if (exact) {
    result_df$EUII_Exact <- round(EUII_Exact, 4)
    if (first_order) result_df$Diff_Due_To_Jensen <- round(EUII_Exact - EUII_1st_Order, 4)
  }
  
  if (separate) {
    result_df$EUII_H0 <- round(EUII_H0, 4)
    result_df$EUII_H1 <- round(EUII_H1, 4)
  }
  
  if (second_order) {
    result_df$EUII_2nd_Order <- round(EUII_2nd_Order, 4)
    if (exact) result_df$Diff_with_2nd <- round(EUII_Exact - EUII_2nd_Order, 4)
  }
  
  return(result_df)
}

# 
#Plot
#

plot_euii_comparison <- function(table_total, 
                                 N1 = 30, 
                                 N2 = 60, 
                                 prior_N = 0,
                                 priors = c(0.01, 0.1, 0.5),
                                 first_order = FALSE,
                                 strict_futility = FALSE, 
                                 exact = FALSE, 
                                 separate = FALSE, 
                                 second_order = FALSE) {
  
  # take into account prior_N
  N1 <- N1 + prior_N
  N2 <- N2 + prior_N
  
  # compute EUII for each prior
  all_results <- bind_rows(lapply(priors, function(p) {
    res <- compute_euii(table_total = table_total, 
                        N1 = N1, 
                        N2 = N2, 
                        Pr_H1 = p,
                        first_order = first_order,
                        strict_futility = strict_futility,
                        exact = exact,
                        separate = separate,
                        second_order = second_order)
    res$Pr_H1 <- p #create new col
    return(res)
  }))
  
  # Create labels for prior
  prior_labels <- paste0("Pr(H1) = ", priors * 100, "%")
  all_results <- all_results %>%
    mutate(Prior_Label = factor(Pr_H1, levels = priors, labels = prior_labels)) #create factor
  
  cols_to_pivot <- c()
  
  if (first_order) cols_to_pivot <- c("EUII_1st_Order")
  if (second_order) cols_to_pivot <- c(cols_to_pivot, "EUII_2nd_Order")
  if (exact) cols_to_pivot <- c(cols_to_pivot, "EUII_Exact")
  if (separate) cols_to_pivot <- c(cols_to_pivot, "EUII_H0", "EUII_H1")
  
  # pivot data for ggplot
  plot_data_multi <- pivot_longer(all_results, 
                                  cols = all_of(cols_to_pivot),
                                  names_to = "Approximation",
                                  values_to = "EUII")
  
  # Clean the legend labels  (e.g., "EUII_1st_Order" -> "1st Order")
  plot_data_multi$Approximation <- gsub("EUII_", "", plot_data_multi$Approximation)
  plot_data_multi$Approximation <- gsub("_", " ", plot_data_multi$Approximation)
  
  # 
  p <- ggplot(plot_data_multi, aes(x = Delta, y = EUII, 
                                   color = Approximation, 
                                   linetype = Prior_Label)) +
    geom_line(linewidth = 1) +
    geom_point(aes(shape = Prior_Label), size = 3) + 
    theme_bw() +
    labs(x = "Treatment Effect",
         y = "EUII",
         color = "Approximation",          
         shape = "Prior Probability",      
         linetype = "Prior Probability") + 
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12),
      legend.position = "right",
      legend.key.width = unit(2, "cm") 
    )
  
  return(p)
}



# Creating the design from Example 3
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

result <- gsb(design=design,simulation=simulation)

table_total <- tab(result, what = "cumulative all")



res <- compute_euii(table_total, N1 = 30, N2 = 60, prior_N = 0, Pr_H1 = 0.1, first_order = TRUE, second_order = FALSE, exact = FALSE)

plot_euii_comparison(table_total, first_order = TRUE, second_order= TRUE )

