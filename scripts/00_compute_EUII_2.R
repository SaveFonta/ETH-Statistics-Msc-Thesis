library(dplyr)
library(tidyr)
library(ggplot2)
library(gsbDesign)






compute_euii <- function(table_total,
                         N         = c(30, 60), 
                         prior_N          = 0,
                         Pr_H1            = 0.1,
                         first_order      = FALSE,
                         strict_futility  = FALSE,
                         exact            = TRUE,
                         separate         = FALSE,
                         second_order     = FALSE,
                         null_value       = 0  ) {
  
  #Never gonna use prior_N, but I want to keep it in case I want to do the comparison with the bayesian design, where I can use it to add the prior sample size to the physical sample size.
  #N <- N + prior_N
  
  K     <- length(N) #how many stops? 
  delta <- table_total[, "delta"]
  
  # Extract exact probabilities from the 'cumulative all' table for K stages
  suc_cols <- paste0("stage", 1:K, ".suc")
  fut_cols <- paste0("stage", 1:K, ".fut")
  
  #note that if we use method = "numerical integration", the col names also have .cum
  if (!all(suc_cols %in% colnames(table_total))) {
    suc_cols <- paste0("stage", 1:K, ".cum.suc")
    fut_cols <- paste0("stage", 1:K, ".cum.fut")
  }
  if (!all(suc_cols %in% colnames(table_total))) {
      stop("Could not find either '.suc' or '.cum.suc' columns in the provided table.\nAre you sure it is the output table from tab(..., what = 'cumulative all')?")
    }



  
  P_S_cum_mat <- as.matrix(table_total[, suc_cols, drop = FALSE])  # cumulative success at each stage 
  P_F_cum_mat <- as.matrix(table_total[, fut_cols, drop = FALSE])  # cumulative futility at each stage
  
  # We need the probs of stopping *exactly* at stage k!! Since table_total gives the cumulative
  P_S_exact <- cbind(P_S_cum_mat[, 1],
                   if (K > 1) P_S_cum_mat[, 2:K, drop = FALSE] - P_S_cum_mat[, 1:(K-1), drop = FALSE])
  #what I am doing is:
  #                                          1             2            3                        K
  # P(stopp for a success exactly at i) = col 1   col 2 - col 1  col 3 - col2      ...    col K - col K-1
  #

  #so we have P (success at k) = P(stopp k and +)
  
  #We need also the same for futility
  P_F_exact <- cbind(P_F_cum_mat[, 1],
                   if (K > 1) P_F_cum_mat[, 2:K, drop = FALSE] - P_F_cum_mat[, 1:(K-1), drop = FALSE])
  
  P_S_cum <- P_S_cum_mat[, K]   #  cumulative success probability at last stage 

  # find row where delta is the null, usuall is is zero but for non inferiority can be different
  h0_index <- which(round(delta, 5) == round(null_value, 5))
  
  if (length(h0_index) == 0) stop(sprintf("The null_value %f was not found in the 'delta' grid of table_total. Check your truth_grid.", null_value))
  

  T1E    <- P_S_cum[h0_index]
  Power  <- P_S_cum[-h0_index]
  LR_pos <- Power / T1E
  
  P_F1 <- P_F_exact[, 1]  #prob of futility at first step
  
  if (strict_futility) { #only futility counts....
    P_neg_cum <- P_F_cum_mat[, K] # Cumulative P (futility) = P (- | ...)
    P_neg_exact <- P_F_exact
    
    Fut0   <- P_neg_cum[h0_index]
    Fut1   <- P_neg_cum[-h0_index]
    LR_neg <- Fut1 / Fut0  # P(- | H1) / P(- | H0)
    
  } else {
    # Now we put in the same bucket all those that are not + (both ind at last stage and fut).
    P_neg_cum <- 1 - P_S_cum #here we have first element (1-alpha) and the others are (Power)
    
    # 
    P_neg_exact      <- P_F_exact #need to use it P(stop at k | - , H-)
    P_neg_exact[, K] <- P_neg_cum - rowSums(P_F_exact[, -K, drop = FALSE]) #--> we are computing P (stop at K | -) = P(total Negative) - sum(from k=1 to K-1) P(stop early for futility at stage K)
    #basically we have the total prob of not success (-), and we also have prob of failure at each step, so the mass remaining after we accounted for failure at each step except for the last, is the 
    #mass of failure in the end.
    #We could have done the same also by using P_neg_exact annd at the end adding it the P_ind[, K], but like that we dont use it, I like this way more  

    LR_neg <- (1 - Power) / (1 - T1E) #easy
  }
  
  DOR <- LR_pos / LR_neg
  
  
  # E(N+) and E(N-) (Equations 16 and 17 from the paper)
  # E[N | +, Hx] = Sum(N_k * P(Stop at stage k AND +)) / P(+) = Sum(N_k * P(Stop at stage k | +)
  E_N_cond_plus  <- drop(P_S_exact  %*% N) / P_S_cum
  E_N_cond_minus <- drop(P_neg_exact %*% N) / P_neg_cum
  

  # Split into H0 (row 1) and H1 (rows 2:end)
  E_N0_plus  <- E_N_cond_plus[h0_index]    # E[N|+, H0]
  E_N0_minus <- E_N_cond_minus[h0_index]   # E[N|-, H0]
  E_N1_plus  <- E_N_cond_plus[-h0_index]   # E[N|+, H1]
  E_N1_minus <- E_N_cond_minus[-h0_index]  # E[N|-, H1]
  
  
  # Define prior and Odds
  Pr_H0 <- 1 - Pr_H1
  O_H1  <- Pr_H1 / (1 - Pr_H1)
  
  O_H1_given_sig    <- LR_pos * O_H1
  O_H1_given_nonsig <- LR_neg * O_H1
  
  # Do back transform to get posteriors
  Pr_H1_given_sig <- O_H1_given_sig / (O_H1_given_sig + 1)
  Pr_H0_given_sig <- 1 - Pr_H1_given_sig
  
  Pr_H1_given_nonsig <- O_H1_given_nonsig / (O_H1_given_nonsig + 1)
  Pr_H0_given_nonsig <- 1 - Pr_H1_given_nonsig
  
  
  # E[N|+] = E[N|H_0, +] P[H_0 | +] + E[N | H_1, +] P[H_1 | +]
  E_N_plus  <- E_N0_plus  * Pr_H0_given_sig    + E_N1_plus  * Pr_H1_given_sig
  E_N_minus <- E_N0_minus * Pr_H0_given_nonsig + E_N1_minus * Pr_H1_given_nonsig
  
  
  # First-Order EUII
  EUII_1st_Order <- (LR_pos ^ (1 / E_N_plus)) / (LR_neg ^ (1 / E_N_minus))
  
  
  # -------------------------------------------------------------------
  # E[1/N | +, Hx] = sum_k (1/N_k) * P(stop at k and + | Hx) / P(+)
  # -------------------------------------------------------------------
  E_invN_cond_plus  <- drop(P_S_exact  %*% (1 / N)) / P_S_cum
  E_invN_cond_minus <- drop(P_neg_exact %*% (1 / N)) / P_neg_cum
  
  E_invN0_plus  <- E_invN_cond_plus[h0_index]    # E[1/N|+, H0]
  E_invN0_minus <- E_invN_cond_minus[h0_index]   # E[1/N|-, H0]
  E_invN1_plus  <- E_invN_cond_plus[-h0_index]   # E[1/N|+, H1]
  E_invN1_minus <- E_invN_cond_minus[-h0_index]  # E[1/N|-, H1]
  
  E_invN_plus  <- E_invN0_plus  * Pr_H0_given_sig    + E_invN1_plus  * Pr_H1_given_sig
  E_invN_minus <- E_invN0_minus * Pr_H0_given_nonsig + E_invN1_minus * Pr_H1_given_nonsig
  
  
  # Exact EUII
  EUII_Exact <- (LR_pos ^ E_invN_plus) / (LR_neg ^ E_invN_minus)
  
  # Another adding: we can remove the choice of how much I believe in P(H_0), just need to
  # create two different EUII
  EUII_H0 <- (LR_pos ^ (1 / E_N0_plus)) / (LR_neg ^ (1 / E_N0_minus))
  EUII_H1 <- (LR_pos ^ (1 / E_N1_plus)) / (LR_neg ^ (1 / E_N1_minus))
  
  
  # Second order computations:
  E_N_plus_squared  <- as.numeric(drop(P_S_exact  %*% N^2)) / P_S_cum
  E_N_minus_squared <- as.numeric(drop(P_neg_exact %*% N^2)) / P_neg_cum
  
  E_N0_plus_squared  <- E_N_plus_squared[h0_index]
  E_N1_plus_squared  <- E_N_plus_squared[-h0_index]
  
  E_N0_minus_squared <- E_N_minus_squared[h0_index]
  E_N1_minus_squared <- E_N_minus_squared[-h0_index]
  
  Var_N0_plus  <- pmax(0, E_N0_plus_squared  - E_N0_plus^2)
  Var_N1_plus  <- pmax(0, E_N1_plus_squared  - E_N1_plus^2)
  Var_N0_minus <- pmax(0, E_N0_minus_squared - E_N0_minus^2)
  Var_N1_minus <- pmax(0, E_N1_minus_squared - E_N1_minus^2)
  
  # Law of total probability for variances (Equation 18)
  Var_N_plus  <- Var_N0_plus  * Pr_H0_given_sig    + Var_N1_plus  * Pr_H1_given_sig    +
    (E_N0_plus  - E_N_plus)^2  * Pr_H0_given_sig    +
    (E_N1_plus  - E_N_plus)^2  * Pr_H1_given_sig
  
  Var_N_minus <- Var_N0_minus * Pr_H0_given_nonsig + Var_N1_minus * Pr_H1_given_nonsig +
    (E_N0_minus - E_N_minus)^2 * Pr_H0_given_nonsig +
    (E_N1_minus - E_N_minus)^2 * Pr_H1_given_nonsig

  Var_N_plus  <- pmax(0, Var_N_plus)
  Var_N_minus <- pmax(0, Var_N_minus)
  
  # Coefficient of variation
  CV_N_plus  <- sqrt(Var_N_plus)  / E_N_plus
  CV_N_minus <- sqrt(Var_N_minus) / E_N_minus
  
  # Second-order EUII approximation (Equation 15)
  EUII_2nd_Order <- (LR_pos ^ ((1 + CV_N_plus^2)  / E_N_plus)) /
    (LR_neg ^ ((1 + CV_N_minus^2) / E_N_minus))
  
  
  result_df <- data.frame(Delta = delta[-h0_index], DOR = round(DOR, 2))
  
  if (first_order)
    result_df$EUII_1st_Order <- round(EUII_1st_Order, 4)
  
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

  N <- data.frame(Delta = delta[-h0_index], E_N_plus = E_N_plus, E_N_minus = E_N_minus, E_invN_plus = E_invN_plus,  E_invN_minus = E_invN_minus)
  effective_N <- data.frame(
    Delta = delta[-h0_index],
    effective_N_plus = 1 / E_invN_plus,
    effective_N_minus = 1 / E_invN_minus,
    CV_N_plus = CV_N_plus,
    CV_N_minus = CV_N_minus
  )
  DOR_df <- data.frame(Delta = delta[-h0_index], DOR = DOR)
  LR_df <- data.frame(Delta = delta[-h0_index], LR_pos = LR_pos, LR_neg = LR_neg)
  P_H1_posterior <- data.frame(Delta = delta[-h0_index], Pr_H1_given_sig = Pr_H1_given_sig, Pr_H1_given_nonsig = Pr_H1_given_nonsig) 


  results <- list(result_df = result_df, DOR_df = DOR_df, N = N, effective_N = effective_N, LR_df = LR_df, P_H1_posterior = P_H1_posterior)
  
  return(results)
}


##################

# if we have a big sample size of the prior, the bayesian method will always have a smaller EUII
# This makes not much sense! If we want to compare frequentist and bayesian group sequential this would mean
# that we always would choose frequentist due to a higher EUII


plot_euii_comparison <- function(table_total,
                                 N         = c(30, 60),  
                                 prior_N          = 0,
                                 priors           = c(0.01, 0.1, 0.5),
                                 first_order      = FALSE,
                                 strict_futility  = FALSE,
                                 exact            = FALSE,
                                 separate         = FALSE,
                                 second_order     = FALSE) {
  
  # take into account prior_N
  N <- N + prior_N
  
  # compute EUII for each prior
  all_results <- bind_rows(lapply(priors, function(p) {
    res <- compute_euii(table_total     = table_total,
                        N        = N,
                        prior_N         = 0,          # already absorbed above
                        Pr_H1           = p,
                        first_order     = first_order,
                        strict_futility = strict_futility,
                        exact           = exact,
                        separate        = separate,
                        second_order    = second_order)
    res <- res$result_df                 
    res$Pr_H1 <- p  # create new col
    return(res)
  }))


  # Create labels for prior
  prior_labels <- paste0("Pr(H1) = ", priors * 100, "%")
  all_results  <- all_results %>%
    mutate(Prior_Label = factor(Pr_H1, levels = priors, labels = prior_labels))  # create fact
  
  cols_to_pivot <- c()
  
  if (first_order)  cols_to_pivot <- c("EUII_1st_Order")
  if (second_order) cols_to_pivot <- c(cols_to_pivot, "EUII_2nd_Order")
  if (exact)        cols_to_pivot <- c(cols_to_pivot, "EUII_Exact")
  if (separate)     cols_to_pivot <- c(cols_to_pivot, "EUII_H0", "EUII_H1")
  
  # pivot data for ggplot
  plot_data_multi <- pivot_longer(all_results,
                                  cols      = all_of(cols_to_pivot),
                                  names_to  = "Approximation",
                                  values_to = "EUII")
  
  # Clean the legend labels (e.g., "EUII_1st_Order" -> "1st Order")
  plot_data_multi$Approximation <- gsub("EUII_", "", plot_data_multi$Approximation)
  plot_data_multi$Approximation <- gsub("_",     " ", plot_data_multi$Approximation)
  
  p <- ggplot(plot_data_multi, aes(x = Delta, y = EUII,
                                   color    = Approximation,
                                   linetype = Prior_Label)) +
    geom_line(linewidth = 1) +
    geom_point(aes(shape = Prior_Label), size = 3) +
    theme_bw() +
    labs(x        = "Treatment Effect",
         y        = "EUII",
         color    = "Approximation",
         shape    = "Prior Probability",
         linetype = "Prior Probability") +
    theme(
      plot.title       = element_text(face = "bold", size = 14),
      axis.text        = element_text(size = 10),
      axis.title       = element_text(size = 12),
      legend.position  = "right",
      legend.key.width = unit(2, "cm")
    )
  
  return(p)
}


# ---- Example calls

library(dplyr)
library(gsbDesign)

# First application fo EUII to bayesian criterion 
design <- gsbDesign(nr.stages = 5,
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
result_new <- gsb(design=design,simulation=simulation)

table_total <- tab(result_new, what = "cumulative all") #use always cumulative all, then also it contains I guess Power at T1E
# cause it basically gives you P (reject | null) and P(reject | H1) that is stage2.suc in this case. 



res <- compute_euii(table_total, N = c(30, 60, 90, 120, 150), prior_N = 200, Pr_H1 = 0.1,
                     first_order = TRUE, second_order = FALSE, exact = TRUE)

# 3-stage (table_total must contain stage3.suc / stage3.fut columns)
# res <- compute_euii(table_total, N = c(30, 60, 90), Pr_H1 = 0.1,
#                     first_order = TRUE, exact = TRUE)

# plot_euii_comparison(table_total, N = c(30, 60), first_order = TRUE, second_order = TRUE)







