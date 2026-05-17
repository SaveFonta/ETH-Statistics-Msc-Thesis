library(dplyr)
library(stringr)
library(purrr)



# -----------------------------------------------------
# FUNCTIONS for Monte Carlo Sequential evaluation
# -----------------------------------------------------
 
 
oc2_seq_mc.normMix <- function(theta_1, theta_2, prior_1, prior_2, 
                                n1_seq, n2_seq, decisions_list,   # recall that n1_seq is cumulative
                                sigma_1 = 88, sigma_2 = 88,       
                                n_sim = 100000, seed = NULL,
                                weight.track = FALSE,
                                type.weight.track = "uninformative.mixture", n_info_comps = 3)  {     # added these 2 params: so we can track the informative weight both if only one uninformative component (given by robustify) 
                                                                                     # or even if we use a mixture as non-informative (e.g. 100 mmixture to approx a t-distri) 
 
  if (!is.null(seed)) set.seed(seed)
 
  if (length(theta_1) == 1) theta_1 <- rep(theta_1, n_sim)
  if (length(theta_2) == 1) theta_2 <- rep(theta_2, n_sim)
 
  # check length 
  if (length(theta_1) != n_sim) stop("Number of theta_1 provided should be either one for classic OC or equal to 
                                      n_sim for avg OC")
  if (length(theta_2) != n_sim) stop("Number of theta_2 provided should be either one for classic OC or equal to 
                                      n_sim for avg OC")
 
  K <- length(n2_seq)
 
  if (length(n2_seq) != length(n1_seq)) stop("length of n1_seq and n2_seq must coincide")
  if (length(n2_seq) != length(decisions_list)) stop("length of n2_seq and decision_list must coincide")

 
  active_idx <- 1:n_sim  # Keep track of trials that haven't stopped
 
  # Initialize a lot of tracking vectors, all of length n_sim 
  # (each element can be seen as an evolving trial) 
  stop_eff      <- rep(FALSE, n_sim)
  stop_fut      <- rep(FALSE, n_sim)
  stage_stopped <- rep(K, n_sim)  # At which stage does that trial stopped? (Default to final stage)
 
  y_2_curr <- rep(0, n_sim)
  y_1_curr <- rep(0, n_sim)
  n_2_prev <- 0
  n_1_prev <- 0
 
  # Vector of probabilities of stopping or futility for each stage
  p_succ_vec <- rep(0, K)
  p_fut_vec  <- rep(0, K)

  w_inf_traj <- rep(NA, K) #vector of average non-informative weight 
 
  for (k in 1:K) {
    if (length(active_idx) == 0) break
 
    n_active <- length(active_idx)
 
    dn_2 <- n2_seq[k] - n_2_prev
    dn_1 <- n1_seq[k] - n_1_prev
 
    # Generate new data for the trials still active
    # to have a vector of length active_idx
    y_1_inc <- if (dn_1 > 0) rnorm(n_active, mean = theta_1[active_idx], sd = sigma_1 / sqrt(dn_1)) else 0
    y_2_inc <- if (dn_2 > 0) rnorm(n_active, mean = theta_2[active_idx], sd = sigma_2 / sqrt(dn_2)) else 0

    # Calculate cumulative means
    # [mean till stage (k-1) * sample size (k-1) + mean of stage k * sample size stage k] / total sample size till stage k
    y_2_curr[active_idx] <- (n_2_prev * y_2_curr[active_idx] + dn_2 * y_2_inc) / n2_seq[k]
    y_1_curr[active_idx] <- (n_1_prev * y_1_curr[active_idx] + dn_1 * y_1_inc) / n1_seq[k]


if (weight.track) { 
    #average y of all the active 
    avg_y_2_curr <- mean(y_2_curr[active_idx]) 
    
    # Update the mixture prior using the average simulated data
    post_avg <- suppressMessages(RBesT::postmix(prior_2, 
                                                m = avg_y_2_curr, 
                                                n = n2_seq[k], 
                                                se = sigma_2 / sqrt(n2_seq[k])))
    
    if (type.weight.track == "robustify")  {
    #  since the robust component is always the last row from robustify
    idx_robust <- nrow(post_avg)
    
    w_robust <- post_avg["w", idx_robust]
    w_inf_traj[k] <- 1 - w_robust }
    else if (type.weight.track == "uninformative.mixture"){
    w_inf_traj[k] <- sum(post_avg["w", 1:n_info_comps])
    }
}


 
    # extract decisions
    d_succ <- decisions_list[[k]]$success
    d_fut  <- decisions_list[[k]]$futility
 
    # Evaluate Efficacy Boundary (if defined for this stage)
    if (!is.null(d_succ)) {
      bnd_eff_fun  <- suppressMessages(RBesT::decision2S_boundary(prior_1, prior_2, n1_seq[k], n2_seq[k], d_succ))

      crit_val_eff <- bnd_eff_fun(y_2_curr[active_idx])
      is_lower_eff <- attr(d_succ, "lower.tail")
      is_succ <- if (is_lower_eff) y_1_curr[active_idx] <= crit_val_eff else y_1_curr[active_idx] > crit_val_eff
    } else {
      is_succ <- rep(FALSE, n_active)
    }
 
    # Evaluate Futility Boundary (if defined for this stage)
    if (!is.null(d_fut)) {
      bnd_fut_fun  <- suppressMessages(RBesT::decision2S_boundary(prior_1, prior_2, n1_seq[k], n2_seq[k], d_fut))

      crit_val_fut <- bnd_fut_fun(y_2_curr[active_idx])
      is_lower_fut <- attr(d_fut, "lower.tail")
      is_fut <- if (is_lower_fut) y_1_curr[active_idx] <= crit_val_fut else y_1_curr[active_idx] > crit_val_fut
    } else {
      is_fut <- rep(FALSE, n_active)
    }
 
    # Record stops, those are indices of active_idx
    just_stopped_eff <- is_succ
    just_stopped_fut <- is_fut & !is_succ  # shouldn't happen that both futility and efficacy are crossed, but anyway
    just_stopped     <- just_stopped_eff | just_stopped_fut
 
    # in the total vector of length n_sim, the ones that stopped now are full.vec[active_idx[stop.now]]
    stop_eff[active_idx[just_stopped_eff]] <- TRUE
    stop_fut[active_idx[just_stopped_fut]] <- TRUE
    stage_stopped[active_idx[just_stopped]] <- k  # at which stage does this simulated trial stop?
 
    # Record prob to stop (marginal over all n_sim trials, matching the analytical function)
    p_succ_vec[k] <- sum(just_stopped_eff) / n_sim
    p_fut_vec[k]  <- sum(just_stopped_fut) / n_sim
 
    # Update active trials
    active_idx <- active_idx[!just_stopped]
    n_2_prev <- n2_seq[k]
    n_1_prev <- n1_seq[k]
  }
 
  # find Exp value for rejection or for not reject
  idx_succ <- which(stop_eff)
  idx_fail <- which(!stop_eff)
 
  EN_1_succ <- if (length(idx_succ) > 0) mean(n1_seq[stage_stopped[idx_succ]]) else NA
  EN_2_succ <- if (length(idx_succ) > 0) mean(n2_seq[stage_stopped[idx_succ]]) else NA
 
  EN_1_fail <- if (length(idx_fail) > 0) mean(n1_seq[stage_stopped[idx_fail]]) else NA
  EN_2_fail <- if (length(idx_fail) > 0) mean(n2_seq[stage_stopped[idx_fail]]) else NA
 
  per_stage_df <- data.frame(
    Stage      = 1:K,
    N_Trt      = n1_seq,
    N_Ctrl     = n2_seq,
    P_Succ     = p_succ_vec,
    P_Fut      = p_fut_vec,
    Cum_P_Succ = cumsum(p_succ_vec),
    Cum_P_Fut  = cumsum(p_fut_vec),
    Exp_info_weight = w_inf_traj
  )
 



  # --------------
  # calculate final metrics 
  # --------------

  power_est <- sum(stop_eff) / n_sim
  fut_est <- sum(stop_fut) / n_sim
  
  # MC errors 
  mce_power <- sqrt((power_est * (1 - power_est)) / n_sim)
  mce_fut   <- sqrt((fut_est * (1 - fut_est)) / n_sim)
  
  mce_en_t  <- sd(n1_seq[stage_stopped]) / sqrt(n_sim)
  mce_en_c  <- sd(n2_seq[stage_stopped]) / sqrt(n_sim)

  # Return OC with MCE included
  return(list(
    Overall = c(
      Power        = power_est,
      MCE_Power    = mce_power,          
      Prob_Fut_seq = fut_est,
      MCE_Fut      = mce_fut,            
      EN_t         = mean(n1_seq[stage_stopped]),
      MCE_EN_t     = mce_en_t,           
      EN_c         = mean(n2_seq[stage_stopped]),
      MCE_EN_c     = mce_en_c,           
      EN_t_Succ    = EN_1_succ,
      EN_c_Succ    = EN_2_succ,
      EN_t_Fail    = EN_1_fail,
      EN_c_Fail    = EN_2_fail
    ),
    Per_Stage = per_stage_df
  ))

}
 
 
avgoc2_seq_mc.normMix <- function(prior_1, prior_2,
                                  n1_seq, n2_seq, decisions_list,   # recall that n1_seq is cumulative
                                  delta = 0, design_prior_c,
                                  sigma_1 = 88, sigma_2 = 88,       
                                  n_sim = 100000, seed = 123) {
 
  set.seed(seed)  
  theta_c_draws <- RBesT::rmix(design_prior_c, n_sim)
  theta_t_draws <- theta_c_draws + delta
 

  return(oc2_seq_mc.normMix(
    theta_1       = theta_t_draws,
    theta_2       = theta_c_draws,
    prior_1       = prior_1,
    prior_2       = prior_2,
    n1_seq        = n1_seq,
    n2_seq        = n2_seq,
    decisions_list = decisions_list,
    sigma_1       = sigma_1,
    sigma_2       = sigma_2,
    n_sim         = n_sim,
    seed          = NULL
  ))
}
 



























 # -----------------
 # FORMAT results and EUII



format.results <- function(res) {
 
  per_stage <- bind_rows(
    lapply(names(res), function(nm) {
      res[[nm]]$Per_Stage |>
        mutate(delta = str_remove(nm, "delta\\.") |> as.numeric())
    })
  )
 
  # ---- Overall results ----
  overall <- bind_rows(
    lapply(names(res), function(nm) {
      as.data.frame(t(res[[nm]]$Overall)) |>
        mutate(delta = str_remove(nm, "delta\\.") |> as.numeric())
    })
  )
 
  K <- max(per_stage$Stage)
 
  # ---- Per-stage table (one row per delta x stage) ----
  
  per_stage_tab <- per_stage |>
    mutate(
      Stage_Label = ifelse(Stage == K, "Final", paste0("Interim_", Stage)),
      P_Ind       = 1 - P_Succ - P_Fut,
      Cum_P_Ind   = 1 - Cum_P_Succ - Cum_P_Fut,
      across(c(P_Succ, P_Fut, P_Ind, Cum_P_Succ, Cum_P_Fut, Cum_P_Ind),
             ~ round(100 * .x, 1))
    ) |>
    select(delta, Stage_Label, N_Trt, N_Ctrl,
           P_Succ, P_Fut, P_Ind,
           Cum_P_Succ, Cum_P_Fut, Cum_P_Ind)  |> 
    pivot_wider(names_from = Stage_Label, 
    values_from = c(
      N_Trt, N_Ctrl,
      P_Succ, P_Fut, P_Ind,
      Cum_P_Succ, Cum_P_Fut, Cum_P_Ind
    ),
    names_sep = "."
    )
 
  # ---- Overall summary (one row per delta) ----
  # Overall power and expected sample sizes.
  overall_tab <- overall |>
    mutate(
      Expected_N = EN_t + EN_c,
      across(where(is.numeric), ~ round(.x, 1))
    ) |>
    select(delta, Power, EN_t, EN_c, Expected_N,
           EN_t_Succ, EN_c_Succ, EN_t_Fail, EN_c_Fail)
 
  return(list(
    per_stage = per_stage_tab,
    overall   = overall_tab
  ))
}
 
 
compute_euii <- function(res) {
 
  # --- build stage-level table ----
  per_stage <- bind_rows(
    lapply(names(res), function(nm) {
      res[[nm]]$Per_Stage |>
        mutate(delta = sub("delta\\.", "", nm) |> as.numeric())
    })
  )
 
  # --- overall table ----
  overall <- bind_rows(
    lapply(names(res), function(nm) {
      as.data.frame(t(res[[nm]]$Overall)) |>
        mutate(delta = str_remove(nm, "delta\\.") |> as.numeric())
    })
  )
 
  # Type I error = power under delta = 0
  t1e <- overall |>
    filter(delta == 0) |>
    pull(Power)
 
  if (length(t1e) == 0) {
    stop("delta = 0 must be included in results to compute EUII (needed for Type I error).")
  }
 
  K <- max(per_stage$Stage)
 
  # ---- compute ----
  # nonsig_prob per stage:
  #   interim stages (k < K): non-significant stops are futility only (P_Fut)
  #   final stage   (k = K): non-significant = residual probability of reaching
  #                           the final stage and not succeeding = 1 - Cum_P_Succ - Cum_P_Fut
  #

  out <- per_stage |>
    mutate(
      N_total     = N_Trt + N_Ctrl,
      nonsig_prob = ifelse(
        Stage == K,
        1 - Cum_P_Succ - Cum_P_Fut,   # final stage: residual
        P_Fut                           # interim stages: futility stops only
      )
    ) |>
    group_by(delta) |>
    summarise(
 
      Power = Cum_P_Succ[Stage == K],
 
      # Likelihood ratios
      LR_pos = Power / t1e,
      LR_neg = (1 - Power) / (1 - t1e),
 
      # E(1/N | significant):
      # sum over all stages of P(succeed at stage k) / Power * 1/N_k
      E_invN_sig =
        sum(P_Succ / Power * (1 / N_total)),
 
      # E(1/N | non-significant):
      # sum over all stages of P(non-sig stop at stage k) / (1-Power) * 1/N_k
      E_invN_nonsig =
        sum(nonsig_prob / (1 - Power) * (1 / N_total)),
 
      .groups = "drop"
 
    ) |>
    mutate(
      EUII =
        LR_pos^E_invN_sig /
        LR_neg^E_invN_nonsig
    )
 
  out
}