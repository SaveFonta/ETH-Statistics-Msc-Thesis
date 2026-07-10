library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(RBesT)


# -----------------------------------------------------
# FUNCTIONS for Monte Carlo Sequential evaluation
# -----------------------------------------------------


# The decision boundary is a function of the ANALYSIS: it is
# built from the priors, the sample sizes and the decision rule, and takes
# no theta and no delta.
build_boundaries <- function(prior_1, prior_2, n1_seq, n2_seq, decisions_list) {
  lapply(seq_along(n2_seq), function(k) {
    d_succ <- decisions_list[[k]]$success
    d_fut  <- decisions_list[[k]]$futility
    list(
      succ = if (!is.null(d_succ))
        suppressMessages(RBesT::decision2S_boundary(prior_1, prior_2, n1_seq[k], n2_seq[k], d_succ)),
      fut  = if (!is.null(d_fut))
        suppressMessages(RBesT::decision2S_boundary(prior_1, prior_2, n1_seq[k], n2_seq[k], d_fut))
    )
  })
}


oc2_seq_mc.normMix <- function(theta_1, theta_2, prior_1, prior_2,
                                n1_seq, n2_seq, decisions_list,   # recall that n1_seq is cumulative
                                sigma_1 = 88, sigma_2 = 88,
                                n_sim = 100000, seed = NULL,
                                boundaries = NULL,   # from build_boundaries()
                                weight.track = FALSE,
                                type.weight.track = "uninformative.mixture", n_info_comps = 3)  {     # added these 2 params: so we can track the informative weight both if only one uninformative component (given by robustify)
                                                                                     # or even if we use a mixture as non-informative (e.g. 100 mmixture to approx a t-distri)
 
  if (!is.null(seed)) set.seed(seed)
 
  # If we are computing conditional, there is only one theta.
  # Instead for the average metric we will provide a vector of theta sampled from pi_D
  if (length(theta_1) == 1) theta_1 <- rep(theta_1, n_sim)
  if (length(theta_2) == 1) theta_2 <- rep(theta_2, n_sim)
 
  # check length 
  if (length(theta_1) != n_sim) stop("Number of theta_1 provided should be either one for classic OC or equal to 
                                      n_sim for avg OC")
  if (length(theta_2) != n_sim) stop("Number of theta_2 provided should be either one for classic OC or equal to 
                                      n_sim for avg OC")
 
  K <- length(n2_seq) # number of stages
 
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
    # d_succ is a single decision2S object. For a dual criterion it carries
    # several criteria via vector pc/qc, e.g.
    #   decision2S(pc = c(0.95, 0.50), qc = c(0, -50), lower.tail = TRUE)
    d_succ <- decisions_list[[k]]$success
    d_fut  <- decisions_list[[k]]$futility

    # Evaluate Efficacy Boundary (if defined for this stage)
    # Reuse the prebuilt boundary when one was supplied, otherwise build it here.
    if (!is.null(d_succ)) {
      bnd_eff_fun  <- if (!is.null(boundaries)) boundaries[[k]]$succ else
        suppressMessages(RBesT::decision2S_boundary(prior_1, prior_2, n1_seq[k], n2_seq[k], d_succ))

      crit_val_eff <- bnd_eff_fun(y_2_curr[active_idx]) # find the critical value for each trial using the closure decision2S_boundary
      is_lower_eff <- attr(d_succ, "lower.tail") # lower is better or larger is better?
      is_succ <- if (is_lower_eff) y_1_curr[active_idx] <= crit_val_eff else y_1_curr[active_idx] > crit_val_eff
    } else {
      is_succ <- rep(FALSE, n_active)
    }
 
    # Evaluate Futility Boundary (if defined for this stage)
    if (!is.null(d_fut)) {
      bnd_fut_fun  <- if (!is.null(boundaries)) boundaries[[k]]$fut else
        suppressMessages(RBesT::decision2S_boundary(prior_1, prior_2, n1_seq[k], n2_seq[k], d_fut))

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
 
 ## Expected samp size given trial succeded
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
  
  ## MC errors 
  # mce_power <- sqrt((power_est * (1 - power_est)) / n_sim)
  # mce_fut   <- sqrt((fut_est * (1 - fut_est)) / n_sim)
  
  # mce_en_t  <- sd(n1_seq[stage_stopped]) / sqrt(n_sim)
  #mce_en_c  <- sd(n2_seq[stage_stopped]) / sqrt(n_sim)

  # Return OC (MC errors are available above, currently switched off)
  return(list(
    Overall = c(
      Power        = power_est,
      # MCE_Power    = mce_power,          
      Prob_Fut_seq = fut_est,
      # MCE_Fut      = mce_fut,            
      EN_t         = mean(n1_seq[stage_stopped]),
      # MCE_EN_t     = mce_en_t,           
      EN_c         = mean(n2_seq[stage_stopped]),
      # MCE_EN_c     = mce_en_c,           
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

  # Determine direction from the stage 1 success decision. It is a single
  # decision2S object. 

  ## NOTE: cannot have a design with no efficacy at the interim (decisions_list[[1]]$success == NULL)
  # Also, we are assuming that the direction for success is identical for every stage

  is_lower_tail <- attr(decisions_list[[1]][["success"]], "lower.tail")

  # Built once and shared by every delta: the boundary does not depend on theta or delta. 
  boundaries <- build_boundaries(prior_1, prior_2, n1_seq, n2_seq, decisions_list)

  theta_c_draws <- RBesT::rmix(design_prior_c, n_sim)

#function to process one specific delta value
# if lower is better, then delta = theta_C - theta_T (like the thesis)
# if lower is worse, then delta = theta_T - theta_C
  process <- function(d) {
   theta_t_draws <- if (is_lower_tail) theta_c_draws - d  else theta_c_draws + d
  
  oc2_seq_mc.normMix(
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
    seed          = NULL,
    boundaries    = boundaries
  )
  }


  if (length(delta) == 1) { #only one delta provided
    return(process(delta))
  } else {
    res_list <- lapply(delta, process)
    names(res_list) <- paste0("delta.",delta)
    return(res_list)
  }
}
# 

# -----------------------------------------------------
# EUII
# -----------------------------------------------------

compute_euii <- function(res, prior_H1 = c(0.01, 0.1, 0.5), eps = 1e-8) {
 
  # --- build stage-level table ----
  per_stage <- bind_rows(
    lapply(names(res), function(nm) { #names(res) is all the different true deltas 
      res[[nm]]$Per_Stage |>
        mutate(delta = sub("delta\\.", "", nm) |> as.numeric())
    })
  )
 

   K <- max(per_stage$Stage)

  per_stage_proc <- per_stage |>
    mutate(
      N_total     = N_Trt + N_Ctrl,
      #if we are at last stage, all that is not a futility stop or a success, is Undefeined outcome, but I decided to say it is a fail of the test (non reject)
      nonsig_prob = ifelse(Stage == K, 1 - Cum_P_Succ - Cum_P_Fut, P_Fut) 
    )


  pointwise_metrics <- per_stage_proc |>
    group_by(delta) |>
    summarise(
      Power = Cum_P_Succ[Stage == K],
      #  E(1/N | sign) and E(1/N | nonsign)
      E_invN_sig_cond = sum(P_Succ / Power * (1 / N_total)),
      E_invN_nonsig_cond = sum(nonsig_prob / (1 - Power) * (1 / N_total)),
      .groups = "drop"
    )

  # Power exactly 0 or 1 makes E(1/N | sign) or E(1/N | nonsign) undefined (0/0),
  bad <- pointwise_metrics$Power %in% c(0, 1)
  if (any(bad)) {
    cat("WARNING in compute_euii: Power is exactly 0 or 1 at delta =",
        paste(pointwise_metrics$delta[bad], collapse = ", "),
        "  the conditional E(1/N) is undefined there (increase n_sim).\n")
  }

 null <- pointwise_metrics |> filter(delta == 0)
  if (nrow(null) == 0) {
    stop("delta = 0 must be included to compute the EUII.")
  }
 
  T1E <- null$Power
  E_invN_sig_null <- null$E_invN_sig_cond
  E_invN_nonsig_null <- null$E_invN_nonsig_cond
 

 out <- lapply(prior_H1, function(x) {
    pointwise_metrics |>
    mutate(
      # Clamp Power and T1E away from 0 and 1, so every log below is finite.
      # The warning above fires exactly when this clamp happens.
      Power_c = pmin(pmax(Power, eps), 1 - eps),
      T1E_c   = pmin(pmax(T1E,   eps), 1 - eps),

      # Likelihood ratios on the log scale. 
      log_LR_pos = log(Power_c)    - log(T1E_c),           # log LR+ = log(Power / T1E)
      log_LR_neg = log1p(-Power_c) - log1p(-T1E_c),        # log LR- = log((1-Power) / (1-T1E))
      LR_pos     = exp(log_LR_pos),                        # kept for reporting only
      LR_neg     = exp(log_LR_neg),

      # Posterior odds of H_1 on the log scale: log O(H_1 | outcome) = log LR + logit(prior_H1).
      # since Posterior Odds = LR * Prior Odds
      log_post_odds_sig = log_LR_pos + qlogis(x),   # log O(H_1 | sign)
      Pr_H1_sig         = plogis( log_post_odds_sig),  # P(H_1 | sign)
      Pr_H0_sig         = plogis(-log_post_odds_sig),  # P(H_0 | sign)

      log_post_odds_nonsig = log_LR_neg + qlogis(x),   # log O(H_1 | non sign)
      Pr_H1_nonsig         = plogis( log_post_odds_nonsig), # P(H_1 | nonsign)
      Pr_H0_nonsig         = plogis(-log_post_odds_nonsig),


      E_invN_sig = (E_invN_sig_null * Pr_H0_sig) + (E_invN_sig_cond * Pr_H1_sig),
      E_invN_nonsig = (E_invN_nonsig_null * Pr_H0_nonsig) + (E_invN_nonsig_cond * Pr_H1_nonsig)
    ) |>
    mutate(
        log_EUII = (E_invN_sig * log_LR_pos) - (E_invN_nonsig * log_LR_neg), # log scale
        EUII     = exp(log_EUII)
            )  |>
    select ("Delta" = delta,  LR_pos, LR_neg, E_invN_sig, E_invN_nonsig, log_EUII, EUII)
 }
 )

 names(out) <- as.character(prior_H1)
  return(out)
}



# -----------------------------------------------------
# QUICK EXAMPLE
# -----------------------------------------------------


if (FALSE) {

  sigma   <- 88
  p_MAP   <- mixnorm(
    c(0.4848, -52.457, 21.154), c(0.4598, -47.465, 7.843), c(0.0554, -50.355, 48.164),
    sigma = sigma, param = "ms"
  )
  p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")

  # Efficacy: dual criterion carried by a single decision2S (vector pc/qc).
  sign.crit <- decision2S(pc = c(0.95, 0.50), qc = c(0, -50), lower.tail = TRUE)
  # Futility: Pr(delta < 40 | data) >= 0.90, i.e. Pr(theta_T - theta_C > -40) >= 0.90.
  fut.crit  <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)

  decisions_fut <- list(
    list(success = sign.crit, futility = fut.crit),   # stage 1: efficacy + futility
    list(success = sign.crit, futility = NULL)        # stage 2: final, efficacy only
  )

  n1_seq <- c(30, 60)   # treatment arm, cumulative
  n2_seq <- c(20, 40)   # control arm, cumulative

  # Conditional OC: a scalar theta fixes theta_C.
  res_cond <- oc2_seq_mc.normMix(
    theta_1 = -50, theta_2 = -50,
    prior_1 = p_vague, prior_2 = p_MAP,
    n1_seq  = n1_seq,  n2_seq  = n2_seq,
    decisions_list = decisions_fut,
    n_sim = 1000, seed = 1
  )
  res_cond$Overall
  res_cond$Per_Stage

  # Average OC: theta_C is drawn from the design prior, one draw per simulated
  # trial. delta = 0 must be in the grid, since compute_euii reads the average
  # T1E from it.
  res_avg <- avgoc2_seq_mc.normMix(
    prior_1 = p_vague, prior_2 = p_MAP,
    n1_seq  = n1_seq,  n2_seq  = n2_seq,
    decisions_list = decisions_fut,
    delta          = c(0, 30, 60),
    design_prior_c = p_MAP,
    n_sim = 1000, seed = 123
  )

  ## --> avgoc2_seq_mc.normMix spirs out a list for every delta. 
  ## Inside that list there are two df: $Overall and $Per_Stage
res_avg$delta.0


  # Sequential EUII. prior_H1 enters only here
  res_euii <- compute_euii(res_avg, prior_H1 = c(0.01, 0.1, 0.5))

  res_euii[["0.01"]]
}





