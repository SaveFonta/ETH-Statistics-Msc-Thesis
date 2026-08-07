library(dplyr)
library(RBesT)



# -----------------------------------------------------
# EXACT METHODS (fixed design, closed-form / numerical integration)
# -----------------------------------------------------

avgoc2S.normMix <- function(
  prior1, prior2, n1, n2, decision, delta, design_prior2,   eps = 1e-6, Ngrid = 10, ...
) {

  # Creates OC
  oc_fun <- RBesT::oc2S(prior1, prior2, n1, n2, decision, eps = eps, ...)

    # warm up cache
  lim2_range <- RBesT::qmix(design_prior2, c(eps / 2, 1 - eps / 2))
  lim1_range <- lim2_range + delta

  # we know y_1 ~ N(theta_1, SEM_1), so in theory we should warm up cache for y_1 not for theta1
  # this is not necessary cause inside oc2s this is alraedy handled

  invisible(oc_fun(lim1_range, lim2_range))


  # output function
  design_fun <- function(delta_new = delta, design_prior2_new = design_prior2) {

    res <- RBesT:::integrate_density_log(
      log_integrand = function(x) {
        log(oc_fun(x + delta_new, x))
      },
      mix = design_prior2_new,
      Lplower = RBesT::logit(eps / 2),
      Lpupper = RBesT::logit(1 - eps / 2)
    )

    return(res)
  }

  return(design_fun)
}


# -----------------------------------------------------
# FUNCTIONS for Monte Carlo Sequential evaluation
# -----------------------------------------------------



# --- build_boundaries ----
# The decision boundary is a function of the ANALYSIS: it is
# built from the priors, the sample sizes and the decision rule, and takes
# no theta and no delta.

# Useful such that in each list of boundaries, instead of the decision2S object, there is the 
# function that takes as input a value of y2 and gives the boundary!! 
# Since it is a closure, cache a lot and makes things faster 

build_boundaries <- function(prior_1, prior_2, n1_seq, n2_seq, decisions_list) { 
  
  lapply(seq_along(n1_seq), function(stage) {
    
    individual_suc <- decisions_list[[stage]]$success 
    individual_fut <- decisions_list[[stage]]$futility

    n1 <- n1_seq[stage]
    n2 <- n2_seq[stage]
    
    list(
        succ = if (!is.null(individual_suc)) {
          suppressMessages(RBesT::decision2S_boundary(prior1 = prior_1, prior2 = prior_2, n1 = n1, n2 = n2, decision = individual_suc)) # <-- Fixed prior_1 / prior_2
        } else {NULL},
        
        fut = if (!is.null(individual_fut)) {
          suppressMessages(RBesT::decision2S_boundary(prior1 = prior_1, prior2 = prior_2, n1 = n1, n2 = n2, decision = individual_fut)) # <-- Fixed prior_1 / prior_2
        } else {NULL}
    )
  })
}






oc2_seq_mc.normMix <- function(theta_1, theta_2, prior_1, prior_2,
                                n1_seq, n2_seq, decisions_list,   # recall that n1_seq is cumulative
                                sigma_1 = 88, sigma_2 = 88,
                                n_sim = 1e5, seed = NULL,
                                boundaries = NULL,   # from build_boundaries()
                                weight.track = FALSE,
                                type.weight.track = "uninformative.mixture", n_info_comps = 3)  {     # added these 2 params: so we can track the informative weight both if only one uninformative component (given by robustify)
                                                                                     # or even if we use a mixture as non-informative (e.g. 100 mmixture to approx a t-distri)
 
  if (!is.null(seed)) set.seed(seed)
 
  # If we are computing conditional, there is only one theta for each of the n_sim simulated trials
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
  


  # Return OC 
  return(list(
    Overall = c(
      Power        = power_est,
      Prob_Fut_seq = fut_est,
      EN_t         = mean(n1_seq[stage_stopped]),
      EN_c         = mean(n2_seq[stage_stopped]),
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
                                  n_sim = 1e5, seed = 123) {
 
  set.seed(seed)  

  # Determine direction from the stage 1 success decision. It is a single
  # decision2S object. 

  ## NOTE: cannot have a design with no efficacy at the first interim (decisions_list[[1]]$success == NULL)
  # Also, we are assuming that the direction for success is identical for every stage (which makes sense....)

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
  per_stage <- lapply(names(res), function(nm) {  # names(res) is all the different true deltas 
      res[[nm]]$Per_Stage |>
        mutate(delta = sub("delta\\.", "", nm) |> as.numeric())
    })  |> bind_rows()
 

   K <- max(per_stage$Stage)

  per_stage_proc <- per_stage |>
    mutate(
      N_total     = N_Trt + N_Ctrl,
      # if we are at last stage, all that is not a futility stop or a success, is Undefined outcome, 
      # anyway, according to the thesis, we follow the convention that it is a fail of the test (non reject)
      nonsig_prob = ifelse(Stage == K, 1 - Cum_P_Succ - Cum_P_Fut, P_Fut) 
    )


  pointwise_metrics <- per_stage_proc |>
    group_by(delta) |>
    summarise(
      Power = Cum_P_Succ[Stage == K],
      #  E(1/N | +) and E(1/N | -)
      E_invN_sig_cond = sum(P_Succ / Power * (1 / N_total)),
      E_invN_nonsig_cond = sum(nonsig_prob / (1 - Power) * (1 / N_total)),
      .groups = "drop"
    )

  # Power exactly 0 or 1 makes E(1/N | sign) or E(1/N | nonsign) undefined (0/0),
  bad <- pointwise_metrics$Power %in% c(0, 1)
  if (any(bad)) {
    cat("WARNING in compute_euii: Power is exactly 0 or 1 at delta =",
        paste(pointwise_metrics$delta[bad], collapse = ", "),
        "  the conditional E(1/N) is undefined there (try to increase n_sim).\n")
  }

 null <- pointwise_metrics |> filter(delta == 0)
  if (nrow(null) == 0) {
    stop("delta = 0 must be included to compute the EUII.")
  }
 
  T1E <- null$Power
  E_invN_sig_null <- null$E_invN_sig_cond
  E_invN_nonsig_null <- null$E_invN_nonsig_cond

  # T1E as a proper column (recycled scalar), so the select below can keep it.
  # Power is already a column from the summarise above.
  pointwise_metrics <- pointwise_metrics |> mutate(T1E = null$Power)


 out <- lapply(prior_H1, function(x) {
    pointwise_metrics |>
    mutate(
      # Clamp Power and T1E away from 0 and 1, so every log below is finite.
      # The warning above fires exactly when this clamp would happen, so I don't think this will
      # ever fire
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


    # Law of total exp: E[1/N | +] = E[1/N | +, H0] P[H_0 | +] + E[1/N | +, H1] P[H_1 | +] 
      E_invN_sig = (E_invN_sig_null * Pr_H0_sig) + (E_invN_sig_cond * Pr_H1_sig), 
      E_invN_nonsig = (E_invN_nonsig_null * Pr_H0_nonsig) + (E_invN_nonsig_cond * Pr_H1_nonsig)
    ) |>
    mutate(
        log_EUII = (E_invN_sig * log_LR_pos) - (E_invN_nonsig * log_LR_neg), # log scale
        EUII     = exp(log_EUII)
            )  |>
    select ("Delta" = delta,  Power, T1E, LR_pos, LR_neg, E_invN_sig, E_invN_nonsig, log_EUII, EUII)
 }
 )

 names(out) <- as.character(prior_H1)
  return(out)
}



# -----
# One MC sweep over every (analysis prior x design prior) combination in `combos`.

# Outputs list of avgoc2_seq_mc.normMix() results, one per row of `combos` (each
#         itself is a list varying deltas), named "<Analysis_Prior> | <Design_Prior>"
#         so a specific combo can be indexed directly, e.g. sweeps$SC[["MAP | Dirac (-50)"]].
run_sequential_sweep <- function(decisions_list, combos, analysis_priors, design_priors,
                                 prior_1, n1_seq, n2_seq, delta_values,
                                 n_sim_avg = 1e6, mc_seed = 123, cores = 6) {

  if (.Platform$OS.type != "unix") cores <- 1L   # mclapply cannot fork on Windows

  cat("\n", nrow(combos), " prior combinations, n_sim = ", n_sim_avg,
      ", cores = ", cores, "\n", sep = "")

  res <- parallel::mclapply(seq_len(nrow(combos)), function(i) {
    avgoc2_seq_mc.normMix(
      prior_1        = prior_1,
      prior_2        = analysis_priors[[combos$Analysis_Prior[i]]],
      n1_seq         = n1_seq,
      n2_seq         = n2_seq,
      decisions_list = decisions_list,
      delta          = delta_values,
      design_prior_c = design_priors[[combos$Design_Prior[i]]],
      n_sim          = n_sim_avg,
      seed           = mc_seed
    )
  }, mc.cores = cores)

  # Name each element by its combo, e.g. "MAP | Dirac (-50)", so can extract like (sweeps$SC[["MAP | Dirac (-50)"]]) 
  names(res) <- paste(combos$Analysis_Prior, combos$Design_Prior, sep = " | ")

  # mclapply reports failures as try-error elements
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) {
    cat("MC failed for the following prior combinations:\n")
    cat(paste0("   ", combos$Analysis_Prior[bad], " | ", combos$Design_Prior[bad], "\n"), sep = "")
    stop("MC failed for ", sum(bad), " of ", nrow(combos), " prior combinations.")
  }
  res
}


# extract every delta's result within every combo, then stacks the rows.

# `extract` returns one data.frame per delta, so summarise outputs
# single data.frame stacking extract()'s rows across every delta and every
# combo, with Analysis_Prior/Design_Prior columns attached.
summarise_sweep <- function(sweeped, combos, extract) {

  list_res <- lapply(seq_along(sweeped), function(i) {
    res <- sweeped[[i]]
    lapply(names(res), function(d) extract(res[[d]], d)) |> bind_rows()
  })

  lapply(seq_len(nrow(combos)), function(combo) {
    list_res[[combo]] |>
      mutate(
        Analysis_Prior = combos$Analysis_Prior[[combo]],
        Design_Prior   = combos$Design_Prior[[combo]]
      )
  }) |> bind_rows()
}


# Extractors to use in summarise_sweep(): pull one delta's result down to the row(s)
# needed for the total OC table (df_oc) or the per-stage table (df_stage).
# Output: -> one-row data.frame: delta, Power, EN.
extract_sweep_total <- function(res_d, d) {
  ov <- res_d$Overall
  data.frame(delta = sub("delta\\.", "", d) |> as.numeric(),
             Power = ov[["Power"]],
             EN    = ov[["EN_t"]] + ov[["EN_c"]])   # expected total sample size
}
# Output: -> one row per stage: Stage, P_Succ, P_Fut, Cum_P_Succ, Cum_P_Fut, delta.
extract_sweep_per_stage <- function(res_d, d) {
  res_d$Per_Stage |>
    select(Stage, P_Succ, P_Fut, Cum_P_Succ, Cum_P_Fut) |>
    mutate(delta = sub("delta\\.", "", d))
}


# compute EUII (every prior_H1) across all combos in one sweep.
# Output a data.frame with EUII/LR/E[1/N] columns for every combo x prior_H1
# combination (see compute_euii() for the column list).
sweep_euii_table <- function(sweep, combos, prior_H1, omega_map, prior_names, dprior_names) {
  bind_rows(lapply(seq_along(sweep), function(i) {
    e <- compute_euii(sweep[[i]], prior_H1 = prior_H1)
    bind_rows(lapply(names(e), function(g) {
      e[[g]] |>
        mutate(prior_H1       = as.numeric(g),
               Analysis_Prior = combos$Analysis_Prior[i],
               Design_Prior   = combos$Design_Prior[i],
               omega          = omega_map[[combos$Analysis_Prior[i]]])
    }))
  })) |>
    mutate(Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
           Design_Prior   = factor(Design_Prior,   levels = dprior_names))
}



# Read the two fixed EUII baselines (SC and DC) at one analysis/design
# prior combination.
# Output: data.frame with Delta, EUII_fixed, Base ("SC"/"DC"), one row per delta per design.
read_fixed_baselines <- function(sc_path = "Output/02_fixed_design_SC/euii.RDS",
                                 dc_path = "Output/03_fixed_design_DC/euii.RDS",
                                 analysis_prior = "MAP", design_prior = "MAP") {
  read_one <- function(path, base_label) {
    readRDS(path)$df_euii |>
      filter(Analysis_Prior == analysis_prior, Design_Prior == design_prior) |>
      mutate(Delta = delta, EUII_fixed = EUII, Base = base_label) |>
      select(Delta, EUII_fixed, Base)
  }
  bind_rows(read_one(sc_path, "SC"), read_one(dc_path, "DC"))
}


# One panel per Criterion, coloured by Analysis_Prior, faceted layout shared by
# the avgT1E and avgPower plots.
oc_facet_layout <- function(p) {
  p +
    geom_point(size = 2.2) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~Criterion, nrow = 1) +
    scale_linetype_vague() +
    scale_color_prior() +
    my_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}


# Summarise a value column into a reference line (prior_H1 = ph1_ref) plus a
# min/max ribbon over the rest of the prior_H1 grid, grouped by `group_cols`.
value_band <- function(df, value_col, group_cols, ph1_ref = 0.5) {
  df |>
    group_by(across(all_of(group_cols))) |>
    summarise(lo  = min(.data[[value_col]]),
              hi  = max(.data[[value_col]]),
              mid = .data[[value_col]][prior_H1 == ph1_ref],
              .groups = "drop")
}




# -----------------------------------------------------
# Threshold calibration for Bayesian group sequential designs
#
# Find the common posterior probability threshold p, 
# such that the avgT1E of the sequential design equals a target. 
# decision_builder(p) handles the parameter p. By default p is like Pocock constant way
#
# To use uniroot, every evaluation reuses the same seed, so 
# T1E(p) is a monotone function of p (decreasing with p)
#
# criterion = "SC": success is Pr(delta > q | data) >= p.
# criterion = "DC": success is Pr(delta > q | data) >= p AND
#                   Pr(delta > DV | data) >= p_DV, with the DV part held fixed.
#   For the DC the target may be unattainable, since the DV caps the T1E from above
# 
# lower_is_better = TRUE (default): lower outcome values are better, so
#   delta = theta_C - theta_T as in the thesis. FALSE: delta = theta_T - theta_C.
# sigma_1, sigma_2: known sampling standard deviations (treatment, control).
# delta_null: true effect at which the T1E is computed (usually 0).

# decisions_builder:
#   A function(p) returning a decisions_list of length K, exactly like the
#   one passed to avgoc2_seq_mc.normMix. When given, it overrides criterion,
#   q, DV, p_DV and lower_is_better: only p is calibrated, everything else is
#   whatever the builder does with it, useful for example with unequal tresholds, 
#  or futility rules or if you want stage specific criteria --> SEE BOTTOM OF THIS FUNCTION for some examples


# -----------------------------------------------------

calibrate_threshold <- function(target_t1e,
                                n1_seq, n2_seq,
                                prior_cntrl,            # control analysis prior
                                prior_treat,            # set to vague usually
                                design_prior = prior_cntrl,  # aligned by default, would be weird to calibrate our design assuming design prior different from an prior
                                criterion    = c("SC", "DC"),
                                q  = 0,   # NV (null value)
                                DV = 50,
                                p_DV = 0.5, # usually the second criterion is a DC, so the p is 0.5
                                p_interval = c(0.6, 0.999),
                                sigma_1 = 88,
                                sigma_2 = 88,
                                delta_null = 0,
                                n_sim = 1e5,
                                seed  = 123,
                                tol   = 1e-4,
                                lower_is_better = TRUE,
                                decisions_builder = NULL) {  
  criterion <- match.arg(criterion)
  K <- length(n2_seq)

  # decisions list for a given threshold p, efficacy at every look, no futility, in case decision_builder not provided
  make_decisions <- if (!is.null(decisions_builder)) decisions_builder else function(p) {
    crit <- if (criterion == "SC") {
      if (lower_is_better) {
        RBesT::decision2S(pc = p, qc = -q, lower.tail = TRUE)
      } else {
        RBesT::decision2S(pc = p, qc = q, lower.tail = FALSE)
      }
    } else { # DC
      if (lower_is_better) {
        RBesT::decision2S(pc = c(p, p_DV), qc = c(-q, -DV), lower.tail = TRUE)
      } else {
        RBesT::decision2S(pc = c(p, p_DV), qc = c(q, DV), lower.tail = FALSE)
      }
    }
    lapply(seq_len(K), function(k) list(success = crit, futility = NULL)) # don't calibrate for futility and default is Pocock way
  }



  # avgT1E at threshold p
  t1e_at <- function(p) {
    res <- avgoc2_seq_mc.normMix(
      prior_1        = prior_treat,
      prior_2        = prior_cntrl,
      n1_seq         = n1_seq,
      n2_seq         = n2_seq,
      decisions_list = make_decisions(p),
      delta          = delta_null,
      design_prior_c = design_prior,
      sigma_1        = sigma_1,
      sigma_2        = sigma_2,
      n_sim          = n_sim,
      seed           = seed #setting the seed is crucial s.t. this is a function
    )
    res$Overall[["Power"]]
  }

  t1e_lo <- t1e_at(p_interval[1])   # most liberal threshold: largest T1E
  t1e_hi <- t1e_at(p_interval[2])   # most strict threshold: smallest T1E
  
  ## As we increase p, the t1e will always become smaller 

  # if the largest t1e we can find is still smaller than the target, there is no chance
  if (t1e_lo < target_t1e) {
    warning("target T1E ", target_t1e, " above what p = ", p_interval[1],
            " reaches (", round(t1e_lo, 5), "), widen p_interval")
    out <- list(p = p_interval[1], achieved_t1e = t1e_lo,
                target_t1e = target_t1e, attained = FALSE,
                criterion = criterion, decisions = make_decisions(p_interval[1]))
    return(out)
  }

  # if the smaller t1e we can find is still larger than the target, no chance.
  if (t1e_hi > target_t1e) {
    warning("target T1E ", target_t1e, " below what p = ", p_interval[2],
            " reaches (", round(t1e_hi, 5), "), widen p_interval")
    out <- list(p = p_interval[2], achieved_t1e = t1e_hi,
                target_t1e = target_t1e, attained = FALSE,
                criterion = criterion, decisions = make_decisions(p_interval[2]))
    return(out)
  }
  
  root <- uniroot(function(p) t1e_at(p) - target_t1e,
                  interval = p_interval, tol = tol)
  p_star   <- root$root
  achieved <- t1e_at(p_star)

  list(p = p_star, achieved_t1e = achieved, target_t1e = target_t1e,
       attained = TRUE, criterion = criterion, decisions = make_decisions(p_star))
}

# -----------------------------------------------------
# Joint calibration of threshold and sample size
#
# Find the smallest n_max (Note that we keep allocation ratio and informatio ratio constant) 
# and the common threshold p such that
#   avgT1E   = target_t1e     (via calibrate_threshold, at every candidate n)
#   avgPower >= target_power  (at delta = delta_power)

# Basically we try to find the n_max to reach a specific avgPower. Then with that n_max
# we adjust p to also achieve that specific avgT1E

# -----------------------------------------------------

calibrate_design <- function(target_t1e, target_power, delta_power,
                             n1_base, n2_base,          
                             prior_cntrl, prior_treat,
                             design_prior = prior_cntrl,
                             criterion = c("SC", "DC"),
                             q = 0, DV = 50, p_DV = 0.5,
                             p_interval = c(0.6, 0.999),
                             sigma_1 = 88, sigma_2 = 88,
                             delta_null = 0,
                             n_sim = 1e5, seed = 123,
                             tol = 1e-4,
                             n2_max_cap = 64 * max(n2_base),
                             lower_is_better = TRUE,
                             decisions_builder = NULL) {

  criterion <- match.arg(criterion)

  # Find initial parameters:
  allo_ratio  <- n1_base[length(n1_base)] / n2_base[length(n2_base)]
  info_ratio <- n2_base / n2_base[length(n2_base)] 

 ## IDEA: keep allocation ratio (e.g. 2:1) and information ratio (e.g. c(0.20,0.80,1)) constant, just vary n_max

  
  ## --- Define functions ---

  # given the n_2max (n_max of control), use allo_ratio and info_ratio to compute the new schedule of looks 
  # the result is a list with n_1 and n_2 schedules for each interim 
  schedule <- function(n_2max) {
    n2 <- pmax(1, round(info_ratio * n_2max)) 
    n2 <- cummax(n2) # so that if by rounding c(10,9,20) --> c(10,10,20)
    if (any(diff(n2) <= 0)) n2 <- n2[1] + cumsum(c(0, pmax(1, diff(n2)))) #if we have c(10,10,20) --> c(10,11,20), forces every interim to collect at least one patient
    n1 <- pmax(1, round(allo_ratio * n2))
    list(n1 = n1, n2 = n2)
  }

  # Given a (n_2max), first generate the schedule using the function
  # then calibrate p at this schedule to target avgT1E
  # then return achieved power at delta_power

  eval_at <- function(n_2max) {
    ns  <- schedule(n_2max)

    # First we calibrate this schedule for the desired t1e
    cal <- calibrate_threshold(target_t1e,
                               n1_seq = ns$n1, n2_seq = ns$n2,
                               prior_cntrl = prior_cntrl, prior_treat = prior_treat,
                               design_prior = design_prior,
                               criterion = criterion, q = q, DV = DV, p_DV = p_DV,
                               p_interval = p_interval,
                               sigma_1 = sigma_1, sigma_2 = sigma_2,
                               delta_null = delta_null,
                               n_sim = n_sim, seed = seed, tol = tol,
                               lower_is_better = lower_is_better,
                               decisions_builder = decisions_builder)
    if (!cal$attained) {
      stop("at n_2max = ", n_2max, " the T1E target ", target_t1e,
           " is unattainable (achievable ", round(cal$achieved_t1e, 5),
           ").")
    }

    # Now we find the value of avgPower at that specific Delta
    pw <- avgoc2_seq_mc.normMix(
      prior_1 = prior_treat, prior_2 = prior_cntrl,
      n1_seq = ns$n1, n2_seq = ns$n2,
      decisions_list = cal$decisions,
      delta = delta_power, design_prior_c = design_prior,
      sigma_1 = sigma_1, sigma_2 = sigma_2,
      n_sim = n_sim, seed = seed
    )$Overall[["Power"]]

    # out
    list(power = pw, cal = cal, ns = ns)
  }

  # --- Use of functions ----

  # we start from the provided n_max. (last element of n2_base)
  n_max_curr <- n2_base[length(n2_base)] # n_2max provided in the call. It is the value for which we start 

  # calibrate the t1e and see which power value we obtain
  e_curr <- eval_at(n_max_curr) 

  

  ## IDEA: we have eval_at(n_max)[["power"]] which outputs the avgPower
  #  The result is an increasing function in n_max. 
  # So we need to find the root of the function:  e_curr(n_max)[["power"]] - targetPower

  # To do that, we can find the interval where this root lies. So the n_max for which 
  # avgPower is larger than target and n_max for which avgPow is smaller than target.  
  

  ## If the current power value is larger, we can halve the n2max till we are below the target power 
  # at that point we have the value for n2max that can be a lower bound for the target one  

  if (e_curr$power >= target_power) {          # search downward bound
    # Upper bound:
    hi <- n_max_curr; e_hi <- e_curr

    ## Halve the n_max (while always staying above a sample size of 2)
    lo <- max(2, floor(n_max_curr / 2))
    e_lo <- eval_at(lo)
    
    ## keep going to halve if we still are above
    while (e_lo$power >= target_power && lo > 2) {

      # old lower bound becomes upper, andd we will repeat the halving
      hi <- lo; e_hi <- e_lo
      lo <- max(2, floor(lo / 2)); e_lo <- eval_at(lo)
    }

    if (e_lo$power >= target_power) { 
      stop("Target power is reached even at the minimum sample size of 2. Check your inputs (e.g., delta_power may be too large).") 
    }

  } else {                                     # search upward for an upper bound
    lo <- n_max_curr

    # double n_max
    hi <- 2 * n_max_curr; e_hi <- eval_at(hi)

    while (e_hi$power < target_power) {
      lo <- hi
      hi <- 2 * hi
      if (hi > n2_max_cap) stop("power target not reached below n2_max_cap = ", n2_max_cap)
      e_hi <- eval_at(hi)
    }
  }



  # ---- integer bisection
  # Now that we have lower and upper bound for n2max
  while (hi - lo > 1) {
    mid <- floor((lo + hi) / 2)
    e_mid <- eval_at(mid)
    if (e_mid$power >= target_power) { 
      hi <- mid; e_hi <- e_mid 
      } else lo <- mid
  }

  out <- list(n1_seq = e_hi$ns$n1, n2_seq = e_hi$ns$n2,
              p = e_hi$cal$p,
              achieved_t1e   = e_hi$cal$achieved_t1e,
              achieved_power = e_hi$power,
              target_t1e = target_t1e, target_power = target_power,
              delta_power = delta_power, criterion = criterion,
              decisions = e_hi$cal$decisions)
  out
}

