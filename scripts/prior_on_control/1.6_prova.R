##############################################################################
# Bayesian Sequential Clinical Trial – Monte Carlo Simulation
# Strictly separates: Reality (DGP / Design Prior) vs. Belief (Analysis Model)
##############################################################################

simulate_bayesian_trial <- function(
    n_sims    = 10000L,
    n_looks   = c(20, 40, 60, 80, 100),  # cumulative N per arm at each look

    # ---- Reality (DGP) -------------------------------------------------------
    mu_d      = 0,      # design prior mean for true control response
    sigma_d   = 0,      # design prior SD  (0 = fixed true control mean)
    delta_true = 2,     # true treatment effect  Δ_true
    sigma_c   = 5,      # true patient-level SD, control arm
    sigma_t   = 5,      # true patient-level SD, treatment arm

    # ---- Belief (Analysis Model) --------------------------------------------
    mu_a      = 0,      # control historical-data prior mean
    sigma_a   = 10,     # control historical-data prior SD
    mu_0      = 0,      # treatment-effect prior mean
    tau       = 10,     # treatment-effect prior SD

    # ---- Decision Boundaries -------------------------------------------------
    # Efficacy: stop for success  if P(Δ > q_eff[k] | data) > p_eff[k]
    q_eff     = c(2, 2, 2, 2, 1.5),
    p_eff     = c(0.99, 0.99, 0.975, 0.975, 0.95),

    # Futility: stop for futility if P(Δ < q_fut[k] | data) > p_fut[k]
    #           NA at a look means no futility check at that look
    q_fut     = c(0.5, 0.5, 1.0, 1.0, NA),
    p_fut     = c(0.90, 0.90, 0.95, 0.95, NA)
) {

  # --------------------------------------------------------------------------
  # 0. Input validation
  # --------------------------------------------------------------------------
  n_stages <- length(n_looks)
  stopifnot(
    length(q_eff) == n_stages, length(p_eff) == n_stages,
    length(q_fut) == n_stages, length(p_fut) == n_stages,
    all(diff(n_looks) > 0),
    sigma_d  >= 0, sigma_c > 0, sigma_t > 0,
    sigma_a  > 0,  tau     > 0
  )

  # Pre-compute squared SDs (used repeatedly in inner loop)
  sigma_a2 <- sigma_a^2
  tau2      <- tau^2
  sigma_c2  <- sigma_c^2
  sigma_t2  <- sigma_t^2

  # --------------------------------------------------------------------------
  # 1. Pre-allocate result storage
  # --------------------------------------------------------------------------
  outcomes      <- integer(n_sims)           # 1=Success, 2=Futility, 3=Inconclusive
  stopping_N    <- integer(n_sims)           # total N (both arms) at stopping
  stopping_k    <- integer(n_sims)           # stage index at stopping

  # Named outcome codes
  SUCCESS     <- 1L
  FUTILITY    <- 2L
  INCONCLUSIVE <- 3L

  # --------------------------------------------------------------------------
  # 2. Main simulation loop
  # --------------------------------------------------------------------------
  for (i in seq_len(n_sims)) {

    # --- Reality: draw true means for this run ---
    theta_c_true <- if (sigma_d > 0) rnorm(1L, mu_d, sigma_d) else mu_d
    theta_t_true <- theta_c_true + delta_true

    outcome_i <- INCONCLUSIVE
    k_stop    <- n_stages
    N_stop    <- 2L * n_looks[n_stages]

    # --- Inner loop over stages ---
    for (k in seq_len(n_stages)) {

      n_ck <- n_looks[k]          # cumulative control N at stage k
      n_tk <- n_looks[k]          # cumulative treatment N at stage k (same schedule)

      # ---- Reality: generate observed cumulative sample means ----
      y_ck <- rnorm(1L, theta_c_true, sqrt(sigma_c2 / n_ck))
      y_tk <- rnorm(1L, theta_t_true, sqrt(sigma_t2 / n_tk))

      # ---- Belief: Bayesian update for control mean ----
      sigma_c2_n  <- sigma_c2 / n_ck          # sampling variance of y_ck
      V_pool_k    <- (sigma_a2 * sigma_c2_n) /
                     (sigma_a2 + sigma_c2_n)   # posterior variance of θ_c
      y_pool_k    <- (mu_a * sigma_c2_n + y_ck * sigma_a2) /
                     (sigma_a2 + sigma_c2_n)   # posterior mean  of θ_c

      # ---- Belief: data likelihood for effect Δ ----
      delta_hat_k  <- y_tk - y_pool_k
      V_data_k     <- (sigma_t2 / n_tk) + V_pool_k

      # ---- Belief: final posterior for Δ ----
      V_post_k  <- 1.0 / (1.0 / tau2 + 1.0 / V_data_k)
      mu_post_k <- V_post_k * (mu_0 / tau2 + delta_hat_k / V_data_k)
      sd_post_k <- sqrt(V_post_k)

      # ---- Decision logic ----
      # Efficacy: P(Δ > q_eff[k] | data)
      P_eff <- pnorm(q_eff[k], mean = mu_post_k, sd = sd_post_k,
                     lower.tail = FALSE)

      if (P_eff > p_eff[k]) {
        outcome_i <- SUCCESS
        k_stop    <- k
        N_stop    <- n_ck + n_tk
        break
      }

      # Futility: P(Δ < q_fut[k] | data)  — skip if NA at this look
      if (!is.na(q_fut[k]) && !is.na(p_fut[k])) {
        P_fut <- pnorm(q_fut[k], mean = mu_post_k, sd = sd_post_k,
                       lower.tail = TRUE)

        if (P_fut > p_fut[k]) {
          outcome_i <- FUTILITY
          k_stop    <- k
          N_stop    <- n_ck + n_tk
          break
        }
      }

    }  # end inner stage loop

    outcomes[i]   <- outcome_i
    stopping_k[i] <- k_stop
    stopping_N[i] <- N_stop

  }  # end simulation loop

  # --------------------------------------------------------------------------
  # 3. Aggregate results
  # --------------------------------------------------------------------------
  is_success  <- outcomes == SUCCESS
  is_futility <- outcomes == FUTILITY

  n_success  <- sum(is_success)
  n_futility <- sum(is_futility)

  # Per-stage stopping probabilities
  Prob_Stop_Stage     <- tabulate(stopping_k, nbins = n_stages) / n_sims
  Prob_Success_Stage  <- tabulate(stopping_k[is_success],  nbins = n_stages) / n_sims
  Prob_Futility_Stage <- tabulate(stopping_k[is_futility], nbins = n_stages) / n_sims

  # Expected sample sizes
  ESS_Total    <- mean(stopping_N)
  ESS_Success  <- if (n_success  > 0) mean(stopping_N[is_success])  else NA_real_
  ESS_Futility <- if (n_futility > 0) mean(stopping_N[is_futility]) else NA_real_

  # --------------------------------------------------------------------------
  # 4. Return named list
  # --------------------------------------------------------------------------
  list(
    # ---- Headline metrics ----
    Power_or_Type1      = n_success / n_sims,
    Prob_Futility       = n_futility / n_sims,
    Prob_Inconclusive   = sum(outcomes == INCONCLUSIVE) / n_sims,

    # ---- Expected sample sizes ----
    ESS_Total           = ESS_Total,
    ESS_Success         = ESS_Success,
    ESS_Futility        = ESS_Futility,

    # ---- Per-stage stopping probabilities ----
    Prob_Stop_Stage     = setNames(Prob_Stop_Stage,     paste0("Stage_", seq_len(n_stages))),
    Prob_Success_Stage  = setNames(Prob_Success_Stage,  paste0("Stage_", seq_len(n_stages))),
    Prob_Futility_Stage = setNames(Prob_Futility_Stage, paste0("Stage_", seq_len(n_stages))),

    # ---- Raw (for diagnostics / further analysis) ----
    raw = data.frame(
      outcome    = factor(outcomes, levels = 1:3,
                          labels = c("Success", "Futility", "Inconclusive")),
      stage      = stopping_k,
      total_N    = stopping_N
    )
  )
}


##############################################################################
# Example usage
##############################################################################

# Uncomment to run:
#
# set.seed(42)
 res <- simulate_bayesian_trial(
   n_sims     = 10000,
   n_looks    = c(20, 40, 60, 80, 100),
   mu_d       = 0,   sigma_d    = 0,
   delta_true = 2,
   sigma_c    = 5,   sigma_t    = 5,
   mu_a       = 0,   sigma_a    = 10,
   mu_0       = 0,   tau        = 10,
   q_eff  = c(2,   2,   2,   2,   1.5),
   p_eff  = c(0.99,0.99,0.975,0.975,0.95),
   q_fut  = c(0.5, 0.5, 1.0, 1.0, NA),
   p_fut  = c(0.90,0.90,0.95,0.95,NA)
 )
#Print headline summary
cat("Power / Type-I error :", round(res$Power_or_Type1, 3), "\n")
cat("Prob Futility        :", round(res$Prob_Futility,   3), "\n")
 cat("Prob Inconclusive    :", round(res$Prob_Inconclusive,3),"\n")
 cat("ESS (all)            :", round(res$ESS_Total,       1), "\n")
 cat("ESS (success)        :", round(res$ESS_Success,     1), "\n")
 cat("ESS (futility)       :", round(res$ESS_Futility,    1), "\n\n")
 print(round(rbind(
   Stop    = res$Prob_Stop_Stage,
   Success = res$Prob_Success_Stage,
   Futility= res$Prob_Futility_Stage
 ), 3))




































































 ##############################################################################
# Bayesian Sequential Clinical Trial – Merged Monte Carlo Simulation
#
# Combines the best of both implementations:
#   ✔ Incremental data generation (correct sequential dependency)
#   ✔ RBesT mixture priors for both arms
#   ✔ Asymmetric allocation (n1_seq ≠ n2_seq)
#   ✔ lower.tail directionality via decision object attributes
#   ✔ Multiple delta support via avgoc wrapper
#   ✔ Active-index vectorisation (fast)
#   ✔ MCE reporting on all headline metrics
#   ✔ Optional informative-weight tracking
##############################################################################


# ============================================================================
# CORE FUNCTION
# ============================================================================

oc2_seq_mc.normMix <- function(
    theta_1,           # true treatment means: scalar or length-n_sim vector
    theta_2,           # true control   means: scalar or length-n_sim vector
    prior_1,           # RBesT normMix prior for treatment arm
    prior_2,           # RBesT normMix prior for control   arm
    n1_seq,            # cumulative treatment sample sizes at each look
    n2_seq,            # cumulative control   sample sizes at each look
    decisions_list,    # list of K elements, each list(success = <dec2S | NULL>,
                       #                                futility = <dec2S | NULL>)
    sigma_1  = 88,     # true patient-level SD, treatment arm
    sigma_2  = 88,     # true patient-level SD, control   arm
    n_sim    = 100000L,
    seed     = NULL,
    # ---- optional informative-weight tracking --------------------------------
    weight.track       = FALSE,
    type.weight.track  = "uninformative.mixture",  # or "robustify"
    n_info_comps       = 3L                        # number of informative components
) {

  # --------------------------------------------------------------------------
  # 0. Setup & validation
  # --------------------------------------------------------------------------
  if (!is.null(seed)) set.seed(seed)

  if (length(theta_1) == 1L) theta_1 <- rep(theta_1, n_sim)
  if (length(theta_2) == 1L) theta_2 <- rep(theta_2, n_sim)

  if (length(theta_1) != n_sim)
    stop("theta_1 must be length 1 (classic OC) or length n_sim (average OC).")
  if (length(theta_2) != n_sim)
    stop("theta_2 must be length 1 (classic OC) or length n_sim (average OC).")

  K <- length(n2_seq)
  if (length(n1_seq)        != K) stop("n1_seq and n2_seq must have the same length.")
  if (length(decisions_list) != K) stop("decisions_list must have the same length as n2_seq.")
  if (any(diff(n1_seq) < 0) || any(diff(n2_seq) < 0))
    stop("n1_seq and n2_seq must be non-decreasing (cumulative).")

  # --------------------------------------------------------------------------
  # 1. Pre-allocate tracking vectors
  # --------------------------------------------------------------------------
  stop_eff      <- rep(FALSE, n_sim)   # did this trial stop for efficacy?
  stop_fut      <- rep(FALSE, n_sim)   # did this trial stop for futility?
  stage_stopped <- rep(K,     n_sim)   # at which stage did it stop? (default = final)

  # Running cumulative means (updated incrementally — key correctness fix)
  y1_curr <- rep(0, n_sim)
  y2_curr <- rep(0, n_sim)

  n1_prev <- 0L
  n2_prev <- 0L

  # Per-stage marginal stopping probabilities
  p_succ_vec <- rep(0, K)
  p_fut_vec  <- rep(0, K)

  # Optional weight trajectory
  w_inf_traj <- rep(NA_real_, K)

  active_idx <- seq_len(n_sim)   # shrinks as trials stop

  # --------------------------------------------------------------------------
  # 2. Sequential stage loop
  # --------------------------------------------------------------------------
  for (k in seq_len(K)) {

    if (length(active_idx) == 0L) break

    n_active <- length(active_idx)

    # ---- Incremental sample sizes at this stage ----------------------------
    dn1 <- n1_seq[k] - n1_prev   # new treatment patients added
    dn2 <- n2_seq[k] - n2_prev   # new control   patients added

    # ---- Reality: generate ONLY the new patients' mean (incremental) -------
    #   y_inc ~ N(theta, sigma^2 / dn)  — sufficient statistic for new batch
    y1_inc <- if (dn1 > 0L)
      rnorm(n_active, mean = theta_1[active_idx], sd = sigma_1 / sqrt(dn1))
    else 0

    y2_inc <- if (dn2 > 0L)
      rnorm(n_active, mean = theta_2[active_idx], sd = sigma_2 / sqrt(dn2))
    else 0

    # ---- Update cumulative means via weighted average ----------------------
    #   y_cum_k = [y_cum_(k-1) * n_(k-1)  +  y_inc_k * dn_k] / n_k
    y1_curr[active_idx] <-
      (n1_prev * y1_curr[active_idx] + dn1 * y1_inc) / n1_seq[k]

    y2_curr[active_idx] <-
      (n2_prev * y2_curr[active_idx] + dn2 * y2_inc) / n2_seq[k]

    # ---- Optional: track informative weight --------------------------------
    if (weight.track) {
      avg_y2 <- mean(y2_curr[active_idx])

      post_avg <- suppressMessages(
        RBesT::postmix(prior_2, m = avg_y2, n = n2_seq[k],
                       se = sigma_2 / sqrt(n2_seq[k]))
      )

      w_inf_traj[k] <- switch(
        type.weight.track,
        "robustify" = {
          # robust (non-informative) component is always the last row from robustify
          idx_robust <- nrow(post_avg)
          1 - post_avg["w", idx_robust]
        },
        "uninformative.mixture" = {
          sum(post_avg["w", seq_len(n_info_comps)])
        },
        stop("type.weight.track must be 'robustify' or 'uninformative.mixture'.")
      )
    }

    # ---- Belief: evaluate decision boundaries via RBesT --------------------
    d_succ <- decisions_list[[k]]$success
    d_fut  <- decisions_list[[k]]$futility

    # -- Efficacy ------------------------------------------------------------
    if (!is.null(d_succ)) {
      # decision2S_boundary returns a function: y2 -> critical value on y1
      bnd_eff <- suppressMessages(
        RBesT::decision2S_boundary(prior_1, prior_2,
                                   n1_seq[k], n2_seq[k], d_succ)
      )
      crit_eff     <- bnd_eff(y2_curr[active_idx])
      lower_eff    <- attr(d_succ, "lower.tail")   # directionality from decision object
      is_succ <- if (lower_eff)
        y1_curr[active_idx] <= crit_eff   # e.g. lower response is better
      else
        y1_curr[active_idx] >  crit_eff   # e.g. higher response is better
    } else {
      is_succ <- rep(FALSE, n_active)
    }

    # -- Futility ------------------------------------------------------------
    if (!is.null(d_fut)) {
      bnd_fut <- suppressMessages(
        RBesT::decision2S_boundary(prior_1, prior_2,
                                   n1_seq[k], n2_seq[k], d_fut)
      )
      crit_fut  <- bnd_fut(y2_curr[active_idx])
      lower_fut <- attr(d_fut, "lower.tail")
      is_fut <- if (lower_fut)
        y1_curr[active_idx] <= crit_fut
      else
        y1_curr[active_idx] >  crit_fut
    } else {
      is_fut <- rep(FALSE, n_active)
    }

    # -- Resolve conflicts: efficacy takes priority if both crossed ----------
    just_eff <- is_succ
    just_fut <- is_fut & !is_succ
    just_any <- just_eff | just_fut

    # -- Record outcomes in full n_sim-length vectors -----------------------
    stop_eff[active_idx[just_eff]]  <- TRUE
    stop_fut[active_idx[just_fut]]  <- TRUE
    stage_stopped[active_idx[just_any]] <- k

    # -- Marginal per-stage probabilities (denominator = n_sim throughout) --
    p_succ_vec[k] <- sum(just_eff) / n_sim
    p_fut_vec[k]  <- sum(just_fut) / n_sim

    # -- Shrink active set ---------------------------------------------------
    active_idx <- active_idx[!just_any]
    n1_prev    <- n1_seq[k]
    n2_prev    <- n2_seq[k]

  }  # end stage loop

  # --------------------------------------------------------------------------
  # 3. Aggregate results
  # --------------------------------------------------------------------------
  idx_succ <- which(stop_eff)
  idx_fail <- which(!stop_eff)   # futility + inconclusive

  EN_1_succ <- if (length(idx_succ) > 0) mean(n1_seq[stage_stopped[idx_succ]]) else NA_real_
  EN_2_succ <- if (length(idx_succ) > 0) mean(n2_seq[stage_stopped[idx_succ]]) else NA_real_
  EN_1_fail <- if (length(idx_fail) > 0) mean(n1_seq[stage_stopped[idx_fail]]) else NA_real_
  EN_2_fail <- if (length(idx_fail) > 0) mean(n2_seq[stage_stopped[idx_fail]]) else NA_real_

  power_est <- sum(stop_eff) / n_sim
  fut_est   <- sum(stop_fut) / n_sim

  # Monte Carlo standard errors
  mce_power <- sqrt(power_est * (1 - power_est) / n_sim)
  mce_fut   <- sqrt(fut_est   * (1 - fut_est)   / n_sim)
  mce_en1   <- sd(n1_seq[stage_stopped]) / sqrt(n_sim)
  mce_en2   <- sd(n2_seq[stage_stopped]) / sqrt(n_sim)

  per_stage_df <- data.frame(
    Stage           = seq_len(K),
    N_Trt           = n1_seq,
    N_Ctrl          = n2_seq,
    P_Succ          = p_succ_vec,
    P_Fut           = p_fut_vec,
    Cum_P_Succ      = cumsum(p_succ_vec),
    Cum_P_Fut       = cumsum(p_fut_vec),
    Exp_info_weight = w_inf_traj
  )

  list(
    Overall = c(
      Power        = power_est,
      MCE_Power    = mce_power,
      Prob_Fut_seq = fut_est,
      MCE_Fut      = mce_fut,
      EN_t         = mean(n1_seq[stage_stopped]),
      MCE_EN_t     = mce_en1,
      EN_c         = mean(n2_seq[stage_stopped]),
      MCE_EN_c     = mce_en2,
      EN_t_Succ    = EN_1_succ,
      EN_c_Succ    = EN_2_succ,
      EN_t_Fail    = EN_1_fail,
      EN_c_Fail    = EN_2_fail
    ),
    Per_Stage = per_stage_df
  )
}


# ============================================================================
# AVERAGE OC WRAPPER  (design prior over theta_c, vectorised over delta)
# ============================================================================

avgoc2_seq_mc.normMix <- function(
    prior_1,           # RBesT normMix analysis prior, treatment arm
    prior_2,           # RBesT normMix analysis prior, control   arm
    n1_seq,            # cumulative treatment sample sizes
    n2_seq,            # cumulative control   sample sizes
    decisions_list,    # list of K decision objects (see core function)
    delta          = 0,          # scalar or vector of true treatment effects
    design_prior_c,              # RBesT normMix design prior for true theta_c
    sigma_1        = 88,
    sigma_2        = 88,
    n_sim          = 100000L,
    seed           = 123L,
    # ---- weight tracking forwarded to core ---------------------------------
    weight.track      = FALSE,
    type.weight.track = "uninformative.mixture",
    n_info_comps      = 3L
) {

  set.seed(seed)

  # Directionality: read once from the first non-null success decision
  first_succ <- Filter(Negate(is.null), lapply(decisions_list, `[[`, "success"))[[1]]
  is_lower   <- attr(first_succ, "lower.tail")

  # Draw true control means once; reused across all delta values
  theta_c_draws <- RBesT::rmix(design_prior_c, n_sim)

  # Helper: run core simulator for a single delta value
  run_one <- function(d) {
    theta_t_draws <- if (is_lower) theta_c_draws - d else theta_c_draws + d

    oc2_seq_mc.normMix(
      theta_1           = theta_t_draws,
      theta_2           = theta_c_draws,
      prior_1           = prior_1,
      prior_2           = prior_2,
      n1_seq            = n1_seq,
      n2_seq            = n2_seq,
      decisions_list    = decisions_list,
      sigma_1           = sigma_1,
      sigma_2           = sigma_2,
      n_sim             = n_sim,
      seed              = NULL,          # seed already set above
      weight.track      = weight.track,
      type.weight.track = type.weight.track,
      n_info_comps      = n_info_comps
    )
  }

  if (length(delta) == 1L) {
    run_one(delta)
  } else {
    res <- lapply(delta, run_one)
    names(res) <- paste0("delta.", delta)
    res
  }
}


##############################################################################
# Example usage (uncomment to run — requires RBesT)
##############################################################################

# library(RBesT)
#
# --- Priors (single-normal for simplicity; use mixnorm() for mixtures) ---

#
sigma <- 88
prior_trt <- mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")
prior_ctrl <- mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")

# # --- Decision criteria (higher response = better) ---
 dec_eff <- decision2S(pc = 0.975, qc = 0, lower.tail = FALSE)
 dec_fut <- decision2S(pc = 0.10,  qc = 0, lower.tail = FALSE)

 decisions <- lapply(1:5, function(k) {
   list(success = if (k >= 2) dec_eff else NULL,
        futility = if (k <= 4) dec_fut else NULL)
 })
#
# # --- Classic OC (fixed true means) ---
 set.seed(42)
 res_classic <- oc2_seq_mc.normMix(
   theta_1        = 10,   # true treatment mean
   theta_2        = 0,    # true control mean
   prior_1        = prior_trt,
   prior_2        = prior_ctrl,
   n1_seq         = c(20, 40, 60, 80, 100),
   n2_seq         = c(10, 20, 30, 40,  50),   # 2:1 allocation
   decisions_list = decisions,
   sigma_1        = 88,
   sigma_2        = 88,
   n_sim          = 100
 )
 print(round(res_classic$Overall, 4))
 print(res_classic$Per_Stage)
#
# # --- Average OC over multiple deltas ---
# design_prior <- mixnorm(c(1, 0, 30), param = "mn")
#
# res_avg <- avgoc2_seq_mc.normMix(
#   prior_1        = prior_trt,
#   prior_2        = prior_ctrl,
#   n1_seq         = c(20, 40, 60, 80, 100),
#   n2_seq         = c(10, 20, 30, 40,  50),
#   decisions_list = decisions,
#   delta          = c(0, 5, 10, 15),
#   design_prior_c = design_prior,
#   n_sim          = 10000,
#   seed           = 42
# )
# lapply(res_avg, function(r) round(r$Overall, 4))