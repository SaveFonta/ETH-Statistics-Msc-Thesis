library(RBesT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(checkmate)  
library(assertthat) 



#' Average Bayesian OC where theta1 = theta2 + delta
#'
#' Computes E_{theta2 ~ design_prior2}[ OC(theta2 + delta, theta2) ]
#'
#' @param prior1   normMix prior for arm 1
#' @param prior2   normMix prior for arm 2
#' @param n1       sample size arm 1 (treatment)
#' @param n2       sample size arm 2 (placebo)
#' @param decision a decision2S object
#' @param delta    numeric scalar: theta1 = theta2 + delta
#' @param design_prior2     normMix distribution that acts ad design distribution of theta2
#' @param sigma1   known SD for arm 1 (defaults to prior reference scale)
#' @param sigma2   known SD for arm 2 (defaults to prior reference scale)
#' @param ... everything that goes inside oc2Sfun
#'
#' @return numeric scalar: the averaged OC
avgoc2S <- function(
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








































# decisions should be passed as a list the length the number of stages
# so let's say in a two stage with futility we have 
# succ.early <- decision2S(0.95, 0, lower.tail = FALSE)
# fut.early <- decision2S(0.20, 0, lower.tail = TRUE)
# succ.final <- decision2S (0.9, 0, lower.tail = FALSE)

# decision <- list ( list (success = succ.early, futility = fut.early), success = succ_final)


# n1 and n2 are vector of CUMULATIVE sample sizes for each stages (c(interim, final))
# this mean NOT c(50,50), but c(50,100)

oc2S_seq.dual.normMix <- function(
  prior1, prior2,
  n1, n2,
  decisions,
  sigma1, sigma2,
  eps = 1e-6, Ngrid = 10
) {

  ## 1. SETUP
  if (missing(sigma1)) sigma1 <- RBesT::sigma(prior1)
  if (missing(sigma2)) sigma2 <- RBesT::sigma(prior2)

  if (length(n1) != 2 || length(n2) != 2) {
    stop("This implementation currently supports exactly 2 stages: provide length-2 cumulative n1 and n2.")
  }
  if (length(decisions) != 2) {
    stop("This implementation currently supports exactly 2 stages: provide 2 stage decision entries.")
  }
  if (any(diff(n1) <= 0) || any(diff(n2) <= 0)) {
    stop("n1 and n2 must be strictly increasing cumulative sample sizes.")
  }

  RBesT::sigma(prior1) <- sigma1
  RBesT::sigma(prior2) <- sigma2


  ## Check decision list and standardise it to list(success=, futility=)
  check_decision <- function(d, k) {

    # if it is a single decision object, it is only a success rule
    if (methods::is(d, "decision2S")) {
      return(list(success = d, futility = NULL))
    }

    if (!is.list(d)) {
      stop(paste0("Error at Stage ", k, ": decision must be a 'decision2S' object or a list"))
    }

    if (is.null(d$success)) {
      stop(paste0("Error at Stage ", k, ": the list is missing the mandatory 'success' rule."))
    }

    if (!methods::is(d$success, "decision2S")) {
      stop(paste0("Error at Stage ", k, ": The 'success' rule must be a valid RBesT decision2S object."))
    }

    return(list(success = d$success, futility = d$futility))
  }

  dec <- lapply(seq_along(decisions), function(k) check_decision(decisions[[k]], k))


  ## create the decision boundaries:
  ## given a stage k and a decision object, returns the boundary function
  ## that maps y2 -> critical y1 value. Returns NULL if no rule exists.
  make_bnd <- function(k, d_obj) {
    if (is.null(d_obj)) return(NULL)
    RBesT::decision2S_boundary(
      prior1, prior2, n1[k], n2[k], d_obj,
      sigma1, sigma2, eps, Ngrid
    )
  }

  bnd_succ <- lapply(1:2, function(k) make_bnd(k, dec[[k]]$success))
  bnd_fut  <- lapply(1:2, function(k) make_bnd(k, dec[[k]]$futility))

  # Now we have:
  # bnd_succ[[1]](y2) = critical y1 for early success at interim given placebo mean y2
  # bnd_fut[[1]](y2)  = critical y1 for futility stop at interim given placebo mean y2
  # bnd_succ[[2]](y2) = critical y1 for final success given placebo mean y2


  ## Extract the tail direction from each decision object.
  ## This encodes whether "success if y1 > boundary" (lower.tail = FALSE, higher is better)
  ## or "success if y1 < boundary" (lower.tail = TRUE, lower is better).

  lower.tail_succ <- lapply(1:2, function(k) attr(dec[[k]]$success,  "lower.tail"))
  lower.tail_fut  <- lapply(1:2, function(k) {
    if (is.null(dec[[k]]$futility)) NULL else attr(dec[[k]]$futility, "lower.tail")
  })

  if ((!is.null(lower.tail_fut[[1]]) && lower.tail_succ[[1]] == lower.tail_fut[[1]]) |
      (!is.null(lower.tail_fut[[2]]) && lower.tail_succ[[2]] == lower.tail_fut[[2]])) {
    stop("The success criteria and futility criteria must be on opposite tails, 
                                                        ensure to control it using 'lower.tail' in 'decision2S'")
  } 

  ## Marginal standard errors at each stage
  ## (used for integration limits and sampling distributions)
  mSEM1 <- sigma1 / sqrt(n1)   # c(sem1_stage1, sem1_stage2)
  mSEM2 <- sigma2 / sqrt(n2)

  ## Conditional standard errors of the incremental data between stages
  ## These answer: given interim data y_1p, what is the SD of the final mean y_1?
  ## Derived from: y_final = (n_interim * y_interim + n_increment * y_new) / n_final
  cSEM1 <- c(sigma1 / sqrt(n1[1]), sigma1 * sqrt(n1[2] - n1[1]) / n1[2])
  cSEM2 <- c(sigma2 / sqrt(n2[1]), sigma2 * sqrt(n2[2] - n2[1]) / n2[2])


  # small note -> cSEM1[1] == mSEM1[1]



  ## Conditional means of the final stage statistic given interim data
  ## These answer: if my interim mean is y_p and the true param is theta,
  ## what is E[y_final | y_interim = y_p, theta]?
  cmean1 <- function(y1p, theta1) (n1[1] * y1p + (n1[2] - n1[1]) * theta1) / n1[2]
  cmean2 <- function(y2p, theta2) (n2[1] * y2p + (n2[2] - n2[1]) * theta2) / n2[2]



  ## Returns OC for a specific (theta1, theta2) pair
  design_fun <- function(theta1, theta2) {

    ## --- CACHE  ---
    for (k in 1:2) {
      lim2_k <- qnorm(c(eps / 2, 1 - eps / 2), theta2, mSEM2[k])
      lim1_k <- qnorm(c(eps / 2, 1 - eps / 2), theta1, mSEM1[k])
      if (!is.null(bnd_succ[[k]])) bnd_succ[[k]](lim2_k, lim1 = lim1_k)
      if (!is.null(bnd_fut[[k]]))  bnd_fut[[k]](lim2_k,  lim1 = lim1_k)
    }


    ## --- STAGE 1 (Interim) ---
    ## Marginal sampling distribution of the interim placebo mean
    ## y2_1 = N(theta2, mSEM2[1])
    mix_y2_1 <- RBesT::mixnorm(c(1, theta2, mSEM2[1]), sigma = mSEM2[1])

    ## P(early success at interim) = E_{y2_1}[ P(y1_1 > (or <) D_1(y2_1) | theta1) ]
    # where D_1(y2_1) = bnd_succ[[1]](y2_1)
    p_succ_1 <- RBesT:::integrate_density_log(
        log_integrand = function(x) {
          pnorm(bnd_succ[[1]](x), theta1, cSEM1[1], 
          lower.tail = lower.tail_succ[[1]], log.p = TRUE)
        },
      mix     = mix_y2_1,
      Lplower = RBesT:::logit(eps / 2),
      Lpupper = RBesT:::logit(1 - eps / 2)
    )

    ## P(futility stop at interim) = E_{y2_1}[ P(y1_1 < bnd_fut[[1]](y2_1) | theta1) ]
    p_fut_1 <- 0
    if (!is.null(bnd_fut[[1]])) {
      p_fut_1 <- RBesT:::integrate_density_log(
        log_integrand = function(x) {
            pnorm(
              bnd_fut[[1]](x), theta1, cSEM1[1],
              lower.tail = lower.tail_fut[[1]], log.p = TRUE
            )
          }
        ,
        mix     = mix_y2_1,
        Lplower = RBesT:::logit(eps / 2),
        Lpupper = RBesT:::logit(1 - eps / 2)
      )
    }

    ## P(continue to stage 2) = 1 - P(early stop)
    p_cont_1 <- 1 - p_succ_1 - p_fut_1


    ## --- STAGE 2 (Final) ---
    ## P(success at final) = E_{y2_1}[ E_{y1_1 in continuation region C_1}[
    ##                          E_{y2_2 | y2_1}[ P(y1_2 > (or <) bnd_succ[[2]](y2_2) | y1_1, theta1) ] ] ]
    ##
    ## Outer integral over y2_1 ~ N(theta2, mSEM2[1]) 
    ## Middle integral over y1_1 in continuation region 
    ## Inner integral over y2_2 | y2_1 

   p_succ_2 <- RBesT:::integrate_density_log(
      log_integrand = function(y2_1_vec) {
        sapply(y2_1_vec, function(y2_1) {

          ## define continutation boundaries for y1_1 given y2_1.
          c_s <- bnd_succ[[1]](y2_1)
          
          if (!is.null(bnd_fut[[1]])) {
            c_f <- bnd_fut[[1]](y2_1)
          } else {
            c_f <- ifelse(lower.tail_succ[[1]] == FALSE, -Inf, Inf)
          }


          if (lower.tail_succ[[1]] == FALSE) {
            # Higher is better: Continue if c_f < y1 < c_s
            if (c_f >= c_s) return(-Inf) # No continuation region (returns log(0))
            c_lo <- c_f
            c_hi <- c_s
          } else {
            # Lower is better: Continue if c_s < y1 < c_f
            if (c_s >= c_f) return(-Inf) 
            c_lo <- c_s
            c_hi <- c_f
          }

          ## Define Mixture and probability limits
          mix_y1_1 <- RBesT::mixnorm(c(1, theta1, cSEM1[1]), sigma = cSEM1[1])

          p_lo  <- pnorm(c_lo, theta1, cSEM1[1])
          p_hi  <- pnorm(c_hi, theta1, cSEM1[1])
        
          ## Prevent exactly 0 or 1 
          p_lo <- max(min(p_lo, 1 - 1e-12), 1e-12) 
          p_hi <- max(min(p_hi, 1 - 1e-12), 1e-12)
          

          ## integrate_density wants logit scale as extreme
          Lp_lo <- RBesT::logit(p_lo)
          Lp_hi <- RBesT::logit(p_hi)

          ## Middle integral
          inner <- RBesT:::integrate_density_log(
            log_integrand = function(y1_1_vec) {
              sapply(y1_1_vec, function(y1_1) {
                
                mu1 <- cmean1(y1_1, theta1)
                mu2 <- cmean2(y2_1, theta2)
                mix_y2_2 <- RBesT::mixnorm(c(1, mu2, cSEM2[2]), sigma = cSEM2[2])

                ## Calculate probability of Stage 2 success
                succ_2 <- RBesT:::integrate_density_log(
                  log_integrand = function(y2_2_vec) {
                    sapply(y2_2_vec, function(y2_2) {
                      pnorm(
                        bnd_succ[[2]](y2_2), mu1, cSEM1[2],
                        lower.tail = lower.tail_succ[[2]], log.p = TRUE
                      )
                    })
                  },
                  mix     = mix_y2_2,
                  Lplower = RBesT::logit(eps / 2),
                  Lpupper = RBesT::logit(1 - eps / 2)
                )

                ## succ_2 is a natural probability. 
                ## We must log it before handing it up to the middle log_integrand
                if (succ_2 <= 0) return(-Inf)
                log(succ_2)
              })
            },
            mix     = mix_y1_1,
            Lplower = Lp_lo,
            Lpupper = Lp_hi
          )


          ## Return log for outer integrate_density_log
          ## The inner tryCatch returned a natural probability, so we log it here
          if (inner <= 0) return(-Inf)
          log(inner)
        })
      },
      mix     = mix_y2_1,
      Lplower = RBesT::logit(eps / 2),
      Lpupper = RBesT::logit(1 - eps / 2)
    )

    ## --- OPERATING CHARACTERISTICS ---
    p_total  <- p_succ_1 + p_succ_2   # overall power
    p_no_rej <- 1 - p_total          # overall non-rejection probability

    ## --- EXPECTED SAMPLE SIZES ---
    ## E[N] = n_interim * P(stop at interim) + n_final * P(reach final)
    EN1 <- n1[1] * (p_succ_1 + p_fut_1) + n1[2] * p_cont_1
    EN2 <- n2[1] * (p_succ_1 + p_fut_1) + n2[2] * p_cont_1

    ## E[N | +]
    EN1_rej <- if (p_total > 0) (n1[1] * p_succ_1 + n1[2] * p_succ_2) / p_total else NA_real_
    EN2_rej <- if (p_total > 0) (n2[1] * p_succ_1 + n2[2] * p_succ_2) / p_total else NA_real_

    ## E[N | -]
    ## p_cont_1 - p_succ_2 = P(reach final but fail)
    p_fail_2 <- max(0, p_cont_1 - p_succ_2)

    EN1_norej <- if (p_no_rej > 0) (n1[1] * p_fut_1 + n1[2] * (p_fail_2)) / p_no_rej else NA_real_
    EN2_norej <- if (p_no_rej > 0) (n2[1] * p_fut_1 + n2[2] * (p_fail_2)) / p_no_rej else NA_real_

    return(data.frame(
      Theta1           = theta1,
      Theta2           = theta2,
      P_Stop_Eff_Stg1  = p_succ_1,
      P_Stop_Fut_Stg1  = p_fut_1,
      P_Continue_Stg1  = p_cont_1,
      P_Success_Stg2   = p_succ_2,
      Total_Power      = p_total,
      EN_Trt           = EN1,
      EN_Pbo           = EN2,
      EN_Trt_Rej       = EN1_rej,
      EN_Pbo_Rej       = EN2_rej,
      EN_Trt_NoRej     = EN1_norej,
      EN_Pbo_NoRej     = EN2_norej
    ))
  }

  return(design_fun)
}











avgoc2S_seq.normMix <- function(
  prior1, prior2, n1, n2, decisions, delta, design_prior2, eps = 1e-6, Ngrid = 10, ...
) {
  
  # Creates OC
  oc_fun <- oc2S_seq.dual.normMix(prior1, prior2, n1, n2, decisions, eps = eps, ...)

  # output function
  design_fun <- function(delta_new = delta, design_prior2_new = design_prior2) {

    res <- RBesT:::integrate_density_log(
      log_integrand = function(x) {
        # Force scalar evaluation over the vector of grid points 'x'
        sapply(x, function(xi) {
          # Evaluate OC for the scalar xi
          res_df <- oc_fun(xi + delta_new, xi)
          # Extract the total power probability and take the log
          log(res_df$Total_Power)
        })
      },
      mix = design_prior2_new,
      Lplower = RBesT::logit(eps / 2),
      Lpupper = RBesT::logit(1 - eps / 2)
    )
    
    return(res)
  }
  
  return(design_fun)
}

















avgoc2S_seq2.normMix <- function(
  prior1, prior2, n1, n2, decisions, delta, design_prior2, eps = 1e-6, ...
) {
  
  # 1. Create the conditional OC function
  oc_fun <- oc2S_seq.dual.normMix(prior1, prior2, n1, n2, decisions, eps = eps, ...)
  
  # 2. Output function
  design_eval_fun <- function(delta_new = delta, design_prior2_new = design_prior2) {
    
    # Define integration limits based on the quantiles of the design prior
    lower_bnd <- RBesT::qmix(design_prior2_new, eps / 2)
    upper_bnd <- RBesT::qmix(design_prior2_new, 1 - eps / 2)
    
    # Helper function to compute the expected value of ANY specific conditional metric
    # It takes an 'extractor' function to pull the exact metric from oc_fun's output
    calc_expectation <- function(metric_extractor) {
      integrand <- function(x_vec) {
        sapply(x_vec, function(xi) {
          # Evaluate conditional OC for a given theta_c (xi) and theta_t (xi + delta_new)
          cond_res <- oc_fun(xi + delta_new, xi)
          
          # Extract the specific metric value (e.g., Power, EN_t, P_Succ)
          metric_val <- metric_extractor(cond_res)
          
          # Multiply by the design prior density p_D(theta_c)
          prior_dens <- RBesT::dmix(design_prior2_new, xi)
          
          return(metric_val * prior_dens)
        })
      }
      
      # Integrate using standard stats::integrate
      res <- stats::integrate(integrand, lower = lower_bnd, upper = upper_bnd, 
                              rel.tol = .Machine$double.eps^0.25)$value
      return(res)
    }
    
    # --- CALCULATE OVERALL METRICS ---
    # We pass an anonymous function to extract the exact column/metric we want to average
    avg_Power       <- calc_expectation(function(res) res$Total_Power)
    avg_Prob_Fut    <- calc_expectation(function(res) res$P_Stop_Fut_Stg1)
    avg_EN_t        <- calc_expectation(function(res) res$EN_Trt)
    avg_EN_c        <- calc_expectation(function(res) res$EN_Pbo)
    avg_EN_t_Succ   <- calc_expectation(function(res) res$EN_Trt_Rej)
    avg_EN_c_Succ   <- calc_expectation(function(res) res$EN_Pbo_Rej)
    avg_EN_t_Fail   <- calc_expectation(function(res) res$EN_Trt_NoRej)
    avg_EN_c_Fail   <- calc_expectation(function(res) res$EN_Pbo_NoRej)
    
    overall_res <- data.frame(
      Power = avg_Power,
      Prob_Fut_seq = avg_Prob_Fut,
      EN_t = avg_EN_t,
      EN_c = avg_EN_c,
      EN_t_Succ = avg_EN_t_Succ,
      EN_c_Succ = avg_EN_c_Succ,
      EN_t_Fail = avg_EN_t_Fail,
      EN_c_Fail = avg_EN_c_Fail
    )
    
    # --- CALCULATE PER-STAGE METRICS ---
    # Evaluate oc_fun once at the prior mean just to extract the structural template 
    # (number of stages, fixed N sizes)
    prior_mean <- RBesT::summary(design_prior2_new)["mean"]
    base_res <- oc_fun(prior_mean + delta_new, prior_mean)
    per_stage_res <- data.frame(
      Stage  = 1:2,
      N_Trt  = n1,
      N_Ctrl = n2,
      P_Succ = c(
        calc_expectation(function(res) res$P_Stop_Eff_Stg1),
        calc_expectation(function(res) res$P_Success_Stg2)
      ),
      P_Fut = c(
        calc_expectation(function(res) res$P_Stop_Fut_Stg1),
        0
      )
    )
    
    # Compute the average cumulative probabilities step-by-step
    per_stage_res$Cum_P_Succ <- cumsum(per_stage_res$P_Succ)
    per_stage_res$Cum_P_Fut  <- cumsum(per_stage_res$P_Fut)
    
    return(list(Overall = overall_res, Per_Stage = per_stage_res))
  }
  
  return(design_eval_fun)
}
















# MONTE CARLO RESULTS


oc2_seq_mc.normMix <- function(theta_1, theta_2, prior_1, prior_2, 
                                  n1_seq, n2_seq,  decisions_list,  # recall that n1_seq is cumulative
                                  sigma_val = 88, n_sim = 100000, seed = 123) {

  set.seed(seed)
  if (length(theta_1) == 1) theta_1 <- rep(theta_1, n_sim)
  if (length(theta_2) == 1) theta_2 <- rep(theta_2, n_sim)


  # check length 
  if(length(theta_1) != n_sim) stop ("Number of theta_1 provided should be either one for classic OC or equal to 
                                      n_sim for avg OC") 
  if(length(theta_2) != n_sim) stop ("Number of theta_2 provided should be either one for classic OC or equal to 
                                      n_sim for avg OC")



  K <- length(n2_seq)
 
  if(length(n2_seq) != length(n1_seq)) stop ("length of n1_seq and n2_seq must coincide")



  active_idx <- 1:n_sim # Keep track of trials that haven't stopped
  
  # Initialize a lot of tracking vectors, all of length n_sim (each element can be seen as an evolving trial) 
  stop_eff <- rep(FALSE, n_sim)
  stop_fut <- rep(FALSE, n_sim)
  stage_stopped <- rep(K, n_sim) # Default to final stage
  
  y_2_curr <- rep(0, n_sim)
  y_1_curr <- rep(0, n_sim)
  n_2_prev <- 0
  n_1_prev <- 0

  # Vector of probabilities of stopping or futility for each stage
  p_succ_vec <- rep(0, K)
  p_fut_vec  <- rep(0, K)

  
  for (k in 1:K) {
    if (length(active_idx) == 0) break 
    
    n_active <- length(active_idx)


    dn_2 <- n2_seq[k] - n_2_prev
    dn_1 <- n1_seq[k] - n_1_prev
    
    # Generate new data for the trials still active
    # to have a vector of length active_idx 
    y_1_inc <- rnorm(n_active, mean = theta_1[active_idx], sd = sigma_val / sqrt(dn_1))
    y_2_inc <- rnorm(n_active, mean = theta_2[active_idx], sd = sigma_val / sqrt(dn_2))


    # Calculate cumulative means
    # [mean till stage (k-1) * sample size (k-1) + mean of stage k * sample size stage k ] / total sample size till stage k 
    y_2_curr[active_idx] <- (n_2_prev * y_2_curr[active_idx] + dn_2 * y_2_inc) / n2_seq[k]
    y_1_curr[active_idx] <- (n_1_prev * y_1_curr[active_idx] + dn_1 * y_1_inc) / n1_seq[k]

    #extract decisions 
    d_succ <- decisions_list[[k]]$success
    d_fut  <- decisions_list[[k]]$futility
    
    # Evaluate Efficacy Boundary (if defined for this stage)
    if (!is.null(d_succ)) {
      bnd_eff_fun <- decision2S_boundary(prior_1, prior_2, n1_seq[k], n2_seq[k], d_succ)
      crit_val_eff <- bnd_eff_fun(y_2_curr[active_idx])

      is_lower_eff <- attr(d_succ, "lower.tail")

      is_succ <- if(is_lower_eff) y_1_curr[active_idx] <= crit_val_eff else y_1_curr[active_idx] > crit_val_eff
    } else {
      is_succ <- rep(FALSE, n_active)
    }
    
    # Evaluate Futility Boundary (if defined for this stage)
    if (!is.null(d_fut)) {
      bnd_fut_fun <- decision2S_boundary(prior_1, prior_2, n1_seq[k], n2_seq[k], d_fut)
      crit_val_fut <- bnd_fut_fun(y_2_curr[active_idx])
      is_lower_fut <- attr(d_fut, "lower.tail")
      is_fut <- if(is_lower_fut) y_1_curr[active_idx] <= crit_val_fut else y_1_curr[active_idx] > crit_val_fut
    } else {
      is_fut <- rep(FALSE, n_active)
    }
    
    # Record stops, those are indeces of active_idx
    just_stopped_eff <- is_succ
    just_stopped_fut <- is_fut & !is_succ #shouldn't happen that both futility and efficacy are crosse, bu anyway
    just_stopped <- just_stopped_eff | just_stopped_fut
    
    # in the total vector of length n_sim, the ones that stopped now are full.vec[active_idx[stop.now]]
    stop_eff[active_idx[just_stopped_eff]] <- TRUE
    stop_fut[active_idx[just_stopped_fut]] <- TRUE
    stage_stopped[active_idx[just_stopped]] <- k    # at which stage does this simulated trial stopped? 
    
    # Record prob to stop 
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

  EN_1_succ <- if(length(idx_succ) > 0) mean(n1_seq[stage_stopped[idx_succ]]) else NA
  EN_2_succ <- if(length(idx_succ) > 0) mean(n2_seq[stage_stopped[idx_succ]]) else NA
  
  EN_1_fail <- if(length(idx_fail) > 0) mean(n1_seq[stage_stopped[idx_fail]]) else NA
  EN_2_fail <- if(length(idx_fail) > 0) mean(n2_seq[stage_stopped[idx_fail]]) else NA




  per_stage_df <- data.frame(
    Stage      = 1:K,
    N_Trt      = n1_seq,
    N_Ctrl     = n2_seq,
    P_Succ     = p_succ_vec,
    P_Fut      = p_fut_vec,
    Cum_P_Succ = cumsum(p_succ_vec),
    Cum_P_Fut  = cumsum(p_fut_vec)
  )
  
    # Return OC
  return(list(
    Overall = c( 
    Power = sum(stop_eff) / n_sim,
    Prob_Fut_seq = sum(stop_fut) / n_sim,
    EN_t = mean(n1_seq[stage_stopped]),
    EN_c = mean(n2_seq[stage_stopped]),
    EN_t_Succ = EN_1_succ,
    EN_c_Succ = EN_2_succ,
    EN_t_Fail = EN_1_fail,
    EN_c_Fail = EN_2_fail
  ), 
  Per_Stage = per_stage_df) )
}



avgoc2_seq_mc.normMix <- function(prior_1, prior_2, 
                                  n1_seq, n2_seq,  decisions_list,  # recall that n1_seq is cumulative
                                  delta = 0, design_prior_c, sigma_val = 88, n_sim = 100000, seed = 123) {
set.seed(seed)
theta_c_draws <- rmix(design_prior_c, n_sim)
theta_t_draws <- theta_c_draws + delta 


res_avgT1E <- oc2_seq_mc.normMix(
  theta_1 = theta_t_draws, 
  theta_2 = theta_c_draws, 
  prior_1 = prior_1, prior_2 = prior_2, 
  n1_seq = n1_seq, n2_seq = n2_seq, 
  decisions_list = decisions_list, n_sim = n_sim, seed = seed
)

return(res_avgT1E)

}













# 1. Define Priors and Sample Sizes
sigma_val <- 88
prior_t <- mixnorm(c(1, 0, 100), sigma = sigma_val)
prior_c <- mixnorm(c(1, 0, 100), sigma = sigma_val)
design_prior_c <- prior_c
n_seq <- c(50, 100) # Cumulative: 50 at interim, 100 at final

# 2. Define Decisions
succ_early <- decision2S(0.95, 0, lower.tail = FALSE)
fut_early  <- decision2S(0.20, 0, lower.tail = TRUE)
succ_final <- decision2S(0.975, 0, lower.tail = FALSE)

decisions <- list(
  list(success = succ_early, futility = fut_early),
  list(success = succ_final)
)

# 3. Analytic Classic OC (Targeting Pointwise Theta)
theta_t_null <- 0
theta_c_null <- 0

oc_analytic_fun <- oc2S_seq.dual.normMix(
  prior1 = prior_t, prior2 = prior_c,
  n1 = n_seq, n2 = n_seq,
  decisions = decisions,
  sigma1 = sigma_val, sigma2 = sigma_val
)
res_classic_analytic <- oc_analytic_fun(theta_t_null, theta_c_null)

# 4. Monte Carlo Classic OC
res_classic_mc <- oc2_seq_mc.normMix(
  theta_1 = theta_t_null, theta_2 = theta_c_null,
  prior_1 = prior_t, prior_2 = prior_c,
  n1_seq = n_seq, n2_seq = n_seq,
  decisions_list = decisions, sigma_val = sigma_val, n_sim = 1e5
)




# 5. Analytic Average Sequential Metrics
# This builds the evaluation function for our Bayesian sequential design
avg_seq_fun <- avgoc2S_seq2.normMix(
  prior1 = prior_t, 
  prior2 = prior_c,
  n1 = n_seq, 
  n2 = n_seq,
  decisions = decisions,
  delta = 0, # Delta null for calculating at1e_seq
  design_prior2 = prior_c,
  sigma1 = sigma_val, 
  sigma2 = sigma_val
)

# Execute the function to get the comprehensive list of Overall and Per_Stage metrics
res_avg_analytic <- avg_seq_fun()

# Print the analytic results to verify the new structure
print("--- Analytic Average Sequential Metrics ---")
print(res_avg_analytic)


# 6. Monte Carlo Average OC
# (Make sure your avgoc2_seq_mc.normMix function is also updated to return 
# the expected sample sizes and per-stage stopping probabilities so you can 
# perform a direct 1:1 validation against the analytic results.)
res_avg_mc <- avgoc2_seq_mc.normMix(
  prior_1 = prior_t, 
  prior_2 = prior_c,
  n1_seq = n_seq, 
  n2_seq = n_seq,
  decisions_list = decisions,
  delta = 0, 
  design_prior_c = design_prior_c,
  sigma_val = sigma_val, 
  n_sim = 1e7, 
  seed = 11
)

# Print the MC results for validation
print("--- Monte Carlo Average Sequential Metrics ---")
print(res_avg_mc)


saveRDS(list(res_classic_analytic, res_classic_mc, res_avg_analytic,res_avg_mc), file = "data/sequential_at1e_results.rds")


