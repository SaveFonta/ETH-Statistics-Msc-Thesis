library(RBesT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(checkmate)  
library(assertthat) 



# -----------------------------------------------
# EXACT METHODS
# -----------------------------------------------



avgoc2S.normMix <- function(
  prior1, prior2, n1, n2, decision, delta, design_prior2,   eps = 1e-6, Ngrid = 10
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








































# -------------------------------------------------------
# SEQUENTIAL DESIGN (only one interim analysis allowed)
# -------------------------------------------------------


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
    EN1_rej <- (n1[1] * p_succ_1 + n1[2] * p_succ_2) / p_total
    EN2_rej <- (n2[1] * p_succ_1 + n2[2] * p_succ_2) / p_total

    ## E[N | -]
    ## p_cont_1 - p_succ_2 = P(reach final but fail)
    p_fail_2 <- max(0, p_cont_1 - p_succ_2)

    EN1_norej <- (n1[1] * p_fut_1 + n1[2] * (p_fail_2)) / p_no_rej
    EN2_norej <- (n2[1] * p_fut_1 + n2[2] * (p_fail_2)) / p_no_rej

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









# This works, but only gives the avgT1E, NO sample sizes 


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














# This should also give the sample sizes


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
    avg_Power       <- calc_expectation(function(res) res$Overall$Power)
    avg_Prob_Fut    <- calc_expectation(function(res) res$Overall$Prob_Fut_seq)
    avg_EN_t        <- calc_expectation(function(res) res$Overall$EN_t)
    avg_EN_c        <- calc_expectation(function(res) res$Overall$EN_c)
    avg_EN_t_Succ   <- calc_expectation(function(res) res$Overall$EN_t_Succ)
    avg_EN_c_Succ   <- calc_expectation(function(res) res$Overall$EN_c_Succ)
    avg_EN_t_Fail   <- calc_expectation(function(res) res$Overall$EN_t_Fail)
    avg_EN_c_Fail   <- calc_expectation(function(res) res$Overall$EN_c_Fail)
    
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
    n_stages <- nrow(base_res$Per_Stage)
    
    per_stage_list <- lapply(1:n_stages, function(k) {
      avg_P_Succ <- calc_expectation(function(res) res$Per_Stage$P_Succ[k])
      avg_P_Fut  <- calc_expectation(function(res) res$Per_Stage$P_Fut[k])
      
      data.frame(
        Stage  = k,
        N_Trt  = base_res$Per_Stage$N_Trt[k],  # Fixed cumulative sample sizes
        N_Ctrl = base_res$Per_Stage$N_Ctrl[k],
        P_Succ = avg_P_Succ,
        P_Fut  = avg_P_Fut
      )
    })
    
    per_stage_res <- do.call(rbind, per_stage_list)
    
    # Compute the average cumulative probabilities step-by-step
    per_stage_res$Cum_P_Succ <- cumsum(per_stage_res$P_Succ)
    per_stage_res$Cum_P_Fut  <- cumsum(per_stage_res$P_Fut)
    
    return(list(Overall = overall_res, Per_Stage = per_stage_res))
  }
  
  return(design_eval_fun)
}























