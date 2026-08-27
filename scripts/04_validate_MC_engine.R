# =============================================================================
# Sequential Design — verification of the Monte Carlo engine
# =============================================================================
# This file produces NO thesis results. Its only purpose is to check that the
# Monte Carlo engine reproduces the exact analytical sequential operating
# characteristics, so that 05_sequential_one_interim.R can be trusted.
#
#   reference  : oc2S_seq.dual.normMix  (defined below, exact, conditional)
#   under test : oc2_seq_mc.normMix     (00_functions.R)
#
# Only the CONDITIONAL operating characteristics are compared. That is enough:
# the averaged metrics are obtained by drawing theta_C from the design prior,
# which is exact Monte Carlo integration by construction, so once the conditional
# engine is right the averaged one is too.
#
# Each analytical evaluation costs a few seconds, so the grid is deliberately
# small. Expect a couple of minutes.
#
# Run from the Thesis/ root:
#   Rscript scripts/04_validate_MC_engine.R
# =============================================================================

suppressMessages({
  library(RBesT)
  library(dplyr)
})
source("scripts/00_functions.R")   # oc2_seq_mc.normMix


# -----------------------------------------------------------------------------
# oc2S_seq.dual.normMix — exact/analytical reference, validation-only.
#
# Used ONLY in this file, to check oc2_seq_mc.normMix (00_functions.R)
# against a closed-form calculation. It does not feed any thesis figure or
# table, so it lives here rather than in the shared 00_functions.R, which
# should only contain functions the actual results pipeline depends on.
# -----------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Design parameters (identical to 05_sequential_one_interim.R)
# ---------------------------------------------------------------------------
sigma      <- 88
delta_MCID <- 60

n1_seq <- c(20, 40)   # treatment arm, cumulative
n2_seq <- c(10, 20)   # control arm, cumulative
n_stg1 <- n1_seq[1] + n2_seq[1]
n_stg2 <- n1_seq[2] + n2_seq[2]

n_sim <- 200000       # MC replicates per evaluation


# ---------------------------------------------------------------------------
# Priors
# ---------------------------------------------------------------------------
p_MAP <- mixnorm(
  c(0.4848, -52.457, 21.154), c(0.4598, -47.465, 7.843), c(0.0554, -50.355, 48.164),
  sigma = sigma, param = "ms"
)
p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")


# ---------------------------------------------------------------------------
# Decision criteria (identical to 2.3)
# ---------------------------------------------------------------------------
# Efficacy, dual criterion: Pr(delta > 0) > 0.95 AND Pr(delta > 50) > 0.50
sign.crit <- decision2S(pc = c(0.95, 0.50), qc = c(0, -50), lower.tail = TRUE)
# Futility, stage 1 only: Pr(delta < 40 | data) >= 0.90
fut.crit  <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)

# The analytical function accepts a bare decision2S for the final stage.
decisions_an <- list(list(success = sign.crit, futility = fut.crit), sign.crit)
# The MC engine wants the success/futility slots at every stage.
decisions_mc <- list(list(success = sign.crit, futility = fut.crit),
                     list(success = sign.crit, futility = NULL))


# ---------------------------------------------------------------------------
# Tolerances: everything is compared against 4 Monte Carlo standard errors.
# N takes only the two values n_stg1 and n_stg2, so its variance is Bernoulli.
# ---------------------------------------------------------------------------
se_prop <- function(p) sqrt(p * (1 - p) / n_sim)
se_EN   <- function(p_stop1) sqrt(p_stop1 * (1 - p_stop1) / n_sim) * abs(n_stg2 - n_stg1)

report <- function(label, analytic, mc, tol) {
  d  <- mc - analytic
  ok <- abs(d) < tol
  cat("  ", formatC(label, width = -18),
      " analytic =", formatC(analytic, format = "f", digits = 5, width = 9),
      "  MC =",      formatC(mc,       format = "f", digits = 5, width = 9),
      "  diff =",    formatC(d,        format = "f", digits = 5, width = 9),
      " ", if (ok) "PASS" else "**FAIL**", "\n")
  ok
}


# ---------------------------------------------------------------------------
# Conditional OC: analytical versus MC
# ---------------------------------------------------------------------------
cat("\n===============================================================\n")
cat("Conditional OC, analytical versus MC (n_sim =", n_sim, ")\n")
cat("tolerance: 4 Monte Carlo standard errors\n")
cat("===============================================================\n")

theta_c_check   <- c(-70, -50, -30)
analysis_priors <- list("MAP" = p_MAP, "vague" = p_vague)
all_ok <- TRUE

for (pname in names(analysis_priors)) {

  oc_an <- oc2S_seq.dual.normMix(
    prior1 = p_vague, prior2 = analysis_priors[[pname]],
    n1 = n1_seq, n2 = n2_seq, decisions = decisions_an
  )

  for (tc in theta_c_check) {
    # delta = 0 gives the type I error, delta = delta_MCID gives the power.
    for (d in c(0, delta_MCID)) {

      a <- oc_an(tc - d, tc)                       # theta_T = theta_C - delta
      m <- oc2_seq_mc.normMix(
        theta_1 = tc - d, theta_2 = tc,
        prior_1 = p_vague, prior_2 = analysis_priors[[pname]],
        n1_seq  = n1_seq,  n2_seq  = n2_seq,
        decisions_list = decisions_mc, n_sim = n_sim, seed = 1
      )
      ps <- m$Per_Stage
      p_stop1 <- a$P_Stop_Eff_Stg1 + a$P_Stop_Fut_Stg1

      cat("\n", pname, " | theta_C =", tc, " | delta =", d, "\n")
      all_ok <- report("P_Stop_Eff_Stg1", a$P_Stop_Eff_Stg1, ps$P_Succ[1],
                       4 * se_prop(a$P_Stop_Eff_Stg1)) && all_ok
      all_ok <- report("P_Stop_Fut_Stg1", a$P_Stop_Fut_Stg1, ps$P_Fut[1],
                       4 * se_prop(a$P_Stop_Fut_Stg1)) && all_ok
      all_ok <- report("Total_Power",     a$Total_Power,     m$Overall[["Power"]],
                       4 * se_prop(a$Total_Power)) && all_ok
      all_ok <- report("EN (Trt + Ctrl)", a$EN_Trt + a$EN_Pbo,
                       m$Overall[["EN_t"]] + m$Overall[["EN_c"]],
                       4 * se_EN(p_stop1)) && all_ok
    }
  }
}

cat("\n\n===============================================================\n")
if (all_ok) {
  cat("ALL CHECKS PASSED: the MC engine matches the analytical reference.\n")
} else {
  cat("SOME CHECKS FAILED: inspect the lines marked **FAIL** above.\n")
}
cat("===============================================================\n")
