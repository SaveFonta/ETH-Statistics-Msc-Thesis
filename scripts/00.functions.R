


#' Average Bayesian OC where theta1 = theta2 + delta
#'
#' Computes E_{theta2 ~ mix2}[ OC(theta2 + delta, theta2) ]
#'
#' @param prior1   normMix prior for arm 1
#' @param prior2   normMix prior for arm 2
#' @param n1       sample size arm 1 (treatment)
#' @param n2       sample size arm 2 (placebo)
#' @param decision a decision2S object
#' @param delta    numeric scalar: theta1 = theta2 + delta
#' @param mix2     normMix distribution over which theta2 is averaged
#' @param sigma1   known SD for arm 1 (defaults to prior reference scale)
#' @param sigma2   known SD for arm 2 (defaults to prior reference scale)
#' @param eps      numerical tolerance (tails cut-off probability)
#' @param Ngrid    grid fineness for boundary calculation
#'
#' @return numeric scalar: the averaged OC
avgoc2S.normMix <- function(
  prior1, #analysis prior
  prior2,
  n1,
  n2,
  decision,
  delta,        # <-- theta1 = theta2 + delta
  mix2,         # <-- distribution of theta2
  sigma1,
  sigma2,
  eps   = 1e-6,
  Ngrid = 10,
  ...
) {
  ## 1. SETUP SIGMAS
      ## distributions of the means of the data generating distributions
  ## for now we assume that the underlying standard deviation
  ## matches the respective reference scales

  if (missing(sigma1)) {
    sigma1 <- RBesT::sigma(prior1)
    message("Using default prior1 reference scale ", sigma1)
  }
  if (missing(sigma2)) {
    sigma2 <- RBesT::sigma(prior2)
    message("Using default prior2 reference scale ", sigma2)
  }
  assert_number(sigma1, lower = 0)
  assert_number(sigma2, lower = 0)

  sigma(prior1) <- sigma1
  sigma(prior2) <- sigma2


  ## 2. GET DECISION BOUNDARY
  crit_y1 <- RBesT:::decision2S_boundary(
    prior1, prior2, n1, n2, decision,
    sigma1, sigma2, eps, Ngrid
  )


  ## 3. SPECIFIC SETUP
  sem1 <- sigma1 / sqrt(n1)
  sem2 <- if (n2 == 0) sigma2 / sqrt(1e-1) else sigma2 / sqrt(n2)


  ## 4. DEFINE CORE EVALUATION FUNCTION
  # here I defined freq inside each if, should be a bit faster 

  if (is(decision, "decision2S_1sided")) {
    assert_function(crit_y1)
    lower.tail <- attr(decision, "lower.tail")

    freq <- function(theta1, theta2) {

      ## Set seach area limits for crit_y1 (closure of decision_2S_boundary)
      lim1 <- qnorm(c(eps / 2, 1 - eps / 2), theta1, sem1)
      if (n2 == 0) {
        pnorm(crit_y1(theta2), theta1, sem1, lower.tail = lower.tail)
      } else {
        ## Double arm case: Integrate over normal sampling distribution
        # we have our theta2 which is a Mixture of Normals. But y_2 is a N(theta_2 , SEM_2)
        # so first we evaluate the inner integral 

        RBesT:::integrate_density_log(
          function(x) {
            pnorm(
              crit_y1(x, lim1 = lim1),
              theta1, sem1,
              lower.tail = lower.tail,
              log.p = TRUE
            )
          },
          mixnorm(c(1, theta2, sem2), sigma = sem2), # just a normal N(theta2, sem2)
          logit(eps / 2),
          logit(1 - eps / 2)
        )
      }
    }

  } else {
    # 2-sided boundary
    assert_list(crit_y1, len = 2, types = "function")
    crit_y1_lower_or_equal_than <- crit_y1$lower_or_equal_than
    crit_y1_higher_than          <- crit_y1$higher_than
    
    # in this case the trial to be a success must end up in b_hi < y_1 < b_lo
    freq <- function(theta1, theta2) {
      lim1 <- qnorm(c(eps / 2, 1 - eps / 2), theta1, sem1)
      if (n2 == 0) {
        b_lo <- crit_y1_lower_or_equal_than(theta2)
        b_hi <- crit_y1_higher_than(theta2)
        if (b_lo <= b_hi) return(0)
        pnorm(b_lo, theta1, sem1, lower.tail = TRUE) -
          pnorm(b_hi, theta1, sem1, lower.tail = TRUE)
      } else {
        integrand <- function(x) {
          b_lo <- crit_y1_lower_or_equal_than(x, lim1 = lim1)
          b_hi <- crit_y1_higher_than(x, lim1 = lim1)
          ifelse(
            b_lo <= b_hi,
            -Inf,
            log(
              pnorm(b_lo, theta1, sem1, lower.tail = TRUE) -
                pnorm(b_hi, theta1, sem1, lower.tail = TRUE)
            )
          )
        }
        RBesT:::integrate_density_log(
          integrand,
          mixnorm(c(1, theta2, sem2), sigma = sem2),
          logit(eps / 2),
          logit(1 - eps / 2)
        )
      }
    }
  }


  Vfreq <- Vectorize (freq)

  ## ---- outer integral over theta2 ~ mix2 ------------------------------
  ## theta1 = theta2 + delta, so the integrand is freq(x + delta, x)

  design_fun <- function(delta_new = delta, mix2_new = mix2) { 

  assert_number(delta_new)                       # delta can be any real number

  lim2_range <- qmix(mix2_new, c(eps / 2, 1 - eps / 2)) #contains min and max possible values of y_2
  lim1_range <- lim2_range + delta_new           # theta_2 + delta

  # but y_1 ~ N(theta_1, SEM_1), so we need to add another layer for the cache to where to look for in case of need to extend
  lim1_cache <- c(
    qnorm(eps / 2,     lim1_range[1], sem1),
    qnorm(1 - eps / 2, lim1_range[2], sem1)
  )


  ## Warm up the boundary cache over the range 
  if (is(decision, "decision2S_1sided")) {
    crit_y1(lim2_range, lim1 = lim1_cache)
  } else {
    crit_y1$lower_or_equal_than(lim2_range, lim1 = lim1_cache)
    crit_y1$higher_than(lim2_range,         lim1 = lim1_cache)
  }
  
  ## integral [OC(x + delta, x) p(x | mix2) dx]
  RBesT:::integrate_density_log(
    log_integrand = function(x) log(Vfreq(x + delta_new, x)),
    mix   = mix2_new,
    Lplower = logit(eps / 2),
    Lpupper = logit(1 - eps / 2)
  )
  }
  

  # so basically we return this function where we fix the design (including the anaysis prior) and change delta and the mixture 
  # design_fun is just the metric (5) of the paper from Best ,
  # where we can control the real distribution of theta_c and delta
  return(design_fun)
}




