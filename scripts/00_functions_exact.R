library(RBesT)


# -----------------------------------------------
# EXACT METHODS
# -----------------------------------------------



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

