







#' internal function used for integration of densities which appears
#' to be much more stable from -Inf to +Inf in the logit space while
#' the density to be integrated recieves inputs from 0 to 1 such that
#' the inverse distribution function must be used. The integral solved
#' is int_x dmix(mix,x) integrand(x) where integrand must be given as
#' log and we integrate over the support of mix.
#'
#' integrate density in logit space and split by component such
#' that the quantile function of each component is used. This
#' ensures that the R implementation of the quantile function is
#' always used.
#'
#' @param log_integrand function to integrate over which must return the log(f)
#' @param mix density over which to integrate
#' @param Lplower logit of lower cumulative density
#' @param Lpupper logit of upper cumulative density
#'
#' @keywords internal
integrate_density_log <- function(
  log_integrand,
  mix,
  Lplower = -Inf,
  Lpupper = Inf,
  eps = getOption("RBesT.integrate_prob_eps", 1E-6)
) {
  .integrand_comp_logit <- function(mix_comp) {
    function(l) { # input a logit value
      u <- inv_logit(l) #trasform a logit back to Probability in [0,1]
      
      # We are doing change of variable, so we must multiply by the Jacobian
      lp <- log_inv_logit(l) 
      lnp <- log_inv_logit(-l)   # lp  + lnp  = log (derivative of inverse logit)
      exp(lp + lnp + log_integrand(qmix(mix_comp, u)))
    }
  }

  Nc <- ncol(mix)

  ## integrate by component of mix separatley to increase precision
  ## when the density is not 0 at the boundaries integration, then
  ## the integration is performed on the natural scale. The check
  ## for that is done on the identity scale to avoid numerical
  ## issues.

  lower <- inv_logit(Lplower)
  upper <- inv_logit(Lpupper)
  return(sum(
    vapply(
      1:Nc,
      function(comp) {
        mix_comp <- mix[[comp, rescale = TRUE]]
        fn_integrand_comp_logit <- .integrand_comp_logit(mix_comp)
        if (all(!is.na(fn_integrand_comp_logit(c(Lplower, Lpupper))))) {
          return(.integrate(fn_integrand_comp_logit, Lplower, Lpupper))
        }
        lower_comp <- ifelse(
          Lplower == -Inf,
          qmix(mix_comp, eps),
          qmix(mix_comp, lower)
        )
        upper_comp <- ifelse(
          Lpupper == Inf,
          qmix(mix_comp, 1 - eps),
          qmix(mix_comp, upper)
        )
        return(.integrate(
          function(x) exp(log_integrand(x) + dmix(mix_comp, x, log = TRUE)),
          lower_comp,
          upper_comp
        ))
      },
      c(0.1)
    ) *
      mix[1, ]
  ))
}

integrate_density <- function(
  integrand,
  mix,
  Lplower = -Inf,
  Lpupper = Inf,
  eps = getOption("RBesT.integrate_prob_eps", 1E-6)
) {
  .integrand_comp_logit <- function(mix_comp) {
    function(l) {
      u <- inv_logit(l)
      lp <- log_inv_logit(l)
      lnp <- log_inv_logit(-l)
      exp(lp + lnp) * integrand(qmix(mix_comp, u))
    }
  }
  Nc <- ncol(mix)

  lower <- inv_logit(Lplower)
  upper <- inv_logit(Lpupper)

  return(sum(
    vapply(
      1:Nc,
      function(comp) {
        mix_comp <- mix[[comp, rescale = TRUE]]
        ## ensure that the integrand is defined at the boundaries...
        fn_integrand_comp_logit <- .integrand_comp_logit(mix_comp)
        if (all(!is.na(fn_integrand_comp_logit(c(Lplower, Lpupper))))) {
          return(.integrate(fn_integrand_comp_logit, Lplower, Lpupper))
        }
        ## ... otherwise we avoid the boundaries by eps prob density:
        lower_comp <- ifelse(
          Lplower == -Inf,
          qmix(mix_comp, eps),
          qmix(mix_comp, lower)
        )
        upper_comp <- ifelse(
          Lpupper == Inf,
          qmix(mix_comp, 1 - eps),
          qmix(mix_comp, upper)
        )
        return(.integrate(
          function(x) integrand(x) * dmix(mix_comp, x),
          lower_comp,
          upper_comp
        ))
      },
      c(0.1)
    ) *
      mix[1, ]
  ))
}

.integrate <- function(integrand, lower, upper) {
  integrate_args_user <- getOption("RBesT.integrate_args", list())
  args <- modifyList(
    list(
      lower = lower,
      upper = upper,
      rel.tol = .Machine$double.eps^0.25,
      abs.tol = .Machine$double.eps^0.25,
      subdivisions = 1000,
      stop.on.error = TRUE
    ),
    integrate_args_user
  )

  integrate(
    integrand,
    lower = args$lower,
    upper = args$upper,
    rel.tol = args$rel.tol,
    abs.tol = args$abs.tol,
    subdivisions = args$subdivisions,
    stop.on.error = args$stop.on.error
  )$value
}





# -----------------------------------
# Now I just unify the two integrating function, stupid not to

integrate_density <- function(
  integrand,
  mix,
  Lplower = -Inf,
  Lpupper = Inf,
  log = FALSE,
  eps = getOption("RBesT.integrate_prob_eps", 1E-6)
) {
  .integrand_comp_logit <- function(mix_comp) {
    function(l) {
      u <- inv_logit(l)
      lp <- log_inv_logit(l)
      lnp <- log_inv_logit(-l)
      if (log) {
        exp(lp + lnp + integrand(qmix(mix_comp, u)))
      } else {
        exp(lp + lnp) * integrand(qmix(mix_comp, u))
      }
    }
  }

  Nc <- ncol(mix)
  lower <- inv_logit(Lplower)
  upper <- inv_logit(Lpupper)

  return(sum(
    vapply(
      1:Nc,
      function(comp) {
        mix_comp <- mix[[comp, rescale = TRUE]]
        fn_integrand_comp_logit <- .integrand_comp_logit(mix_comp)
        if (all(!is.na(fn_integrand_comp_logit(c(Lplower, Lpupper))))) {
          return(.integrate(fn_integrand_comp_logit, Lplower, Lpupper))
        }
        lower_comp <- ifelse(
          Lplower == -Inf,
          qmix(mix_comp, eps),
          qmix(mix_comp, lower)
        )
        upper_comp <- ifelse(
          Lpupper == Inf,
          qmix(mix_comp, 1 - eps),
          qmix(mix_comp, upper)
        )
        if (log) {
          return(.integrate(
            function(x) exp(integrand(x) + dmix(mix_comp, x, log = TRUE)),
            lower_comp,
            upper_comp
          ))
        } else {
          return(.integrate(
            function(x) integrand(x) * dmix(mix_comp, x),
            lower_comp,
            upper_comp
          ))
        }
      },
      c(0.1)
    ) *
      mix[1, ]
  ))
}













## returns a function object which is the decision boundary. That is
## the function finds at a regular grid between llim1 and ulim1 the
## roots of the decision function and returns an interpolation
## function object

# this is the actual function that created the grid and searches with uniroot
solve_boundary2S_normMix <- function(
  decision,
  mix1,
  mix2,
  n1,
  n2,
  lim1,
  lim2,
  delta2
) {
  assert_class(decision, "decision2S_atomic")

  grid <- seq(lim2[1], lim2[2], length = diff(lim2) / delta2)

  sigma1 <- sigma(mix1)
  sigma2 <- sigma(mix2)

  sem1 <- sigma1 / sqrt(n1)
  scale1 <- sigma1 / (n1^0.25)

  cond_decisionStep <- function(post2) {
    fn <- function(m1) {
      decision(postmix(mix1, m = m1, se = sem1), post2) - 0.75
    }
    Vectorize(fn)
  }

  Neval <- length(grid)
  # cat("Calculating boundary from", lim2[1], "to", lim2[2], "with", Neval, "points\n")
  tol <- min(delta2 / 100, .Machine$double.eps^0.25)
  ## cat("Using tolerance", tol, "\n")
  crit <- rep(NA, times = Neval)
  for (i in 1:Neval) {
    if (n2 == 0) {
      post2 <- mix2
    } else {
      post2 <- postmix(mix2, m = grid[i], se = sigma2 / sqrt(n2))
    }
    ind_fun <- cond_decisionStep(post2)
    dec_bounds <- ind_fun(lim1)
    ## if decision function is not different at boundaries, lim1
    ## is too narrow and we then enlarge
    while (prod(dec_bounds) > 0) {
      w <- diff(lim1)
      lim1 <- c(lim1[1] - w / 2, lim1[2] + w / 2)
      dec_bounds <- ind_fun(lim1)
    }
    y1c <- uniroot(
      ind_fun,
      lim1,
      f.lower = dec_bounds[1],
      f.upper = dec_bounds[2],
      tol = tol
    )$root
    crit[i] <- y1c
    ## set lim1 tightly around the current critical value and use the
    ## last boundary limits to not shrink too fast
    lim1 <- c(mean(lim1[1], y1c - 2 * scale1), mean(y1c + 2 * scale1, lim1[2]))
  }

  cbind(grid, crit)
}

# this just looks at the 'decision' object and send it to 1sided or 2sided functions 
decision2S_boundary.normMix <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps = 1e-6,
  Ngrid = 10,
  ...
) {
  ## distributions of the means of the data generating distributions
  ## for now we assume that the underlying standard deviation
  ## matches the respective reference scales
  if (missing(sigma1)) {
    sigma1 <- RBesT::sigma(prior1)
    message("Using default prior 1 reference scale ", sigma1)
  }
  if (missing(sigma2)) {
    sigma2 <- RBesT::sigma(prior2)
    message("Using default prior 2 reference scale ", sigma2)
  }

  if (is(decision, "decision2S_2sided")) {
    decision2S_boundary_normMix_2sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      sigma1,
      sigma2,
      eps,
      Ngrid,
      ...
    )
  } else {
    decision2S_boundary_normMix_1sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      sigma1,
      sigma2,
      eps,
      Ngrid,
      ...
    )
  }
}

decision2S_boundary_normMix_2sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps,
  Ngrid,
  ...
) {
  crit_lower <- decision2S_boundary_normMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    lower(decision),
    sigma1,
    sigma2,
    eps,
    Ngrid,
    ...
  )
  crit_upper <- decision2S_boundary_normMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    upper(decision),
    sigma1,
    sigma2,
    eps,
    Ngrid,
    ...
  )
  list(lower_or_equal_than = crit_lower, higher_than = crit_upper)
}




# function tat actually performs the whole calcualtion and cache using <<-
decision2S_boundary_normMix_atomic <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps,
  Ngrid,
  ...
) {
  assert_class(decision, "decision2S_atomic")

  assert_number(sigma1, lower = 0)
  assert_number(sigma2, lower = 0)

  sem1 <- sigma1 / sqrt(n1)
  sem2 <- sigma2 / sqrt(n2)

  sigma(prior1) <- sigma1
  sigma(prior2) <- sigma2

  ## only n2 can be zero
  assert_that(n1 > 0)
  assert_that(n2 >= 0)

  if (n2 == 0) {
    sem2 <- sigma(prior2) / sqrt(1E-1)
  }

  ## change the reference scale of the prior such that the prior
  ## represents the distribution of the respective means
  mean_prior1 <- prior1
  sigma(mean_prior1) <- sem1
  ## mean_prior2 <- prior2
  ## sigma(mean_prior2) <- sem2

  ## discretization step-size
  delta2 <- sem2 / Ngrid

  ## for the case of mix1 and mix2 having just 1 component, then one
  ## can prove that the decision boundary is a linear function.
  ## Hence we only calculate a very rough grid and apply linear
  ## interpolation.

  linear_boundary <- FALSE
  if (ncol(prior1) == 1 && ncol(prior2) == 1) {
    linear_boundary <- TRUE
    ## we could relax this even further
    delta2 <- sigma2 / Ngrid
  }

  ## the boundary function depends only on the samples sizes n1, n2,
  ## the priors and the decision, but not the assumed truths

  clim2 <- c(Inf, -Inf)

  ## the boundary function which gives conditional on the second
  ## variable the critical value where the decision changes
  boundary <- NA
  boundary_discrete <- matrix(NA, nrow = 0, ncol = 2)

  decision_boundary <- function(y2, lim1) {
    lim2 <- range(y2)

    ## check if boundary function must be recomputed
    if (lim2[1] < clim2[1] | lim2[2] > clim2[2]) {
      new_lim2 <- clim2
      ## note: the <<- assignment is needed to set the variable in the enclosure
      if (missing(lim1)) {
        lim1 <- qmix(mean_prior1, c(eps / 2, 1 - eps / 2))
      }
      if (nrow(boundary_discrete) == 0) {
        ## boundary hasn't been calculated before, do it all
        boundary_discrete <<- solve_boundary2S_normMix(
          decision,
          prior1,
          prior2,
          n1,
          n2,
          lim1,
          lim2,
          delta2
        )
        new_lim2 <- lim2
      } else {
        if (lim2[1] < clim2[1]) {
          ## the lower bound is not low enough... only add the region which is missing
          new_left_lim2 <- min(lim2[1], clim2[1] - 2 * delta2)
          boundary_extra <- solve_boundary2S_normMix(
            decision,
            prior1,
            prior2,
            n1,
            n2,
            lim1,
            c(new_left_lim2, clim2[1] - delta2),
            delta2
          )
          new_lim2[1] <- new_left_lim2
          boundary_discrete <<- rbind(boundary_extra, boundary_discrete)
        }
        if (lim2[2] > clim2[2]) {
          ## the upper bound is not large enough.. again only add what's missing
          new_right_lim2 <- max(lim2[2], clim2[2] + 2 * delta2)
          boundary_extra <- solve_boundary2S_normMix(
            decision,
            prior1,
            prior2,
            n1,
            n2,
            lim1,
            c(clim2[2] + delta2, new_right_lim2),
            delta2
          )
          new_lim2[2] <- new_right_lim2
          boundary_discrete <<- rbind(boundary_discrete, boundary_extra)
        }
      }
      ## only for debugging
      ## assert_that(all(order(boundary_discrete[,1]) == 1:nrow(boundary_discrete)), msg="x grid must stay ordered!")
      if (linear_boundary) {
        boundary <<- approxfun(
          boundary_discrete[, 1],
          boundary_discrete[, 2],
          rule = 2
        )
      } else {
        boundary <<- splinefun(boundary_discrete[, 1], boundary_discrete[, 2])
      }
      clim2 <<- new_lim2
    }

    return(boundary(y2))
  }

  decision_boundary
}






































































#' internal function used for integration of densities which appears
#' to be much more stable from -Inf to +Inf in the logit space while
#' the density to be integrated recieves inputs from 0 to 1 such that
#' the inverse distribution function must be used. The integral solved
#' is int_x dmix(mix,x) integrand(x) where integrand must be given as
#' log and we integrate over the support of mix.
#'
#' integrate density in logit space and split by component such
#' that the quantile function of each component is used. This
#' ensures that the R implementation of the quantile function is
#' always used.
#'
#' @param log_integrand function to integrate over which must return the log(f)
#' @param mix density over which to integrate
#' @param Lplower logit of lower cumulative density
#' @param Lpupper logit of upper cumulative density
#'
#' @keywords internal
integrate_density_log <- function(
  log_integrand,
  mix,
  Lplower = -Inf,
  Lpupper = Inf,
  eps = getOption("RBesT.integrate_prob_eps", 1E-6)
) {
  .integrand_comp_logit <- function(mix_comp) {
    function(l) { # input a logit value
      u <- inv_logit(l) #trasform a logit back to Probability in [0,1]
      
      # We are doing change of variable, so we must multiply by the Jacobian
      lp <- log_inv_logit(l) 
      lnp <- log_inv_logit(-l)   # lp  + lnp  = log (derivative of inverse logit)
      exp(lp + lnp + log_integrand(qmix(mix_comp, u)))
    }
  }

  Nc <- ncol(mix)

  ## integrate by component of mix separatley to increase precision
  ## when the density is not 0 at the boundaries integration, then
  ## the integration is performed on the natural scale. The check
  ## for that is done on the identity scale to avoid numerical
  ## issues.

  lower <- inv_logit(Lplower)
  upper <- inv_logit(Lpupper)
  return(sum(
    vapply(
      1:Nc,
      function(comp) {
        mix_comp <- mix[[comp, rescale = TRUE]]
        fn_integrand_comp_logit <- .integrand_comp_logit(mix_comp)
        if (all(!is.na(fn_integrand_comp_logit(c(Lplower, Lpupper))))) {
          return(.integrate(fn_integrand_comp_logit, Lplower, Lpupper))
        }
        lower_comp <- ifelse(
          Lplower == -Inf,
          qmix(mix_comp, eps),
          qmix(mix_comp, lower)
        )
        upper_comp <- ifelse(
          Lpupper == Inf,
          qmix(mix_comp, 1 - eps),
          qmix(mix_comp, upper)
        )
        return(.integrate(
          function(x) exp(log_integrand(x) + dmix(mix_comp, x, log = TRUE)),
          lower_comp,
          upper_comp
        ))
      },
      c(0.1)
    ) *
      mix[1, ]
  ))
}

integrate_density <- function(
  integrand,
  mix,
  Lplower = -Inf,
  Lpupper = Inf,
  eps = getOption("RBesT.integrate_prob_eps", 1E-6)
) {
  .integrand_comp_logit <- function(mix_comp) {
    function(l) {
      u <- inv_logit(l)
      lp <- log_inv_logit(l)
      lnp <- log_inv_logit(-l)
      exp(lp + lnp) * integrand(qmix(mix_comp, u))
    }
  }
  Nc <- ncol(mix)

  lower <- inv_logit(Lplower)
  upper <- inv_logit(Lpupper)

  return(sum(
    vapply(
      1:Nc,
      function(comp) {
        mix_comp <- mix[[comp, rescale = TRUE]]
        ## ensure that the integrand is defined at the boundaries...
        fn_integrand_comp_logit <- .integrand_comp_logit(mix_comp)
        if (all(!is.na(fn_integrand_comp_logit(c(Lplower, Lpupper))))) {
          return(.integrate(fn_integrand_comp_logit, Lplower, Lpupper))
        }
        ## ... otherwise we avoid the boundaries by eps prob density:
        lower_comp <- ifelse(
          Lplower == -Inf,
          qmix(mix_comp, eps),
          qmix(mix_comp, lower)
        )
        upper_comp <- ifelse(
          Lpupper == Inf,
          qmix(mix_comp, 1 - eps),
          qmix(mix_comp, upper)
        )
        return(.integrate(
          function(x) integrand(x) * dmix(mix_comp, x),
          lower_comp,
          upper_comp
        ))
      },
      c(0.1)
    ) *
      mix[1, ]
  ))
}

.integrate <- function(integrand, lower, upper) {
  integrate_args_user <- getOption("RBesT.integrate_args", list())
  args <- modifyList(
    list(
      lower = lower,
      upper = upper,
      rel.tol = .Machine$double.eps^0.25,
      abs.tol = .Machine$double.eps^0.25,
      subdivisions = 1000,
      stop.on.error = TRUE
    ),
    integrate_args_user
  )

  integrate(
    integrand,
    lower = args$lower,
    upper = args$upper,
    rel.tol = args$rel.tol,
    abs.tol = args$abs.tol,
    subdivisions = args$subdivisions,
    stop.on.error = args$stop.on.error
  )$value
}
















## returns a function object which is the decision boundary. That is
## the function finds at a regular grid between llim1 and ulim1 the
## roots of the decision function and returns an interpolation
## function object

# this is the actual function that created the grid and searches with uniroot
solve_boundary2S_normMix <- function(
  decision,
  mix1,
  mix2,
  n1,
  n2,
  lim1,
  lim2,
  delta2
) {
  assert_class(decision, "decision2S_atomic")

  grid <- seq(lim2[1], lim2[2], length = diff(lim2) / delta2)

  sigma1 <- sigma(mix1)
  sigma2 <- sigma(mix2)

  sem1 <- sigma1 / sqrt(n1)
  scale1 <- sigma1 / (n1^0.25)

  cond_decisionStep <- function(post2) {
    fn <- function(m1) {
      decision(postmix(mix1, m = m1, se = sem1), post2) - 0.75
    }
    Vectorize(fn)
  }

  Neval <- length(grid)
  # cat("Calculating boundary from", lim2[1], "to", lim2[2], "with", Neval, "points\n")
  tol <- min(delta2 / 100, .Machine$double.eps^0.25)
  ## cat("Using tolerance", tol, "\n")
  crit <- rep(NA, times = Neval)
  for (i in 1:Neval) {
    if (n2 == 0) {
      post2 <- mix2
    } else {
      post2 <- postmix(mix2, m = grid[i], se = sigma2 / sqrt(n2))
    }
    ind_fun <- cond_decisionStep(post2)
    dec_bounds <- ind_fun(lim1)
    ## if decision function is not different at boundaries, lim1
    ## is too narrow and we then enlarge
    while (prod(dec_bounds) > 0) {
      w <- diff(lim1)
      lim1 <- c(lim1[1] - w / 2, lim1[2] + w / 2)
      dec_bounds <- ind_fun(lim1)
    }
    y1c <- uniroot(
      ind_fun,
      lim1,
      f.lower = dec_bounds[1],
      f.upper = dec_bounds[2],
      tol = tol
    )$root
    crit[i] <- y1c
    ## set lim1 tightly around the current critical value and use the
    ## last boundary limits to not shrink too fast
    lim1 <- c(mean(lim1[1], y1c - 2 * scale1), mean(y1c + 2 * scale1, lim1[2]))
  }

  cbind(grid, crit)
}

# this just looks at the 'decision' object and send it to 1sided or 2sided functions 
decision2S_boundary.normMix <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps = 1e-6,
  Ngrid = 10,
  ...
) {
  ## distributions of the means of the data generating distributions
  ## for now we assume that the underlying standard deviation
  ## matches the respective reference scales
  if (missing(sigma1)) {
    sigma1 <- RBesT::sigma(prior1)
    message("Using default prior 1 reference scale ", sigma1)
  }
  if (missing(sigma2)) {
    sigma2 <- RBesT::sigma(prior2)
    message("Using default prior 2 reference scale ", sigma2)
  }

  if (is(decision, "decision2S_2sided")) {
    decision2S_boundary_normMix_2sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      sigma1,
      sigma2,
      eps,
      Ngrid,
      ...
    )
  } else {
    decision2S_boundary_normMix_1sided(
      prior1,
      prior2,
      n1,
      n2,
      decision,
      sigma1,
      sigma2,
      eps,
      Ngrid,
      ...
    )
  }
}



decision2S_boundary_normMix_1sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps,
  Ngrid,
  ...
) {
  decision <- if (has_lower(decision)) {
    lower(decision)
  } else {
    upper(decision)
  }
  decision2S_boundary_normMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    decision,
    sigma1,
    sigma2,
    eps,
    Ngrid,
    ...
  )
}


decision2S_boundary_normMix_2sided <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps,
  Ngrid,
  ...
) {
  crit_lower <- decision2S_boundary_normMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    lower(decision),
    sigma1,
    sigma2,
    eps,
    Ngrid,
    ...
  )
  crit_upper <- decision2S_boundary_normMix_atomic(
    prior1,
    prior2,
    n1,
    n2,
    upper(decision),
    sigma1,
    sigma2,
    eps,
    Ngrid,
    ...
  )
  list(lower_or_equal_than = crit_lower, higher_than = crit_upper)
}


decision2S_boundary_normMix_atomic <- function(
  prior1,
  prior2,
  n1,
  n2,
  decision,
  sigma1,
  sigma2,
  eps,
  Ngrid,
  ...
) {
  assert_class(decision, "decision2S_atomic")

  assert_number(sigma1, lower = 0)
  assert_number(sigma2, lower = 0)

  sem1 <- sigma1 / sqrt(n1)
  sem2 <- sigma2 / sqrt(n2)

  sigma(prior1) <- sigma1
  sigma(prior2) <- sigma2

  ## only n2 can be zero
  assert_that(n1 > 0)
  assert_that(n2 >= 0)

  if (n2 == 0) {
    sem2 <- sigma(prior2) / sqrt(1E-1)
  }

  ## change the reference scale of the prior such that the prior
  ## represents the distribution of the respective means
  mean_prior1 <- prior1
  sigma(mean_prior1) <- sem1
  ## mean_prior2 <- prior2
  ## sigma(mean_prior2) <- sem2

  ## discretization step-size
  delta2 <- sem2 / Ngrid

  ## for the case of mix1 and mix2 having just 1 component, then one
  ## can prove that the decision boundary is a linear function.
  ## Hence we only calculate a very rough grid and apply linear
  ## interpolation.

  linear_boundary <- FALSE
  if (ncol(prior1) == 1 && ncol(prior2) == 1) {
    linear_boundary <- TRUE
    ## we could relax this even further
    delta2 <- sigma2 / Ngrid
  }

  ## the boundary function depends only on the samples sizes n1, n2,
  ## the priors and the decision, but not the assumed truths

  clim2 <- c(Inf, -Inf)

  ## the boundary function which gives conditional on the second
  ## variable the critical value where the decision changes
  boundary <- NA
  boundary_discrete <- matrix(NA, nrow = 0, ncol = 2)

  decision_boundary <- function(y2, lim1) {
    lim2 <- range(y2)

    ## check if boundary function must be recomputed
    if (lim2[1] < clim2[1] | lim2[2] > clim2[2]) {
      new_lim2 <- clim2
      ## note: the <<- assignment is needed to set the variable in the enclosure
      if (missing(lim1)) {
        lim1 <- qmix(mean_prior1, c(eps / 2, 1 - eps / 2))
      }
      if (nrow(boundary_discrete) == 0) {
        ## boundary hasn't been calculated before, do it all
        boundary_discrete <<- solve_boundary2S_normMix(
          decision,
          prior1,
          prior2,
          n1,
          n2,
          lim1,
          lim2,
          delta2
        )
        new_lim2 <- lim2
      } else {
        if (lim2[1] < clim2[1]) {
          ## the lower bound is not low enough... only add the region which is missing
          new_left_lim2 <- min(lim2[1], clim2[1] - 2 * delta2)
          boundary_extra <- solve_boundary2S_normMix(
            decision,
            prior1,
            prior2,
            n1,
            n2,
            lim1,
            c(new_left_lim2, clim2[1] - delta2),
            delta2
          )
          new_lim2[1] <- new_left_lim2
          boundary_discrete <<- rbind(boundary_extra, boundary_discrete)
        }
        if (lim2[2] > clim2[2]) {
          ## the upper bound is not large enough.. again only add what's missing
          new_right_lim2 <- max(lim2[2], clim2[2] + 2 * delta2)
          boundary_extra <- solve_boundary2S_normMix(
            decision,
            prior1,
            prior2,
            n1,
            n2,
            lim1,
            c(clim2[2] + delta2, new_right_lim2),
            delta2
          )
          new_lim2[2] <- new_right_lim2
          boundary_discrete <<- rbind(boundary_discrete, boundary_extra)
        }
      }
      ## only for debugging
      ## assert_that(all(order(boundary_discrete[,1]) == 1:nrow(boundary_discrete)), msg="x grid must stay ordered!")
      if (linear_boundary) {
        boundary <<- approxfun(
          boundary_discrete[, 1],
          boundary_discrete[, 2],
          rule = 2
        )
      } else {
        boundary <<- splinefun(boundary_discrete[, 1], boundary_discrete[, 2])
      }
      clim2 <<- new_lim2
    }

    return(boundary(y2))
  }

  decision_boundary
}










###############################
# ----------------------------------------------
# INSPECTION of oc2s for normMix case 
# ------------------------------------------------

# NOTE THAT sigma is the sd
# it is the fixed parameter for which we assume the measuremnt noise of the endpoint.
# we found it is 88, so we assume also the DGP has 88. 
# so we are assuming that the true population we are testing has the same variability as the one I told the prior to expect

oc2S.normMix <- function(
  prior1, prior2, n1, n2, decision,
  sigma1, sigma2, eps = 1e-6, Ngrid = 10, ...
) {
  ## 1. SETUP SIGMAS
    ## distributions of the means of the data generating distributions
  ## for now we assume that the underlying standard deviation
  ## matches the respective reference scales


  if (missing(sigma1)) {
    sigma1 <- RBesT::sigma(prior1)
    message("Using default prior 1 reference scale ", sigma1)
  }
  if (missing(sigma2)) {
    sigma2 <- RBesT::sigma(prior2)
    message("Using default prior 2 reference scale ", sigma2)
  }
  assert_number(sigma1, lower = 0)
  assert_number(sigma2, lower = 0)
  
  sigma(prior1) <- sigma1
  sigma(prior2) <- sigma2
  
  ## 2. GET DECISION BOUNDARY
  crit_y1 <- decision2S_boundary(
    prior1, prior2, n1, n2, decision, sigma1, sigma2, eps, Ngrid
  )

  ## 3. SPECIFIC SETUP
  ## oc2S needs explicit standard errors for the sampling distribution
  sem1 <- sigma1 / sqrt(n1)
  if (n2 == 0) {
    sem2 <- sigma(prior2) / sqrt(1E-1)
  } else {
    sem2 <- sigma2 / sqrt(n2)
  }
  
  ## 4. DEFINE CORE EVALUATION FUNCTION
  freq <- if (is(decision, "decision2S_1sided")) {
    assert_function(crit_y1)
    lower.tail <- attr(decision, "lower.tail")
    
    ## oc2S evaluates over fixed parameters (theta1, theta2)
    function(theta1, theta2) {
      
      ## [Step skipped in oc2S: preddist is not used]
      ## since here we just use normal distribution to weight the integral 
      ##
      ##
      ##
      ##
      ##
      ##
      ##
      
      ## Set seach area limits for crit_y1 (closure of decision_2S_boundary)
      lim1 <- qnorm(c(eps / 2, 1 - eps / 2), theta1, sem1)
      ## [Step skipped in oc2S: lim2 is generated in wrapper later]
      ##
      
      if (n2 == 0) {
        ## Single arm case
        ## [Step skipped in oc2S: mean_prior2 is passed as theta2]
        pnorm(
          crit_y1(theta2),
          theta1, sem1,
          lower.tail = lower.tail
        )
      } else {
        ## Double arm case: Integrate over normal sampling distribution
        # this is our integral to evaluate. Basically x = placebo, then you do CDF evaluated at crit_y1(x, lim1 = lim1) = Decision value 

        integrate_density_log(
          function(x) {
            pnorm(    ## CDF of a Normal 
              crit_y1(x, lim1 = lim1),
              theta1, sem1,
              lower.tail = lower.tail,
              log.p = TRUE
            )
          },
          mixnorm(c(1, theta2, sem2), sigma = sem2), # just a normal N(theta2, sem2)
          Lplower = logit(eps / 2),
          Lpupper = logit(1 - eps / 2)
        )
      }
    }
    
  } else {
    assert_list(crit_y1, len = 2, types = "function")
    crit_y1_lower_or_equal_than <- crit_y1$lower_or_equal_than
    crit_y1_higher_than         <- crit_y1$higher_than
    
    ## oc2S evaluates over fixed parameters (theta1, theta2)
    function(theta1, theta2) {
      assert_scalar(theta1)
      assert_scalar(theta2)
      
      ## [Step skipped in oc2S: preddist is not used]
      ##
      ##
      ##
      ##
      ##
      ##
      ##
      
      ## Integration limits 
      lim1 <- qnorm(c(eps / 2, 1 - eps / 2), theta1, sem1)
      ## [Step skipped in oc2S: lim2 is generated in wrapper later]
      ##
      ##

      if (n2 == 0) {
        ## Single arm case
        ## [Step skipped in oc2S: mean_prior2 is passed as theta2]
        bound_lower_or_equal_than <- crit_y1_lower_or_equal_than(theta2)
        bound_higher_than         <- crit_y1_higher_than(theta2)
        
        if (bound_lower_or_equal_than <= bound_higher_than) {
          0
        } else {
          pnorm(bound_lower_or_equal_than, theta1, sem1, lower.tail = TRUE) -
            pnorm(bound_higher_than, theta1, sem1, lower.tail = TRUE)
        }
      } else {
        ## Double arm case: Integrate over normal sampling distribution
        integrand <- function(x) {
          bound_lower_or_equal_than <- crit_y1_lower_or_equal_than(x, lim1 = lim1)
          bound_higher_than         <- crit_y1_higher_than(x, lim1 = lim1)
          
          ifelse(
            bound_lower_or_equal_than <= bound_higher_than,
            -Inf,
            log(
              pnorm(bound_lower_or_equal_than, theta1, sem1, lower.tail = TRUE) -
                pnorm(bound_higher_than, theta1, sem1, lower.tail = TRUE)
            )
          )
        }
        integrate_density_log(
          integrand,
          mixnorm(c(1, theta2, sem2), sigma = sem2),
          logit(eps / 2),
          logit(1 - eps / 2)
        )
      }
    }
  }

  ## 5. WRAPPING AND VECTORIZATION
  Vfreq <- Vectorize(freq)
  

  # actual function that will be returned 

  design_fun <- function(theta1, theta2, y2) {
    if (!missing(y2)) {
      theta2 <- y2
    }
    if (n2 == 0 && missing(theta2)) {
      theta2 <- theta1
    }


    ## What we want to do now is to use at the entire rage of theta2 requested, then calls crit_y1 for the entire range at once
    # let's say theta1 seq(10, 50), then here it finds the absolute lowest possible data point for 10 and the absolute highest for 50 
    # since we need to integrate a Normal centered at those theta values, lim are the limits of the normals centered
    # in theta1 and theta2 that range form the (eps/2)-th percentile to the (1-eps/2)-th percentile

    
    lim1 <- c(qnorm(eps / 2, min(theta1), sem1), qnorm(1 - eps / 2, max(theta1), sem1))

        # same for theta2
    lim2 <- c(qnorm(eps / 2, min(theta2), sem2), qnorm(1 - eps / 2, max(theta2), sem2))
    

    # This is how crit_y1(x, lim1 = lim1) works internally: 
    # First it looks at x (placebo). It checks in the Cache if it already calculated the y_1,c for this Placebo

    # If it knows the answer: it outputs it ignoring lim1
    # If it doesnt: it uses lim1 as the starting boundaries for new calculations and extend the cache 

    if (is(decision, "decision2S_1sided")) {
      crit_y1(lim2, lim1 = lim1)
    } else {
      crit_y1$lower_or_equal_than(lim2, lim1 = lim1)
      crit_y1$higher_than(lim2, lim1 = lim1)
    }
    
    theta_df <- try(data.frame(theta1 = theta1, theta2 = theta2, row.names = NULL))
    if (inherits(theta_df, "try-error")) {
      stop("theta1 and theta2 need to be of same size")
    }
    do.call(Vfreq, theta_df)  # Vfreq(theta1 = theta_df$theta1, theta2 = theta_df$theta2)
  }

  design_fun
}





provo <- function (x) if (x < 0)  return("ok")
provo(c(10,20))
provo <- Vectorize(provo)
provo(c(10,20))















































log_inv_logit <- function(mat) {
  ## - ifelse(is.finite(mat) & (mat < 0), log1p(exp(mat)) - mat, log1p(exp(-mat)))
  ## idx <- is.finite(mat) & (mat < 0)
  idx <- mat < 0
  mat[idx] <- mat[idx] - log1p(exp(mat[idx]))
  mat[!idx] <- -1 * log1p(exp(-mat[!idx]))
  mat
}





































































pos2S.normMix <- function(
  prior1, prior2, n1, n2, decision,
  sigma1, sigma2, eps = 1e-6, Ngrid = 10, ...
) {
  ## 1. SETUP SIGMAS

    ## distributions of the means of the data generating distributions
  ## for now we assume that the underlying standard deviation
  ## matches the respective reference scales


  if (missing(sigma1)) {
    sigma1 <- RBesT::sigma(prior1)
    message("Using default prior 1 reference scale ", sigma1)
  }
  if (missing(sigma2)) {
    sigma2 <- RBesT::sigma(prior2)
    message("Using default prior 2 reference scale ", sigma2)
  }
  assert_number(sigma1, lower = 0)
  assert_number(sigma2, lower = 0)
  
  sigma(prior1) <- sigma1
  sigma(prior2) <- sigma2
  
  ## 2. GET DECISION BOUNDARY
  crit_y1 <- decision2S_boundary(
    prior1, prior2, n1, n2, decision, sigma1, sigma2, eps, Ngrid
  )

  ## 3. SPECIFIC SETUP
  ## [Step skipped in pos2S: preddist handles sem internally]
  ##
  ##
  ##
  ##
  ##
  ##
  
  ## 4. DEFINE CORE EVALUATION FUNCTION
  design_fun <- if (is(decision, "decision2S_1sided")) {
    assert_function(crit_y1)
    lower.tail <- attr(decision, "lower.tail")
    
    ## pos2S evaluates over mixture distributions (mix1, mix2) not true theta values
    function(mix1, mix2) {
      
      # NOTE:
      ## distributions of the means of the data generating distributions
      ## for now we assume that the underlying standard deviation
      ## matches the respective reference scales

      ## Get predictive distributions
      pred_mix1_mean <- preddist(mix1, n = n1, sigma = sigma1)
      if (n2 == 0) {
        pred_mix2_mean <- preddist(mix2, n = 1,  sigma = sigma2)
      } else {
        pred_mix2_mean <- preddist(mix2, n = n2, sigma = sigma2)
      }
      assert_that(inherits(pred_mix1_mean, "normMix"))
      assert_that(inherits(pred_mix2_mean, "normMix"))
      
      ## Integration limits 
      lim1 <- qmix(pred_mix1_mean, c(eps / 2, 1 - eps / 2))
      lim2 <- qmix(pred_mix2_mean, c(eps / 2, 1 - eps / 2))

      # cache the intehgration limits, so we are fast when we actually integrate in the next step
      crit_y1(x = lim2, lim1 = lim1)

      if (n2 == 0) {
        ## Single arm case
        mean_prior2 <- summary(prior2, probs = c())["mean"]
        pmix(
          pred_mix1_mean,
          crit_y1(mean_prior2),
          lower.tail = lower.tail
        )
      } else {
        ## Double arm case: Integrate over predictive distribution
        integrate_density_log(
          function(x) {
            pmix(
              pred_mix1_mean,
              crit_y1(x, lim1 = lim1),
              lower.tail = lower.tail,
              log.p = TRUE
            )
          },
          pred_mix2_mean,
          logit(eps / 2),
          logit(1 - eps / 2)
        )
      }
    }
    
  } else {
    assert_list(crit_y1, len = 2, types = "function")
    crit_y1_lower_or_equal_than <- crit_y1$lower_or_equal_than
    crit_y1_higher_than         <- crit_y1$higher_than
    
    ## pos2S evaluates over mixture distributions (mix1, mix2)
    function(mix1, mix2) {
      
      ## Get predictive distributions
      pred_mix1_mean <- preddist(mix1, n = n1, sigma = sigma1)
      if (n2 == 0) {
        pred_mix2_mean <- preddist(mix2, n = 1,  sigma = sigma2)
      } else {
        pred_mix2_mean <- preddist(mix2, n = n2, sigma = sigma2)
      }
      assert_that(inherits(pred_mix1_mean, "normMix"))
      assert_that(inherits(pred_mix2_mean, "normMix"))
      
      ## Integration limits 
      lim1 <- qmix(pred_mix1_mean, c(eps / 2, 1 - eps / 2))
      lim2 <- qmix(pred_mix2_mean, c(eps / 2, 1 - eps / 2))
      crit_y1_lower_or_equal_than(lim2, lim1)
      crit_y1_higher_than(lim2, lim1)

      if (n2 == 0) {
        ## Single arm case
        mean_prior2 <- summary(prior2, probs = c())["mean"]
        bound_lower_or_equal_than <- crit_y1_lower_or_equal_than(mean_prior2)
        bound_higher_than         <- crit_y1_higher_than(mean_prior2)
        
        if (bound_lower_or_equal_than <= bound_higher_than) {
          0
        } else {
          pmix(pred_mix1_mean, bound_lower_or_equal_than, lower.tail = TRUE) -
            pmix(pred_mix1_mean, bound_higher_than, lower.tail = TRUE)
        }
      } else {
        ## Double arm case: Integrate over predictive distribution
        integrand <- function(x) {
          bound_lower_or_equal_than <- crit_y1_lower_or_equal_than(x, lim1 = lim1)
          bound_higher_than         <- crit_y1_higher_than(x, lim1 = lim1)
          
          ifelse(
            bound_lower_or_equal_than <= bound_higher_than,
            -Inf,
            log(
              pmix(pred_mix1_mean, bound_lower_or_equal_than, lower.tail = TRUE) -
                pmix(pred_mix1_mean, bound_higher_than, lower.tail = TRUE)
            )
          )
        }
        integrate_density_log(
          integrand,
          pred_mix2_mean,
          logit(eps / 2),
          logit(1 - eps / 2)
        )
      }
    }
  }

  ## 5. WRAPPING AND VECTORIZATION
  ## [Step skipped in pos2S: Vectorization handled outside/not needed]
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##
  ##

  ## 6. RETURN FUNCTION
  design_fun
}


































