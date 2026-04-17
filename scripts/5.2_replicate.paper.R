library(RBesT)



# deifne data
sigma <- 88
data <- crohn
data <- transform(data, y.se =sigma /sqrt(n))


# define possible priors
p_MAP <- mixnorm(c(0.51, -51, 19.9), c(0.44, -46.8, 7.6), c(0.05, -54.1, 51.7), sigma = sigma, param = "ms")
p_MAP



# --------------------------
# Can even do it automatically

map_mcmc <- gMAP(cbind(y, y.se) ~ 1 | study,
  weights = n, data = data,
  family = gaussian,
  beta.prior = cbind(0, sigma),
  tau.dist = "HalfNormal", tau.prior = cbind(0, sigma / 2)
)

map <- automixfit(map_mcmc)
print(map)
plot(map)$mix

round(ess(map)) ## default elir method


# Now lets fit only one Normal
one <- mixfit(map_mcmc, Nc = 1) 
round(ess(one)) #almost 20 patients, not bad, similar to the 20 from Gsponer .



# --------------------------------------

p_rob <- robustify(p_MAP, 0.2, mean = -50)

p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma)



# define success criterion 
succ.crit <- decision2S(pc = 0.975, qc = 0, lower.tail = TRUE)





# --------------------
# CREATE OC
# --------------------
n.act <- 40
n.pbo <- 20

oc_vague <- oc2S(prior1 = p_vague, prior2 = p_vague, n1 = n.act, n2 = n.pbo, decision = succ.crit)
oc_MAP <- oc2S(prior1 = p_vague, prior2 = p_MAP, n1 = n.act, n2 = n.pbo, decision = succ.crit)
oc_rob <- oc2S(prior1 = p_vague, prior2 = p_rob, n1 = n.act, n2 = n.pbo, decision = succ.crit)



# -------------------
# Evalute classic T1E
# -------------------
theta_c <- seq(-150, 50)
theta_a <- theta_c

T1E_vague <- oc_vague(theta_c, theta_a)
T1E_MAP <- oc_MAP(theta_c, theta_a)
T1E_rob <- oc_rob(theta_c, theta_a)

df <- data.frame (theta = theta_c, T1E_vague = T1E_vague, T1E_MAP = T1E_MAP, T1E_rob= T1E_rob )

library(dplyr)
library(tidyr)

df_plot <- pivot_longer(df, cols = c(T1E_vague, T1E_MAP, T1E_rob), names_to = "prior")


library(ggplot2)
ggplot(df_plot, aes(x = theta, y = value, color = prior)) + 
    geom_line(linewidth = 1) + 
    geom_vline (xintercept = -50, linetype = "dashed", alpha = 0.5) + 
    geom_hline (yintercept = 0.035, linetype = "dotted") +
    labs (
        title = "Classic T1E", 
x = "True Mean CDAI Change from Baseline",
    y = "T1E",
    color = "Analysis Prior"
    ) +
    ylim(0,1)




# -----------------------------------
# EVALUATE Average T1E
# -----------------------------------

# WRONG APPROACH:
pos_vague <- pos2S(prior1 = p_vague, prior2 = p_vague, n1 = n.act, n2 = n.pbo, decision = succ.crit)
pos_MAP <- pos2S(prior1 = p_vague, prior2 = p_MAP, n1 = n.act, n2 = n.pbo, decision = succ.crit)
pos_rob <- pos2S(prior1 = p_vague, prior2 = p_rob, n1 = n.act, n2 = n.pbo, decision = succ.crit) 


avg.T1E_vague <- pos_vague(p_vague, p_vague)
avg.T1E_MAP <- pos_MAP(p_vague, p_MAP)
avg.T1E_rob <- pos_rob(p_vague, p_rob)




# no this doesnt work, since like this we are assuming the two priors are completely indepedent 
# we need to sample and then plug those values 


# SIMULATION APPROACH
samples_vague <- rmix(p_vague, n = 10000)
samples_MAP   <- rmix(p_MAP, n = 10000)
samples_rob   <- rmix(p_rob, n = 10000)
 



pointwise_T1E_vague <- oc_vague(samples_vague, samples_vague)
pointwise_T1E_MAP   <- oc_MAP(samples_MAP, samples_MAP)
pointwise_T1E_rob   <- oc_rob(samples_rob, samples_rob)

mean(pointwise_T1E_vague)
mean(pointwise_T1E_MAP)
mean(pointwise_T1E_rob)



# Ok that we can simulate, but we can also define the integral!
#our function is integral[P(reject | theta_c = theta_a) * p(theta_c)] d theta_c



# -----------------------------------
# EVALUATE Average T1E for real
# -----------------------------------
average_T1E <- function(analysis.criteria, design.prior, lower = NA, upper = NA, n_nodes = 128, width_max = 2000) {

    
    # Extract the extremes for integration 
    if (is.na(lower) && is.na(upper)) { 

        # find the x values that cover 99.99 prob mass
        bounds <- qmix(design.prior, p = c(1e-5, 1 - 1e-5))
        lower <- bounds[1]
        upper <- bounds[2]
    }

    width <- upper - lower


    #create the function to integrate 
    to.integrate <- function(x) analysis.criteria(x,x) * dmix(design.prior, x)

    # integrate 
    if (width > width_max) {  # Gauss-Legendre 
    gl <- statmod::gauss.quad(n = n_nodes, kind = "legendre")
    
    # Scale nodes to our window
    mid  <- (upper + lower) / 2
    half.width <- (upper - lower) / 2
    x_nodes <- mid + half.width * gl$nodes
    
    # Weighted Summation
    fvals <- to.integrate(x_nodes)
    res_val <- half.width * sum(gl$weights * fvals)
    
    cat("Gauss-Legendre used with", n_nodes, "nodes since width larger ", width_max)
  }

  else { # classic integrate 
    res <- integrate(to.integrate, lower = lower, upper = upper)
    res_val <- res$value
  }
  
  return(res_val)
}






# Lets recreate the table 2
# we have an additional design prior, which is skeptical 
p_skep <- mixnorm(c(1, -90, 25), sigma = 88)


T1E_vague_row <- c(
  average_T1E(oc_vague, p_vague),
  average_T1E(oc_vague, p_skep),
  average_T1E(oc_vague, p_MAP),
  average_T1E(oc_vague, p_rob)
)


T1E_MAP_row <- c(
  average_T1E(oc_MAP, p_vague),
  average_T1E(oc_MAP, p_skep),
  average_T1E(oc_MAP, p_MAP),
  average_T1E(oc_MAP, p_rob)
)

T1E_rob_row <-  c(
  average_T1E(oc_rob, p_vague),
  average_T1E(oc_rob, p_skep),
  average_T1E(oc_rob, p_MAP),
  average_T1E(oc_rob, p_rob)
)

Table_2 <- rbind(
  T1E_vague_row,
  T1E_MAP_row,
  T1E_rob_row
)
Table_2 <- as.data.frame(Table_2)



colnames(Table_2) <- c("Design_Vague", "Design_Skeptical", "Design_MAP", "Design_Robust")
rownames(Table_2) <- c("Analysis_Vague", "Analysis_MAP", "Analysis_Robust")
Table_2_Formatted <- round(Table_2 * 100, 1)


print(Table_2_Formatted)










library(checkmate)   # assert_number, assert_scalar, assert_list, assert_function
library(assertthat)  # assert_that

avgoc2S.normMix(
  prior1 = p_vague,
  prior2 = p_rob,
  n1 = n.act,
  n2 = n.pbo,
  decision = succ.crit,
  delta = 0,        # <-- theta1 = theta2 + delta
  mix2 = p_skep,         # <-- distribution of theta2 ( design prior)
  sigma1 = sigma,
  sigma2 = sigma
) 


prior1 = p_vague, prior2 = p_MAP, n1 = n.act, n2 = n.pbo, decision = succ.crit




































































#' Averaged sequential Bayesian OC with theta1 = theta2 + delta
#'
#' Computes E_{theta2 ~ mix2}[ OC_seq(theta2 + delta, theta2) ] for a K-stage
#' sequential design with optional early stopping for efficacy and/or futility
#' at interim stages.
#'
#' @param prior1     RBesT normMix prior for arm 1
#' @param prior2     RBesT normMix prior for arm 2
#' @param n1         numeric vector (length K): CUMULATIVE arm-1 sizes per stage
#' @param n2         numeric vector (length K): CUMULATIVE arm-2 sizes per stage
#' @param decisions  list of length K. Each element is either:
#'   - a decision2S object (success criterion; no futility stopping at this stage), OR
#'   - a named list list(success = d_s, futility = d_f) for interim stages with both
#'     efficacy and futility stopping. Only 1-sided decisions are supported.
#'   The final stage (k = K) must specify only a success criterion.
#' @param delta      scalar: theta1 := theta2 + delta
#' @param mix2       normMix distribution over theta2
#' @param sigma1     known SD for arm 1 (defaults to prior reference scale)
#' @param sigma2     known SD for arm 2 (defaults to prior reference scale)
#' @param eps        numerical tail cut-off probability (default 1e-6)
#' @param Ngrid      grid fineness for decision2S_boundary (default 10)
#'
#' @return numeric scalar: the averaged sequential OC
#'
#' @details
#' The algorithm is a recursive numerical integration over stages. At stage k,
#' the probability of success at that stage is a 1-D integral (over y2_k). The
#' probability of eventually succeeding after passing through the continue region
#' requires a 2-D integral (outer over y2_k, inner over y1_k in the continue
#' region), with each inner evaluation recursing into stage k+1.
#'
#' The conditional distribution of the cumulative mean at stage k given stage k-1:
#'   Y1_k | Y1_{k-1} = y, theta1  ~  N( (n1[k-1]*y + dn1*theta1) / n1[k],
#'                                        sigma1^2 * dn1 / n1[k]^2 )
#' where dn1 = n1[k] - n1[k-1]. Arm 2 analogously.
#'
#' Continue region convention (both decisions must be 1-sided):
#'   - success = upper tail (lower.tail=FALSE), futility = lower tail (lower.tail=TRUE):
#'       continue iff  crit_futility(y2) < y1 <= crit_success(y2)
#'   - success = lower tail (lower.tail=TRUE), futility = upper tail (lower.tail=FALSE):
#'       continue iff  crit_success(y2) < y1 <= crit_futility(y2)
#'
#' Setting K=1 with no futility rule recovers the single-stage avgoc2S result.
#' Computational cost scales as O(quadrature_points^K) nested integrations per
#' theta2 quadrature point; practical for K <= 3.
avgoc2S_sequential.normMix <- function(
  prior1,
  prior2,
  n1,
  n2,
  decisions,
  delta,
  mix2,
  sigma1,
  sigma2,
  eps   = 1e-6,
  Ngrid = 10,
  ...
) {
  ## ---- input validation -----------------------------------------------
  K <- length(n1)
  stopifnot(
    K >= 1,
    length(n2) == K,
    length(decisions) == K,
    is.numeric(n1), all(n1 > 0), all(diff(n1) > 0),   # strictly increasing cumulative sizes
    is.numeric(n2), all(n2 >= 0),
    is.numeric(delta), length(delta) == 1,
    inherits(mix2, "normMix")
  )

  if (missing(sigma1)) {
    sigma1 <- RBesT::sigma(prior1)
    message("Using default prior1 reference scale ", sigma1)
  }
  if (missing(sigma2)) {
    sigma2 <- RBesT::sigma(prior2)
    message("Using default prior2 reference scale ", sigma2)
  }
  stopifnot(
    is.numeric(sigma1), length(sigma1) == 1, sigma1 > 0,
    is.numeric(sigma2), length(sigma2) == 1, sigma2 > 0
  )

  sigma(prior1) <- sigma1
  sigma(prior2) <- sigma2

  ## ---- parse stage decisions ------------------------------------------
  ## Normalise each stage entry to list(success = d_s, futility = d_f | NULL)
  parse_decision <- function(d) {
    if (is(d, "decision2S")) return(list(success = d, futility = NULL))
    stopifnot(is.list(d), !is.null(d$success))
    list(success = d$success, futility = d$futility)
  }
  dec <- lapply(decisions, parse_decision)

  if (!is.null(dec[[K]]$futility)) {
    warning("Futility criterion at the final stage (k=K) is not supported and will be ignored.")
    dec[[K]]$futility <- NULL
  }

  ## ---- build boundary functions per stage ----------------------------
  make_bnd <- function(k, d_obj) {
    if (is.null(d_obj)) return(NULL)
    decision2S_boundary(prior1, prior2, n1[k], n2[k], d_obj,
                        sigma1, sigma2, eps, Ngrid)
  }
  bnd_success  <- lapply(seq_len(K), function(k) make_bnd(k, dec[[k]]$success))
  bnd_futility <- lapply(seq_len(K), function(k) make_bnd(k, dec[[k]]$futility))

  ## Extract lower.tail from each 1-sided decision (NA if stage has none)
  get_lt <- function(d_obj) {
    if (is.null(d_obj)) return(NA)
    stopifnot(is(d_obj, "decision2S_1sided"))
    attr(d_obj, "lower.tail")
  }
  lt_success  <- sapply(dec, function(d) get_lt(d$success))
  lt_futility <- sapply(dec, function(d) get_lt(d$futility))

  ## ---- conditional distribution helpers ------------------------------
  ## Pre-compute the conditional SD of the cumulative mean at stage k
  ## given stage k-1 (independent-increments property of normal data).
  ##
  ##   sd(Y1_k | Y1_{k-1}, theta1) = sigma1 * sqrt(n1[k] - n1[k-1]) / n1[k]
  ##   sd(Y1_1 | theta1)           = sigma1 / sqrt(n1[1])
  csd1 <- numeric(K)
  csd2 <- numeric(K)
  csd1[1] <- sigma1 / sqrt(n1[1])
  csd2[1] <- sigma2 / sqrt(n2[1])
  if (K > 1) for (k in 2:K) {
    csd1[k] <- sigma1 * sqrt(n1[k] - n1[k-1]) / n1[k]
    csd2[k] <- sigma2 * sqrt(n2[k] - n2[k-1]) / n2[k]
  }

  ## Conditional mean of cumulative arm-1 mean at stage k
  cmean1 <- function(k, y1p, theta1) {
    if (k == 1L) return(theta1)
    (n1[k-1] * y1p + (n1[k] - n1[k-1]) * theta1) / n1[k]
  }
  cmean2 <- function(k, y2p, theta2) {
    if (k == 1L) return(theta2)
    (n2[k-1] * y2p + (n2[k] - n2[k-1]) * theta2) / n2[k]
  }

  ## ---- continue-region helper ----------------------------------------
  ## Returns list(lower, upper): the open interval for Y1_k that triggers
  ## "continue" at stage k, given Y2_k = y2.
  ##
  ## Convention: success and futility boundaries must have opposite lower.tail
  ## directions (one upper-tail, one lower-tail), which covers all practical designs.
  continue_limits <- function(k, y2) {
    c_s  <- bnd_success[[k]](y2)
    lt_s <- lt_success[k]

    if (is.null(bnd_futility[[k]])) {
      ## No futility stopping: continue = complement of success
      if (!lt_s) return(list(lower = -Inf, upper = c_s))  # success = Y1 > c_s
      else       return(list(lower =  c_s, upper =  Inf)) # success = Y1 <= c_s (unusual)
    }

    c_f  <- bnd_futility[[k]](y2)
    lt_f <- lt_futility[k]

    ## Typical: success = upper tail (lt_s=F), futility = lower tail (lt_f=T)
    ##   continue iff  c_f < Y1 <= c_s
    if (!lt_s &&  lt_f) return(list(lower = c_f, upper = c_s))

    ## Mirror: success = lower tail (lt_s=T), futility = upper tail (lt_f=F)
    ##   continue iff  c_s < Y1 <= c_f
    if ( lt_s && !lt_f) return(list(lower = c_s, upper = c_f))

    stop(
      "success and futility decisions at stage ", k,
      " have the same lower.tail direction; cannot determine the continue region."
    )
  }

  ## ---- recursive stage-OC computation --------------------------------
  ##
  ## compute_from_stage(k, y1p, y2p, theta1, theta2)
  ##
  ## Returns P(eventual success | we have reached stage k with cumulative
  ## running means Y1_{k-1} = y1p, Y2_{k-1} = y2p, and true means theta1, theta2).
  ## For k=1 the previous means y1p, y2p are unused (pass NA_real_).
  ##
  ## Structure at each stage k:
  ##   P(eventual success)
  ##     = P(success at k)                        [1-D integral over y2_k]
  ##     + P(continue at k) * P(success later)    [2-D integral: outer y2_k, inner y1_k]
  ##
  ## The 2-D integral has the inner integrand as a recursive call to stage k+1,
  ## giving a naturally nested quadrature that terminates at the final stage.
  compute_from_stage <- function(k, y1p, y2p, theta1, theta2) {
    mu1 <- cmean1(k, y1p, theta1)
    mu2 <- cmean2(k, y2p, theta2)
    sd1 <- csd1[k]
    sd2 <- csd2[k]

    y2_lo <- qnorm(eps / 2,   mu2, sd2)
    y2_hi <- qnorm(1-eps / 2, mu2, sd2)
    y1_lo <- qnorm(eps / 2,   mu1, sd1)
    y1_hi <- qnorm(1-eps / 2, mu1, sd1)

    ## Warm up boundary caches over the integration range for this call
    if (is(dec[[k]]$success,  "decision2S_1sided"))
      bnd_success[[k]](c(y2_lo, y2_hi),  lim1 = c(y1_lo, y1_hi))
    if (!is.null(bnd_futility[[k]]) && is(dec[[k]]$futility, "decision2S_1sided"))
      bnd_futility[[k]](c(y2_lo, y2_hi), lim1 = c(y1_lo, y1_hi))

    lt_s <- lt_success[k]
    tol  <- eps^0.4   # relaxed per nesting level; tighter hurts nested performance

    ## --- P(success at stage k) ----------------------------------------
    ## = integral_{y2} P(Y1_k in success region | mu1, sd1) * phi(y2; mu2, sd2) dy2
    p_success <- integrate(
      function(y2_vec) {
        vapply(y2_vec, function(y2) {
          pnorm(bnd_success[[k]](y2), mu1, sd1, lower.tail = lt_s) *
            dnorm(y2, mu2, sd2)
        }, FUN.VALUE = 0.0)
      },
      lower = y2_lo, upper = y2_hi, rel.tol = tol
    )$value

    ## Final stage: no future stages to continue into
    if (k == K) return(p_success)

    ## --- P(continue at k then eventually succeed) ---------------------
    ## = integral_{y2} [ integral_{y1 in continue(y2)} 
    ##     compute_from_stage(k+1, y1, y2, ...) * phi(y1; mu1, sd1) dy1 ]
    ##   * phi(y2; mu2, sd2) dy2
    p_future <- integrate(
      function(y2_vec) {
        vapply(y2_vec, function(y2) {
          lims     <- continue_limits(k, y2)
          y1_lower <- max(lims$lower, y1_lo)   # clip to numerical support
          y1_upper <- min(lims$upper, y1_hi)

          if (y1_lower >= y1_upper) return(0.0) # continue region has negligible probability

          inner <- tryCatch(
            integrate(
              function(y1_vec) {
                ## recursive call: each y1 triggers a full stage-(k+1) computation
                vapply(y1_vec, function(y1) {
                  compute_from_stage(k + 1L, y1, y2, theta1, theta2)
                }, FUN.VALUE = 0.0) *
                  dnorm(y1_vec, mu1, sd1)
              },
              lower = y1_lower, upper = y1_upper, rel.tol = tol
            )$value,
            error = function(e) 0.0   # integration failure treated as 0
          )

          inner * dnorm(y2, mu2, sd2)
        }, FUN.VALUE = 0.0)
      },
      lower = y2_lo, upper = y2_hi, rel.tol = tol
    )$value

    p_success + p_future
  }

  ## ---- fixed-(theta1, theta2) wrapper --------------------------------
  seq_oc_fixed <- function(theta1, theta2) {
    compute_from_stage(1L, NA_real_, NA_real_, theta1, theta2)
  }

  ## ---- warm up boundary caches over the full theta2 support ----------
  ## This seeds the spline caches for the range that the outer theta2
  ## integration will traverse, avoiding repeated cold-start rebuilds.
  lim2_mix  <- qmix(mix2, c(eps / 2, 1 - eps / 2))
  sem1_all  <- sigma1 / sqrt(n1)
  sem2_all  <- sigma2 / sqrt(n2)
  for (k in seq_len(K)) {
    lim2_k <- c(qnorm(eps / 2,   lim2_mix[1],         sem2_all[k]),
                qnorm(1-eps / 2, lim2_mix[2],         sem2_all[k]))
    lim1_k <- c(qnorm(eps / 2,   lim2_mix[1] + delta, sem1_all[k]),
                qnorm(1-eps / 2, lim2_mix[2] + delta, sem1_all[k]))
    if (is(dec[[k]]$success, "decision2S_1sided"))
      bnd_success[[k]](lim2_k,  lim1 = lim1_k)
    if (!is.null(bnd_futility[[k]]) && is(dec[[k]]$futility, "decision2S_1sided"))
      bnd_futility[[k]](lim2_k, lim1 = lim1_k)
  }

  ## ---- outer average: E_{theta2 ~ mix2}[ seq_oc(theta2+delta, theta2) ]
  integrate_density(
    integrand = function(x) vapply(x, function(th2) seq_oc_fixed(th2 + delta, th2), 0.0),
    mix       = mix2,
    Lplower   = logit(eps / 2),
    Lpupper   = logit(1 - eps / 2)
  )
}




# Two-stage design: interim at n=50, final at n=100
# Interim: stop early for efficacy OR futility
# Final: success/failure only

d_efficacy <- RBesT::decision2S(0.95, 0,  lower.tail = FALSE)  # P(θ₁>θ₂|data) > 0.95
d_futility <- RBesT::decision2S(0.10, 0,  lower.tail = TRUE)   # P(θ₁>θ₂|data) < 0.10
d_final    <- RBesT::decision2S(0.90, 0,  lower.tail = FALSE)

avgoc2S_sequential.normMix(
  prior1      = p_vague,
  prior2      = p_MAP,
  n1          = c(50, 100),
  n2          = c(50, 100),
  decisions   = list(
    list(success = d_efficacy, futility = d_futility),   # interim
    d_final                                              # final
  ),
  delta       = 0,        # null scenario: T1E
  mix2        = p_MAP,
  sigma1      = 88,
  sigma2      = 88
)
