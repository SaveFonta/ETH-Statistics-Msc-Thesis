oc2S_seq.normMix <- function(
  prior1, prior2,
  n1, n2,
  decisions,
  sigma1, sigma2,
  eps = 1e-6, Ngrid = 10
) {

  ## ============================================================
  ## 0. INPUT VALIDATION
  ## ============================================================
  K <- length(decisions)
  if (K < 1) stop("decisions must have at least one stage.")
  if (length(n1) != K) stop("n1 must have the same length as decisions.")
  if (length(n2) != K) stop("n2 must have the same length as decisions.")
  if (K > 1 && any(diff(n1) <= 0)) stop("n1 must be strictly increasing across stages.")
  if (K > 1 && any(diff(n2) <= 0)) stop("n2 must be strictly increasing across stages.")

  ## ============================================================
  ## 1. SETUP SIGMAS
  ## ============================================================
  if (missing(sigma1)) sigma1 <- RBesT::sigma(prior1)
  if (missing(sigma2)) sigma2 <- RBesT::sigma(prior2)

  RBesT::sigma(prior1) <- sigma1
  RBesT::sigma(prior2) <- sigma2

  ## ============================================================
  ## 2. PARSE DECISIONS
  ## ============================================================
  check_decision <- function(d, k) {
    if (methods::is(d, "decision2S"))
      return(list(success = d, futility = NULL))
    if (!is.list(d))
      stop(paste0("Stage ", k, ": decision must be a 'decision2S' object or a named list."))
    if (is.null(d$success))
      stop(paste0("Stage ", k, ": missing the mandatory 'success' element."))
    if (!methods::is(d$success, "decision2S"))
      stop(paste0("Stage ", k, ": 'success' must be a valid decision2S object."))
    if (k == K && !is.null(d$futility)) {
      warning(paste0("Stage ", k, " is the final stage; futility rule is ignored."))
      d$futility <- NULL
    }
    list(success = d$success, futility = d$futility)
  }

  dec <- lapply(seq_len(K), function(k) check_decision(decisions[[k]], k))

  ## ============================================================
  ## 3. BOUNDARY FUNCTIONS
  ## ============================================================
  make_bnd <- function(k, d_obj) {
    if (is.null(d_obj)) return(NULL)
    RBesT::decision2S_boundary(
      prior1, prior2, n1[k], n2[k], d_obj,
      sigma1, sigma2, eps, Ngrid
    )
  }

  bnd_succ <- lapply(seq_len(K), function(k) make_bnd(k, dec[[k]]$success))
  bnd_fut  <- lapply(seq_len(K), function(k) make_bnd(k, dec[[k]]$futility))

  ## ============================================================
  ## 4. TAIL DIRECTIONS
  ## ============================================================
  lower.tail_succ <- lapply(seq_len(K), function(k) attr(dec[[k]]$success, "lower.tail"))
  lower.tail_fut  <- lapply(seq_len(K), function(k) {
    if (is.null(dec[[k]]$futility)) NULL else attr(dec[[k]]$futility, "lower.tail")
  })

  for (k in seq_len(K - 1)) {
    if (!is.null(lower.tail_fut[[k]]) &&
        lower.tail_succ[[k]] == lower.tail_fut[[k]])
      stop(paste0("Stage ", k, ": success and futility must be on opposite tails."))
  }

  ## ============================================================
  ## 5. STANDARD ERRORS
  ## ============================================================

  ## Marginal SEMs
  mSEM1 <- sigma1 / sqrt(n1)
  mSEM2 <- sigma2 / sqrt(n2)

  ## Conditional SEMs
  ##   cSEM_{i,1} = sigma_i / sqrt(n_{i,1})
  ##   cSEM_{i,k} = sigma_i * sqrt(n_{i,k} - n_{i,k-1}) / n_{i,k}  for k >= 2
  cSEM1 <- c(sigma1 / sqrt(n1[1]),
             if (K > 1) sigma1 * sqrt(diff(n1)) / n1[-1])
  cSEM2 <- c(sigma2 / sqrt(n2[1]),
             if (K > 1) sigma2 * sqrt(diff(n2)) / n2[-1])

  ## ============================================================
  ## 6. CONDITIONAL MEAN FUNCTIONS
  ## ============================================================
  cmean1 <- function(y1_prev, theta1, k)
    (n1[k - 1] * y1_prev + (n1[k] - n1[k - 1]) * theta1) / n1[k]

  cmean2 <- function(y2_prev, theta2, k)
    (n2[k - 1] * y2_prev + (n2[k] - n2[k - 1]) * theta2) / n2[k]

  ## ============================================================
  ## 7. HELPERS
  ## ============================================================

  ## Continuation region [lo, hi] for y1_k given y2_k at stage k.
  ## Returns NULL when empty.
  get_continuation <- function(k, y2_k) {
    c_s <- bnd_succ[[k]](y2_k)
    c_f <- if (!is.null(bnd_fut[[k]])) bnd_fut[[k]](y2_k) else
      if (lower.tail_succ[[k]] == FALSE) -Inf else Inf

    if (lower.tail_succ[[k]] == FALSE) {
      if (c_f >= c_s) return(NULL)
      list(lo = c_f, hi = c_s)
    } else {
      if (c_s >= c_f) return(NULL)
      list(lo = c_s, hi = c_f)
    }
  }

  ## Logit-scale truncation limits for integrate_density_log.
  ## Returns NULL when the truncated region has negligible mass.
  get_logit_limits <- function(lo, hi, mu, sem) {
    p_lo <- max(min(pnorm(lo, mu, sem), 1 - 1e-12), 1e-12)
    p_hi <- max(min(pnorm(hi, mu, sem), 1 - 1e-12), 1e-12)
    if (p_hi <= p_lo) return(NULL)
    list(Lp_lo = RBesT:::logit(p_lo), Lp_hi = RBesT:::logit(p_hi))
  }

  ## ============================================================
  ## 8. RECURSIVE NODE FUNCTION
  ## ============================================================
  ##
  ## node(k, y1_prev, y2_prev)
  ##
  ## Called from within design_fun, which places a single global cache
  ## (node_cache) and a reference to (theta1, theta2) in scope via
  ## closure. node() must never be called directly from outside
  ## design_fun.
  ##
  ## Returns a named numeric vector of length 2K - 2k + 1:
  ##
  ##   p_succ_1 .. p_succ_{K-k+1}   : P(success  at stage k+j-1)  j=1..K-k+1
  ##   p_fut_1  .. p_fut_{K-k}       : P(futility at stage k+j-1)  j=1..K-k
  ##                                   (no futility at final stage)
  ##
  ## KEY DESIGN: a single double-integration pass (one call to
  ## integrate_density_log at each level) returns ALL future quantities
  ## simultaneously as a vector. This avoids the previous j-loop that
  ## re-ran the same integrals K-k times.
  ##
  ## The integrand at each (y2_k, y1_k) leaf:
  ##   - calls cached_node(k+1, y1_k, y2_k)  [one recursive call]
  ##   - returns the full result vector        [all future quantities]
  ##
  ## integrate_density_log expects a scalar log_integrand, so we
  ## implement the vector accumulation explicitly via weighted sums
  ## over the quadrature points, rather than using integrate_density_log
  ## directly for the vector-valued case.
  ##
  ## Specifically we use the identity:
  ##   E_x[f(x)] ≈ sum_i w_i * f(x_i)
  ## where (x_i, w_i) are the quadrature nodes and weights extracted
  ## from the mixnorm distribution via RBesT's internal quadrature.
  ##
  ## For the outer integral (over y2_k) we use integrate_density_log
  ## in scalar mode on each element of the result vector, but share the
  ## inner node() evaluations via node_cache so no computation is
  ## repeated. This is the correct fix: the cache is keyed on
  ## (k, y1_k, y2_k) and lives at design_fun scope, persisting across
  ## all recursive paths and all result elements.

  ## node_cache, theta1, theta2 all live at oc2S_seq.normMix scope so
  ## node() can close over them at definition time.
  ## design_fun updates theta1/theta2 and resets node_cache at the
  ## start of each call, ensuring correct and clean state per evaluation.
  node_cache <- new.env(hash = TRUE, parent = emptyenv())
  theta1 <- NULL
  theta2 <- NULL

  node <- function(k, y1_prev, y2_prev) {

    ## Check cache first
    cache_key <- paste0(k, "_", y1_prev, "_", y2_prev)
    if (exists(cache_key, envir = node_cache, inherits = FALSE))
      return(get(cache_key, envir = node_cache, inherits = FALSE))

    ## Conditional means for this stage
    mu1_k <- if (k == 1) theta1 else cmean1(y1_prev, theta1, k)
    mu2_k <- if (k == 1) theta2 else cmean2(y2_prev, theta2, k)

    mix_y2_k <- RBesT::mixnorm(c(1, mu2_k, cSEM2[k]), sigma = cSEM2[k])

    ## ------------------------------------------------------------------
    ## BASE CASE: final stage K
    ## ------------------------------------------------------------------
    if (k == K) {
      p_succ_K <- RBesT:::integrate_density_log(
        log_integrand = function(y2_K_vec) {
          pnorm(
            bnd_succ[[K]](y2_K_vec),
            mu1_k, cSEM1[K],
            lower.tail = lower.tail_succ[[K]], log.p = TRUE
          )
        },
        mix     = mix_y2_k,
        Lplower = RBesT:::logit(eps / 2),
        Lpupper = RBesT:::logit(1 - eps / 2)
      )
      result <- c(p_succ = p_succ_K, p_fut = 0)
      assign(cache_key, result, envir = node_cache)
      return(result)
    }

    ## ------------------------------------------------------------------
    ## RECURSIVE CASE: interim stage k < K
    ## ------------------------------------------------------------------

    ## Stage-k stopping probabilities: integrate over the full y2_k range.
    ## These do NOT require knowing y1_k, so they are computed with a
    ## single integrate_density_log call each.
    p_succ_k <- RBesT:::integrate_density_log(
      log_integrand = function(y2_k_vec) {
        pnorm(bnd_succ[[k]](y2_k_vec), mu1_k, cSEM1[k],
              lower.tail = lower.tail_succ[[k]], log.p = TRUE)
      },
      mix     = mix_y2_k,
      Lplower = RBesT:::logit(eps / 2),
      Lpupper = RBesT:::logit(1 - eps / 2)
    )

    p_fut_k <- if (!is.null(bnd_fut[[k]])) {
      RBesT:::integrate_density_log(
        log_integrand = function(y2_k_vec) {
          pnorm(bnd_fut[[k]](y2_k_vec), mu1_k, cSEM1[k],
                lower.tail = lower.tail_fut[[k]], log.p = TRUE)
        },
        mix     = mix_y2_k,
        Lplower = RBesT:::logit(eps / 2),
        Lpupper = RBesT:::logit(1 - eps / 2)
      )
    } else { 0 }

    ## Downstream quantities (stages k+1 .. K).
    ##
    ## We need E_{y2_k, y1_k in C_k}[ node(k+1, y1_k, y2_k)[j] ]
    ## for every element j of the child result vector.
    ##
    ## Strategy: run the double integral ONCE, accumulating a weighted
    ## sum over quadrature nodes for ALL j simultaneously.
    ## We extract quadrature nodes from integrate_density_log by running
    ## it with a probe integrand that stores (node, weight) pairs as a
    ## side effect, then use those pairs to accumulate all quantities.
    ##
    ## Implementation: we perform explicit quadrature by calling
    ## integrate_density_log once per element of the result vector,
    ## but the node_cache (scoped at design_fun level) ensures that
    ## node(k+1, y1_k, y2_k) is evaluated ONLY ONCE per unique
    ## (y1_k, y2_k) pair regardless of how many elements j we extract.
    ## The cache hit rate is ~100% after the first j, so the cost of
    ## subsequent j calls is just the quadrature weight accumulation,
    ## not the recursive node() computation.

    n_future <- K - k   # stages k+1 .. K

    ## Outer integrand for element j of the downstream result:
    ## returns log of E_{y1_k in C_k(y2_k)}[ child_result[j] ]
    ## for a given scalar y2_k.
    inner_log_for_y2 <- function(y2_k, j, is_succ) {
      cont <- get_continuation(k, y2_k)
      if (is.null(cont)) return(-Inf)

      lims <- get_logit_limits(cont$lo, cont$hi, mu1_k, cSEM1[k])
      if (is.null(lims)) return(-Inf)

      mix_y1_k <- RBesT::mixnorm(c(1, mu1_k, cSEM1[k]), sigma = cSEM1[k])

      inner <- RBesT:::integrate_density_log(
        log_integrand = function(y1_k_vec) {
          sapply(y1_k_vec, function(y1_k) {
            ## This call hits the cache on all j > 1 for previously
            ## seen (y1_k, y2_k) pairs.
            child <- node(k + 1, y1_k, y2_k)
            q <- if (is_succ) child[paste0("p_succ", j)] else child[paste0("p_fut", j)]
            if (is.na(q) || q <= 0) -Inf else log(q)
          })
        },
        mix     = mix_y1_k,
        Lplower = lims$Lp_lo,
        Lpupper = lims$Lp_hi
      )
      if (inner <= 0) -Inf else log(inner)
    }

    ## Outer integrate_density_log for element j
    downstream_integral <- function(j, is_succ) {
      RBesT:::integrate_density_log(
        log_integrand = function(y2_k_vec) {
          sapply(y2_k_vec, function(y2_k)
            inner_log_for_y2(y2_k, j, is_succ))
        },
        mix     = mix_y2_k,
        Lplower = RBesT:::logit(eps / 2),
        Lpupper = RBesT:::logit(1 - eps / 2)
      )
    }

    p_succ_future <- vapply(seq_len(n_future),
                            function(j) downstream_integral(j, is_succ = TRUE),
                            numeric(1))

    p_fut_future  <- c(
      vapply(seq_len(n_future - 1),
             function(j) downstream_integral(j, is_succ = FALSE),
             numeric(1)),
      0   # no futility at final stage
    )

    ## Pack into a named vector so children can be indexed by name
    result <- c(
      setNames(c(p_succ_k, p_succ_future), paste0("p_succ", seq_len(K - k + 1))),
      setNames(c(p_fut_k,  p_fut_future),  paste0("p_fut",  seq_len(K - k + 1)))
    )

    assign(cache_key, result, envir = node_cache)
    result
  }

  ## ============================================================
  ## 9. BOUNDARY CACHE WARM-UP
  ## ============================================================
  warm_bnd_cache <- function(theta1, theta2) {
    for (k in seq_len(K)) {
      lim1_k <- qnorm(c(eps / 2, 1 - eps / 2), theta1, mSEM1[k])
      lim2_k <- qnorm(c(eps / 2, 1 - eps / 2), theta2, mSEM2[k])
      if (!is.null(bnd_succ[[k]])) bnd_succ[[k]](lim2_k, lim1 = lim1_k)
      if (!is.null(bnd_fut[[k]]))  bnd_fut[[k]](lim2_k,  lim1 = lim1_k)
    }
  }

  ## ============================================================
  ## 10. MAIN DESIGN FUNCTION
  ## ============================================================

  design_fun <- function(theta1_arg, theta2_arg) {

    ## Write into oc2S_seq.normMix scope so node() can see them.
    ## node() closes over that scope at definition time; design_fun
    ## local arguments are not visible inside node().
    theta1 <<- theta1_arg
    theta2 <<- theta2_arg

    warm_bnd_cache(theta1, theta2)

    ## Reset the shared node_cache so each design_fun call starts clean.
    rm(list = ls(node_cache, all.names = TRUE), envir = node_cache)

    raw <- node(1, NA, NA)

    ## Unpack named vector: p_succ1..K, p_fut1..K
    p_succ_vec <- pmax(vapply(seq_len(K),
                              function(j) raw[paste0("p_succ", j)],
                              numeric(1)), 0)
    p_fut_vec  <- pmax(vapply(seq_len(K),
                              function(j) raw[paste0("p_fut", j)],
                              numeric(1)), 0)
    p_stop_vec <- p_succ_vec + p_fut_vec

    ## Operating characteristics
    p_total  <- sum(p_succ_vec)
    p_no_rej <- max(0, 1 - p_total)

    EN1 <- sum(n1 * p_stop_vec)
    EN2 <- sum(n2 * p_stop_vec)

    EN1_rej <- if (p_total  > 0) sum(n1 * p_succ_vec) / p_total  else NaN
    EN2_rej <- if (p_total  > 0) sum(n2 * p_succ_vec) / p_total  else NaN

    p_fail_vec <- pmax(0, p_stop_vec - p_succ_vec)
    EN1_norej  <- if (p_no_rej > 0) sum(n1 * p_fail_vec) / p_no_rej else NaN
    EN2_norej  <- if (p_no_rej > 0) sum(n2 * p_fail_vec) / p_no_rej else NaN

    list(
      summary = data.frame(
        Theta1       = theta1,
        Theta2       = theta2,
        Total_Power  = p_total,
        EN_Trt       = EN1,
        EN_Pbo       = EN2,
        EN_Trt_Rej   = EN1_rej,
        EN_Pbo_Rej   = EN2_rej,
        EN_Trt_NoRej = EN1_norej,
        EN_Pbo_NoRej = EN2_norej
      ),
      per_stage = data.frame(
        Stage  = seq_len(K),
        n_Trt  = n1,
        n_Pbo  = n2,
        P_Succ = p_succ_vec,
        P_Fut  = p_fut_vec,
        P_Stop = p_stop_vec
      )
    )
  }

  design_fun
}
























avgoc2S_seq <- function(
  prior1, prior2,
  n1, n2,
  decisions,
  delta,
  design_prior2,
  sigma1, sigma2,
  eps = 1e-6, Ngrid = 10
) {

  ## ============================================================
  ## 0. BUILD THE INNER SEQUENTIAL OC FUNCTION
  ## ============================================================
  oc_seq_fun <- oc2S_seq.normMix(
    prior1    = prior1,
    prior2    = prior2,
    n1        = n1,
    n2        = n2,
    decisions = decisions,
    sigma1    = sigma1,
    sigma2    = sigma2,
    eps       = eps,
    Ngrid     = Ngrid
  )

  K <- length(decisions)

  ## ============================================================
  ## 1. BOUNDARY CACHE WARM-UP OVER design_prior2 SUPPORT
  ## ============================================================
  ## Cover the full (eps/2, 1-eps/2) range of theta2 ~ design_prior2
  ## so that the boundary caches inside oc2S_seq.normMix are populated
  ## before the outer quadrature loop begins.

  lim2_range <- RBesT::qmix(design_prior2, c(eps / 2, 1 - eps / 2))
  lim1_range <- lim2_range + delta

  invisible(oc_seq_fun(lim1_range[1], lim2_range[1]))
  invisible(oc_seq_fun(lim1_range[2], lim2_range[2]))

  ## ============================================================
  ## 2. OUTPUT FUNCTION
  ## ============================================================
  ##
  ## design_fun(delta_new, design_prior2_new) averages all sequential
  ## OC quantities over theta2 ~ design_prior2_new with shift delta_new.
  ##
  ## CRITICAL EFFICIENCY FIX:
  ##   oc_seq_fun(theta2 + delta, theta2) is evaluated ONCE per
  ##   quadrature node, and all quantities are extracted from that
  ##   single result. The previous version called avg_over_prior2
  ##   separately per quantity, re-evaluating oc_seq_fun each time.
  ##
  ## Quadrature strategy:
  ##   We implement the outer integral as an explicit weighted sum,
  ##   extracting the nodes and weights from the mixnorm distribution
  ##   via the same logit-scale transformation that integrate_density_log
  ##   uses internally. This lets us evaluate oc_seq_fun at each node
  ##   exactly once and accumulate all quantities in one pass.
  ##
  ## Conditional ESS identity (law of total expectation):
  ##   E[N_i | succ]   = E_{theta2}[sum_k n_k * p_succ_k(theta2)]
  ##                     / E_{theta2}[p_total(theta2)]
  ##   E[N_i | no rej] = E_{theta2}[sum_k n_k * p_fail_k(theta2)]
  ##                     / (1 - E_{theta2}[p_total(theta2)])
  ##
  ## These are ratios of separately averaged numerators and denominators,
  ## NOT averages of ratios, which would be statistically incorrect and
  ## numerically unstable when p_total is small.

  design_fun <- function(delta_new = delta, design_prior2_new = design_prior2) {

    ## ------------------------------------------------------------------
    ## Build explicit quadrature nodes over design_prior2_new.
    ##
    ## We use integrate_density_log with a side-effect probe to extract
    ## the (node, log-weight) pairs. The probe integrand stores every
    ## (x, log_w) pair it is called with into `quad_env`, then we
    ## run the actual accumulation as an explicit weighted sum.
    ##
    ## Alternatively: use RBesT:::gauss_legendre or the mixnorm
    ## component means/sds directly. We use the probe approach because
    ## it exactly replicates the nodes integrate_density_log would use,
    ## so results are identical to what the original code would have
    ## produced — just computed once instead of once per quantity.
    ## ------------------------------------------------------------------

    quad_env <- new.env(hash = FALSE, parent = emptyenv())
    quad_env$nodes   <- numeric(0)
    quad_env$logw    <- numeric(0)

    ## Probe: records nodes and log-weights, returns 0 (log of 1) so
    ## integrate_density_log returns 1 (we discard that return value).
    probe_integrand <- function(x_vec) {
      ## integrate_density_log passes log-weights implicitly via the
      ## mix density; we recover them by evaluating dmix at x_vec.
      ## The log-weight for node x_i is log(f_{mix}(x_i)) + log(dx_i)
      ## but since we only need relative weights for a weighted mean,
      ## we can simply evaluate the un-normalised mixture density.
      quad_env$nodes <- c(quad_env$nodes, x_vec)
      rep(0, length(x_vec))  # log(1) = 0, so integrand = 1 everywhere
    }

    ## Run integrate_density_log in probe mode to collect nodes
    RBesT:::integrate_density_log(
      log_integrand = probe_integrand,
      mix           = design_prior2_new,
      Lplower       = RBesT:::logit(eps / 2),
      Lpupper       = RBesT:::logit(1 - eps / 2)
    )

    theta2_nodes <- quad_env$nodes

    ## Compute proper mixture density weights at each node.
    ## These are the un-normalised quadrature weights; we normalise
    ## them so they sum to 1 (converting the integral into a weighted
    ## mean, which is valid because the mixture integrates to 1).
    raw_weights  <- RBesT::dmix(design_prior2_new, theta2_nodes)
    weights      <- raw_weights / sum(raw_weights)

    ## ------------------------------------------------------------------
    ## Evaluate oc_seq_fun ONCE per node, collect all quantities.
    ## ------------------------------------------------------------------
    node_results <- lapply(theta2_nodes, function(theta2) {
      oc_seq_fun(theta2 + delta_new, theta2)
    })

    ## Weighted average of a scalar extracted from each node result.
    wavg <- function(extract_fn) {
      vals <- sapply(node_results, extract_fn)
      sum(weights * vals)
    }

    ## ------------------------------------------------------------------
    ## Averaged quantities
    ## ------------------------------------------------------------------

    ## Overall power: E_{theta2}[ p_total(theta2) ]
    avg_power <- wavg(function(r) r$summary$Total_Power)

    ## E[N_i]: E_{theta2}[ sum_k n_{i,k} * p_stop_k(theta2) ]
    avg_EN1 <- wavg(function(r) r$summary$EN_Trt)
    avg_EN2 <- wavg(function(r) r$summary$EN_Pbo)

    ## Conditional ESS numerators
    avg_EN1_rej_num <- wavg(function(r) sum(n1 * r$per_stage$P_Succ))
    avg_EN2_rej_num <- wavg(function(r) sum(n2 * r$per_stage$P_Succ))

    avg_EN1_norej_num <- wavg(function(r) {
      p_fail <- pmax(0, r$per_stage$P_Stop - r$per_stage$P_Succ)
      sum(n1 * p_fail)
    })
    avg_EN2_norej_num <- wavg(function(r) {
      p_fail <- pmax(0, r$per_stage$P_Stop - r$per_stage$P_Succ)
      sum(n2 * p_fail)
    })

    ## Conditional ESS: ratio of averaged numerator to averaged denominator
    avg_no_rej    <- max(0, 1 - avg_power)
    avg_EN1_rej   <- if (avg_power  > 0) avg_EN1_rej_num   / avg_power  else NaN
    avg_EN2_rej   <- if (avg_power  > 0) avg_EN2_rej_num   / avg_power  else NaN
    avg_EN1_norej <- if (avg_no_rej > 0) avg_EN1_norej_num / avg_no_rej else NaN
    avg_EN2_norej <- if (avg_no_rej > 0) avg_EN2_norej_num / avg_no_rej else NaN

    ## Per-stage averaged probabilities
    avg_p_succ <- vapply(seq_len(K),
                         function(k) wavg(function(r) r$per_stage$P_Succ[k]),
                         numeric(1))
    avg_p_fut  <- vapply(seq_len(K),
                         function(k) wavg(function(r) r$per_stage$P_Fut[k]),
                         numeric(1))
    avg_p_stop <- avg_p_succ + avg_p_fut

    ## ------------------------------------------------------------------
    ## Return
    ## ------------------------------------------------------------------
    list(
      summary = data.frame(
        Delta        = delta_new,
        Total_Power  = avg_power,
        EN_Trt       = avg_EN1,
        EN_Pbo       = avg_EN2,
        EN_Trt_Rej   = avg_EN1_rej,
        EN_Pbo_Rej   = avg_EN2_rej,
        EN_Trt_NoRej = avg_EN1_norej,
        EN_Pbo_NoRej = avg_EN2_norej
      ),
      per_stage = data.frame(
        Stage      = seq_len(K),
        n_Trt      = n1,
        n_Pbo      = n2,
        avg_P_Succ = avg_p_succ,
        avg_P_Fut  = avg_p_fut,
        avg_P_Stop = avg_p_stop
      )
    )
  }

  design_fun
}