# =============================================================================
# Usage examples for 00_functions.R
# =============================================================================
# Nothing here runs on source() (everything costly is wrapped in if (FALSE)),
# this file is meant to be read and copy-pasted from, or run block by block.
# Run from the Thesis/ root so the source() below resolves.
#
# What's here:
#   1. Core engine walkthrough: conditional OC, average OC, sequential EUII
#   2. Calibration walkthrough: default threshold, then a custom decisions_builder
#   3. Ad-hoc custom criteria (mixed SC/DC across stages)
#   4. Power-family spending style
#   5. O'Brien-Fleming, whole shape scaled by one parameter p
#   6. Haybittle-Peto, whole shape scaled by one parameter p
#   7. Shi & Yin (2019) style: keep the *nominal* frequentist shape fixed and
#      calibrate a separate interpolation parameter u instead of p directly
#   8. Linking the futility threshold to the efficacy threshold
#   9. Reference-only: comparing raw boundary Z-scores across shapes via gsDesign
#
# Sections 3-9 are about `decisions_builder`: calibrate_threshold() and
# calibrate_design() (in 00_functions.R) by default calibrate a single flat
# threshold p, applied identically at every look (Pocock style). Passing a
# `decisions_builder = function(p) {...}` overrides that default: only p is
# calibrated, but the builder decides how p turns into the actual per-stage
# decisions_list, so any boundary shape (or any custom multi-criterion rule)
# can be calibrated the same way. A decisions_builder must be function(p)
# returning a list of length K (one entry per stage), each entry
# list(success = <decision2S>, futility = <decision2S or NULL>).
#
# None of sections 3-9 is wired into the thesis pipeline yet (see manuscript's
# Future Work: boundary shapes), it is scaffolding for that extension.
# =============================================================================

source("scripts/00_functions.R")   # oc2_seq_mc.normMix, avgoc2_seq_mc.normMix,
                                    # compute_euii, calibrate_threshold, calibrate_design
library(gsDesign)


# =============================================================================
# 1. Core engine walkthrough: conditional OC, average OC, sequential EUII
# =============================================================================
# n_sim is kept small (1000) here since this is just a walkthrough of the
# calling convention, not a real evaluation (the thesis pipeline uses 1e6).
if (FALSE) {
  library(RBesT)
  sigma   <- 88
  p_MAP   <- mixnorm(
    c(0.4848, -52.457, 21.154), c(0.4598, -47.465, 7.843), c(0.0554, -50.355, 48.164),
    sigma = sigma, param = "ms"
  )
  p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")

  # Efficacy: dual criterion carried by a single decision2S (vector pc/qc).
  sign.crit <- decision2S(pc = c(0.95, 0.50), qc = c(0, -50), lower.tail = TRUE)
  # Futility: Pr(delta < 40 | data) >= 0.90, i.e. Pr(theta_T - theta_C > -40) >= 0.90.
  fut.crit  <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)

  decisions_fut <- list(
    list(success = sign.crit, futility = fut.crit),  # stage 1: efficacy + futility
    list(success = sign.crit, futility = NULL)        # stage 2: final, efficacy only
  )

  n1_seq <- c(30, 60)   # treatment arm, cumulative
  n2_seq <- c(20, 40)   # control arm, cumulative

  # Conditional OC: a scalar theta fixes theta_C.
  res_cond <- oc2_seq_mc.normMix(
    theta_1 = -50, theta_2 = -50,
    prior_1 = p_vague, prior_2 = p_MAP,
    n1_seq  = n1_seq,  n2_seq  = n2_seq,
    decisions_list = decisions_fut,
    n_sim = 1000, seed = 1
  )
  res_cond$Overall
  res_cond$Per_Stage

  # Average OC: theta_C is drawn from the design prior, one draw per simulated
  # trial. delta = 0 must be in the grid, since compute_euii reads the average
  # T1E from it.
  res_avg <- avgoc2_seq_mc.normMix(
    prior_1 = p_vague, prior_2 = p_MAP,
    n1_seq  = n1_seq,  n2_seq  = n2_seq,
    decisions_list = decisions_fut,
    delta          = c(0, 30, 60),
    design_prior_c = p_MAP,
    n_sim = 1000, seed = 123
  )

  # avgoc2_seq_mc.normMix returns a list keyed by delta; each entry has
  # $Overall and $Per_Stage, same shape as oc2_seq_mc.normMix's own output.
  res_avg$delta.0

  # Sequential EUII. prior_H1 enters only here, not in the simulation itself.
  res_euii <- compute_euii(res_avg, prior_H1 = c(0.01, 0.1, 0.5))
  res_euii[["0.01"]]



  # Sensitivity checks: how much do avgT1E/avgPower actually move?
  #   (a) vary p, keep the schedule (n1_seq/n2_seq) fixed -> avgT1E(p), avgPower(p)
  #   (b) vary n2max, keep p fixed -> avgT1E(n2max)

  # --- (a) avgT1E(p) and avgPower(p), schedule fixed -------------------------
  pis <- seq(0.94, 0.99, by = 0.01)

  res_p <- lapply(pis, function(p) {
    crit <- decision2S(pc = p, qc = 0, lower.tail = TRUE)
    dec  <- lapply(seq_along(n1_seq), function(k) list(success = crit, futility = NULL))
    avgoc2_seq_mc.normMix(
      prior_1 = p_vague, prior_2 = p_MAP,
      n1_seq  = n1_seq,  n2_seq  = n2_seq,
      decisions_list = dec,
      delta          = c(0, 60),   
      design_prior_c = p_MAP,
      n_sim = 2e4, seed = 123
    )
  })

  data.frame(
    p        = pis,
    avgT1E   = sapply(res_p, function(r) r$delta.0$Overall[["Power"]]),
    avgPower = sapply(res_p, function(r) r$delta.60$Overall[["Power"]])
  )

  # --- (b) avgT1E(n2max), p fixed ---------------------------------------------
  p_fixed <- 0.95

  allo_ratio <- n1_seq[length(n1_seq)] / n2_seq[length(n2_seq)]
  info_ratio <- n2_seq / n2_seq[length(n2_seq)]
  schedule <- function(n2max) {
    n2 <- pmax(1, round(info_ratio * n2max)); n2 <- cummax(n2)
    n1 <- pmax(1, round(allo_ratio * n2))
    list(n1 = n1, n2 = n2)
  }

  n2max_grid <- seq(20, 100, by = 10)
  crit_fixed <- decision2S(pc = p_fixed, qc = 0, lower.tail = TRUE)

  res_n <- lapply(n2max_grid, function(n2max) {
    ns  <- schedule(n2max)
    dec <- lapply(seq_along(ns$n1), function(k) list(success = crit_fixed, futility = NULL))
    avgoc2_seq_mc.normMix(
      prior_1 = p_vague, prior_2 = p_MAP,
      n1_seq  = ns$n1,   n2_seq  = ns$n2,
      decisions_list = dec,
      delta          = c(0, 60),     # scalar delta -> result returned unwrapped, no $delta.0
      design_prior_c = p_MAP,
      n_sim = 2e4, seed = 123
    )
  })

  data.frame(
    n2max  = n2max_grid,
    avgT1E   = sapply(res_n, function(r) r$delta.0$Overall[["Power"]]),
    avgPower = sapply(res_n, function(r) r$delta.60$Overall[["Power"]])
  )

  res_p 
  res_n
  #as p varies, avgT1E varies a lot, and also avgPower varies a bit
  # if instead we keep p fixed and vary sample size, avgT1E fluctuates and avgPower varies a lot. 

  # This mean that changin n_max doesnt influence much 
}


# =============================================================================
# 2. Calibration walkthrough: default threshold, then a custom decisions_builder
# =============================================================================
if (FALSE) {
  library(RBesT)
  sigma   <- 88
  p_MAP   <- mixnorm(
    c(0.4848, -52.457, 21.154), c(0.4598, -47.465, 7.843), c(0.0554, -50.355, 48.164),
    sigma = sigma, param = "ms"
  )
  p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")

  n1_seq <- c(30, 60)   # treatment arm, cumulative
  n2_seq <- c(20, 40)   # control arm, cumulative

  # Default: calibrate a single flat threshold p (Pocock style) to hit avgT1E = 5%.
  resT1E <- calibrate_threshold(0.05, n1_seq = n1_seq, n2_seq = n2_seq,
                                prior_cntrl = p_MAP, prior_treat = p_vague,
                                criterion = "SC")

  sign.crit <- decision2S(pc = resT1E$p, qc = 0, lower.tail = TRUE)


    decisions_fut <- list(
    list(success = sign.crit, futility = NULL), 
    list(success = sign.crit, futility = NULL)        
  )

  res_avg <- avgoc2_seq_mc.normMix(
    prior_1 = p_vague, prior_2 = p_MAP,
    n1_seq  = n1_seq,  n2_seq  = n2_seq,
    decisions_list = decisions_fut,
    delta          = c(0, 30, 60),
    design_prior_c = p_MAP,
    n_sim = 1e6, seed = 123
  )

  # so the next function use avgoc2_seq_mc.normMix a lot of time, everytime, if the avgPower is too
  # low, it inspect new values of n, first with that specific n, it optimizes to find the specific p,
  # thn looks at avgPower and see if it hits the target. 

  # the search works by finding upper bound and lower bound of n2max such that avgPower is there. Then
  # with bisection we find the smallest integer s.t. avgPower > targetpower. 

  # Of course there is a schedule function that gien n2max, keep constant the information rati and the allocation ratio 

 # for example now, we should increase the n2max, then recalibrate p and find avgPower 


  # Joint calibration: smallest sample size that also reaches 90% avgPower at delta = 60, with the threshold
  # p recalibrated at every candidate size.
  res_design <- calibrate_design(target_t1e = 0.05, target_power = 0.9, delta_power = 60,
                                 n1_base = n1_seq, n2_base = n2_seq,
                                 prior_cntrl = p_MAP, prior_treat = p_vague,
                                 criterion = "SC", seed = 123)

## what abou varying the seed? 

seeds <- seq(123, 135)
cores <- if (.Platform$OS.type != "unix") 1L else length(seeds)  # one core per seed

res_seeds <- parallel::mclapply(seq_along(seeds), function(x){
    res_design <- calibrate_design(target_t1e = 0.05, target_power = 0.9, delta_power = 60,
                                 n1_base = n1_seq, n2_base = n2_seq,
                                 prior_cntrl = p_MAP, prior_treat = p_vague,
                                 criterion = "SC", seed = seeds[[x]])
}, mc.cores = cores)

final <- data.frame(
  seed           = seeds,
  p              = sapply(res_seeds, function(r) r$p),
  n_max          = sapply(res_seeds, function(r) max(r$n1_seq) + max(r$n2_seq)),
  achieved_t1e   = sapply(res_seeds, function(r) r$achieved_t1e),
  achieved_power = sapply(res_seeds, function(r) r$achieved_power)
)
print(final)

mean(final$p)





  # What about using a special decision criterion instead of the default flat
  # threshold? Let's say we want the power-family shape (see section 4 below).
  builder_pow <- function(p) {
    p_int <- 1 - (1 - p) * (1 / 2)^2
    list(
      list(success = decision2S(pc = p_int, qc = 0, lower.tail = TRUE), futility = NULL),
      list(success = decision2S(pc = p,     qc = 0, lower.tail = TRUE), futility = NULL)
    )
  }

  resT1E_special <- calibrate_threshold(0.05, n1_seq = n1_seq, n2_seq = n2_seq,
                                        prior_cntrl = p_MAP, prior_treat = p_vague,
                                        criterion = "SC", decisions_builder = builder_pow)

  res_design_special <- calibrate_design(target_t1e = 0.05, target_power = 0.9, delta_power = 60,
                                         n1_base = n1_seq, n2_base = n2_seq,
                                         prior_cntrl = p_MAP, prior_treat = p_vague,
                                         criterion = "SC", decisions_builder = builder_pow)
}


# =============================================================================
# 3. Ad-hoc custom criteria
# =============================================================================

# SC with futility stopping at the interim only.
builder_sc_fut <- function(p) {
  succ <- decision2S(pc = p,    qc = 0,   lower.tail = TRUE)
  fut  <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)
  list(
    list(success = succ, futility = fut),
    list(success = succ, futility = NULL)
  )
}

# SC at the interim, DC at the final look.
builder_sc_then_dc <- function(p) {
  sc <- decision2S(pc = p,         qc = 0,         lower.tail = TRUE)
  dc <- decision2S(pc = c(p, 0.5), qc = c(0, -50), lower.tail = TRUE)
  list(
    list(success = sc, futility = NULL),
    list(success = dc, futility = NULL)
  )
}


# =============================================================================
# 4. Power-family spending style ((Kim-DeMets, 1987))
# =============================================================================
# f(t) = alpha * t^rho, applied here as a shrinkage of the calibrated p at
# earlier looks: p_k = 1 - (1 - p) * t_k^rho, where t_k is the information
# fraction at look k (t_K = 1 at the final look, which therefore always uses
# exactly p). `fracs` are the information fractions, last one = 1.

make_builder_pow <- function(fracs, rho = 2) {
  function(p) {
    p_k <- 1 - (1 - p) * fracs^rho
    lapply(p_k, function(pk) list(
      success  = RBesT::decision2S(pc = pk, qc = 0, lower.tail = TRUE),
      futility = NULL
    ))
  }
}

builder_pow_1int <- make_builder_pow(c(1 / 2, 1))         # one interim
builder_pow_2int <- make_builder_pow(c(1 / 3, 2 / 3, 1))  # two interims


# =============================================================================
# 5. O'Brien-Fleming, whole shape scaled by one parameter p
# =============================================================================
# IN the freq case: 
# Name Z_OF a constant on the Z score scale and c_k the treshold at look k
# then the formula is c_k = Z_OF / sqrt(t_k )where t_k is the information fraction at
# look k 

# Since we look at p, we need to trasform p on the Z scale, then divide by t_k and then 
# go back to p scale:
# p_k = Phi(Phi^-1(p) / sqrt(t_k)), (t_K = 1 at the final look, which therefore always uses exactly p).


make_builder_obf <- function(fracs) {
  function(p) {
    p_k <- pnorm(qnorm(p) / sqrt(fracs))
    lapply(p_k, function(pk) list(
      success  = RBesT::decision2S(pc = pk, qc = 0, lower.tail = TRUE),
      futility = NULL
    ))
  }
}

builder_obf_1int <- make_builder_obf(c(1 / 2, 1))               # one interim
builder_obf_2int <- make_builder_obf(c(1 / 3, 2 / 3, 1))        # two interims
builder_obf_3int <- make_builder_obf(c(1 / 4, 2 / 4, 3 / 4, 1)) # three interims


# =============================================================================
# 6. Haybittle-Peto, whole shape scaled by one parameter p
# =============================================================================
# Near-impossible threshold at every interim
# The final look then uses the standard calibrated p directly.
# `K` is the total number of
# looks (interims + final).

make_builder_hp <- function(K, p_interim = 0.999) {
  function(p) {
    c(
      lapply(seq_len(K - 1), function(k) list(
        success  = RBesT::decision2S(pc = p_interim, qc = 0, lower.tail = TRUE),
        futility = NULL
      )),
      list(list(success = RBesT::decision2S(pc = p, qc = 0, lower.tail = TRUE), futility = NULL))  # final look: standard target
    )
  }
}

builder_hp_1int <- make_builder_hp(2)   # one interim
builder_hp_2int <- make_builder_hp(3)   # two interims


# =============================================================================
# 7. Shi & Yin (2019) style: fixed nominal shape, calibrate an interpolation u
# =============================================================================
# Start from a classical frequentist group-sequential boundary (OBF or
# Pocock, computed once from gsDesign at some nominal alpha), then instead of
# recalibrating the whole shape (sections 5-6), calibrate a single
# interpolation parameter u that slides every look's threshold between that
# nominal boundary and certainty:
#   p_k(u) = Phi(z_k) + (1 - Phi(z_k)) * u
# u = 0 reproduces the raw nominal boundary exactly (already fairly strict,
# e.g. ~0.98-0.997 for a two-look 0.025 design), u = 1 gives threshold = 1
# (impossible). The meaningful range is therefore [0, 1], NOT the default
# p_interval = c(0.6, 0.999) used elsewhere in this file: passing one of
# these builders to calibrate_threshold/calibrate_design without overriding
# p_interval to something like c(0, 1) will search entirely inside the
# "already near-impossible" region and likely misreport the target as
# unattainable. Widen further below 0 if even u = 0 turns out too strict.

# `fracs`: information fractions, last one = 1. `alpha`: the nominal
# frequentist one-sided alpha the raw (uncalibrated) shape targets.
make_builder_sy_obf <- function(fracs, alpha = 0.025) {
  k <- length(fracs)

  # Exact frequentist Z-scores for an O'Brien-Fleming design (one-sided
  # superior, exact "OF" boundary, not the Lan-DeMets approximation).
  gs_obj   <- gsDesign::gsDesign(k = k, test.type = 1, alpha = alpha, timing = fracs, sfu = "OF")
  z_scores <- gs_obj$upper$bound

  p_theoretical <- pnorm(z_scores)

  function(u) {
    # Shi & Yin calibration formula: c_k = Phi(z_k) + {1 - Phi(z_k)} * u
    p_calibrated <- p_theoretical + (1 - p_theoretical) * u
    lapply(p_calibrated, function(pk) {
      list(success = RBesT::decision2S(pc = pk, qc = 0, lower.tail = TRUE), futility = NULL)
    })
  }
}

builder_sy_1int <- make_builder_sy_obf(c(1 / 2, 1))         # 1 interim at 50%, final at 100%
builder_sy_2int <- make_builder_sy_obf(c(1 / 3, 2 / 3, 1))  # 2 interims at 33%/67%, final at 100%

# Same construction, Pocock nominal shape instead of O'Brien-Fleming (flat
# nominal Z-score across looks instead of a rising one).
make_builder_sy_pocock <- function(fracs, alpha = 0.025) {
  k <- length(fracs)

  gs_obj   <- gsDesign::gsDesign(k = k, test.type = 1, alpha = alpha, timing = fracs, sfu = "Pocock")
  z_scores <- gs_obj$upper$bound

  p_theoretical <- pnorm(z_scores)

  function(u) {
    p_calibrated <- p_theoretical + (1 - p_theoretical) * u
    lapply(p_calibrated, function(pk) {
      list(success = RBesT::decision2S(pc = pk, qc = 0, lower.tail = TRUE), futility = NULL)
    })
  }
}

builder_sy_pocock_2int <- make_builder_sy_pocock(c(1 / 3, 2 / 3, 1))  # 2 interims at 33%/67%, final at 100%

# IMPORTANT: for any other classical shape, write a builder the same way
# (fixed nominal Z-scores -> Shi & Yin interpolation), so it's easy to adapt
# to as many looks as needed.


# =============================================================================
# 8. Linking the futility threshold to the efficacy threshold
# =============================================================================
# Not a bad idea: a stricter (larger) efficacy threshold p also loosens the
# futility rule somewhat, capped to stay in a sensible range.

builder_linked <- function(p) {
  succ <- decision2S(pc = p, qc = 0, lower.tail = TRUE)
  fut  <- decision2S(pc = min(0.995, max(0.95, p + 0.1)), qc = -40, lower.tail = FALSE)  # capped, p + 0.1 can pass 1
  list(
    list(success = succ, futility = fut),
    list(success = succ, futility = NULL)
  )
}


# =============================================================================
# 9. Reference only: raw boundary Z-scores side by side, not tied to
#    calibrate_threshold/calibrate_design at all, just for comparing shapes.
# =============================================================================
if (FALSE) {
  k_stages <- 4                       # 4 looks (3 interims + 1 final)
  t_frac   <- c(0.25, 0.5, 0.75, 1)   # information fractions
  alpha    <- 0.05

  # Approx O'Brien-Fleming (Lan-DeMets): 1 - Phi(Phi^-1(alpha) / sqrt(t))
  obf <- gsDesign(k = k_stages, test.type = 1, alpha = alpha, timing = t_frac, sfu = sfLDOF)

  # Exact O'Brien-Fleming
  obf_true <- gsDesign(k = k_stages, test.type = 1, alpha = alpha, timing = t_frac, sfu = "OF")

  # Power family (rho = 2)
  pow2 <- gsDesign(k = k_stages, test.type = 1, alpha = alpha, timing = t_frac, sfu = sfPower, sfupar = 2)

  # Approx Pocock (Lan-DeMets): f = alpha * log(1 + (e-1)*t)
  poc <- gsDesign(k = k_stages, test.type = 1, alpha = alpha, timing = t_frac, sfu = sfLDPocock)
  # equivalently: sfLDPocock(alpha = 0.025, t = 1:4 / 4)$spend

  # Exact Pocock
  poc_true <- gsDesign(k = k_stages, test.type = 1, alpha = alpha, timing = t_frac, sfu = "Pocock")

  # Haybittle-Peto (brick wall at interims, standard at the end).
  # p = 0.001 translates to a Z-score of ~3.09.
  hp <- c(qnorm(1 - 0.001), qnorm(1 - 0.001), qnorm(1 - 0.001), qnorm(1 - alpha))

  z_bounds <- data.frame(
    Fraction       = t_frac,
    OBF            = obf$upper$bound,
    True_OBF       = obf_true$upper$bound,
    Power_R2       = pow2$upper$bound,
    Pocock         = poc$upper$bound,
    True_Pocock    = poc_true$upper$bound,
    Haybittle_Peto = hp
  )

  print(round(z_bounds, 3))
}
