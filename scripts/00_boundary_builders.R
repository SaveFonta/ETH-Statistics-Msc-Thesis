# =============================================================================
# Boundary shape builders for the calibrated sequential comparison
# =============================================================================
# A "builder" is a function(p) returning a decisions_list of length K, as
# expected by calibrate_threshold / calibrate_design. The shape fixes the
# relative strictness across looks, and the single scalar p (or u, for the
# Shi & Yin construction) is what the calibration tunes.
#
# Extracted verbatim from scripts/garbage_functions.R so that
# 09_calibrate_boundary_shapes.R has a dependency it can source without side
# effects: both garbage_functions.R and 00_functions_examples.R execute heavy
# top level mclapply sweeps and a ggsave on source, which is fine for scratch
# work but not for a reproducible pipeline. Definitions are unchanged.
# =============================================================================

library(gsDesign)



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

