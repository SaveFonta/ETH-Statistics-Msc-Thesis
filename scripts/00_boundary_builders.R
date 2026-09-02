# =============================================================================
# Boundary shape builders for the calibrated sequential comparison
# =============================================================================
# A "builder" is a function(p) returning a decisions_list of length K, as
# expected by calibrate_threshold / calibrate_design. 
# =============================================================================



# =============================================================================
# Wang Tsiatis (Wang and Tsiatis, 1987)
# =============================================================================
# the critical value at
# look k is proportional to t_k^(Delta - 1/2) on the test statistic scale.
# Carried onto the posterior probability scale by the same mapping used for
# O'Brien-Fleming: p_k = Phi(Phi^-1(p) * t_k^(Delta - 1/2)), where t_k is the
# information fraction at look k (t_K = 1 at the final look, which therefore
# always uses exactly p). Delta = 1/2 recovers the flat (Pocock) shape and
# Delta = 0 recovers O'Brien-Fleming exactly. `fracs` are the information
# fractions, last one = 1.

make_builder_wt <- function(fracs, Delta = 0.25) {
  function(p) {
    p_k <- pnorm(qnorm(p) * fracs^(Delta - 1/2))
    lapply(p_k, function(pk) list(
      success  = RBesT::decision2S(pc = pk, qc = 0, lower.tail = TRUE),
      futility = NULL
    ))
  }
}

builder_wt_1int <- make_builder_wt(c(1 / 2, 1))         # one interim
builder_wt_2int <- make_builder_wt(c(1 / 3, 2 / 3, 1))  # two interims


# =============================================================================
# O'Brien-Fleming, whole shape scaled by one parameter p
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
#  Haybittle-Peto, whole shape scaled by one parameter p
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


