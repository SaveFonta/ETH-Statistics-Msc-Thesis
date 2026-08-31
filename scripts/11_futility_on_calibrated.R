# =============================================================================
# STEP 3 - Futility bolted onto the already calibrated designs
# =============================================================================
# The designs of 07_calibrate_designs.R are calibrated to a common avgT1E and
# a common avgPower WITHOUT any futility rule. Here we take those finished
# designs, leave n and the efficacy threshold p exactly as calibrated, and add
# the futility rule of 00_shared_setup.R (fut.crit) at every interim look,
# never at the final one. Nothing is recalibrated.
#
# This is deliberately not a like for like comparison of error rates: adding a
# futility rule can only remove trials, so it lowers both avgT1E and avgPower
# relative to the calibrated targets. Those shifts are reported next to the
# EUII rather than calibrated away, since the question here is what a sponsor
# actually buys by bolting futility onto a design that was already signed off.
#
# Only the SC family is evaluated. The DC designs of 07 sit at sample sizes far
# outside the range of a proof of concept trial and at an avgT1E of essentially
# zero, so their EUII is not interpretable (see the discussion in the
# Applications chapter).
#
# Run from the Thesis/ root, AFTER 07_calibrate_designs.R:
#   Rscript scripts/11_futility_on_calibrated.R
# =============================================================================

library(dplyr)
source("scripts/00_shared_setup.R")          # p_MAP, p_vague, fut.crit, ph1_ref
source("scripts/00_functions.R")

out_dir <- "Output/11_futility_on_calibrated"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
n_sim_eval   <- 1e6
delta_values <- c(0, seq(10, 100, by = 10))   # delta = 0 required by compute_euii
prior_H1     <- c(0.01, 0.1, 0.5)

cal <- readRDS("Output/07_calibrate_designs/designs.RDS")
sc_names <- grep("SC$", names(cal$designs_cal), value = TRUE)
cat("Calibrated SC designs found:", paste(sc_names, collapse = ", "), "\n")
cat("Targets: avgT1E =", cal$target_t1e, ", avgPower =", cal$target_power,
    "at delta =", cal$delta_MCID, "\n\n")


# ---------------------------------------------------------------------------
# Build the with futility twin of each calibrated design
#
# d$decisions is a list of length K, one entry per look, each
# list(success = <criterion>, futility = NULL). The efficacy criterion carries
# the calibrated threshold and is left untouched; only the futility slot is
# filled, and only at looks 1..K-1. A design with K = 1 (the fixed design) has
# no interim, so it has no futility twin.
# ---------------------------------------------------------------------------
add_futility <- function(decisions) {
  K <- length(decisions)
  if (K < 2) return(NULL)
  lapply(seq_len(K), function(k) {
    list(success  = decisions[[k]]$success,
         futility = if (k < K) fut.crit else NULL)
  })
}

variants <- list()
for (nm in sc_names) {
  d     <- cal$designs_cal[[nm]]
  looks <- sub(" \\| SC$", "", nm)
  variants[[paste0(nm, " | no futility")]] <-
    list(n1_seq = d$n1_seq, n2_seq = d$n2_seq, decisions = d$decisions,
         Base = nm, Looks = looks, Futility = "no futility", K = length(d$decisions))
  fut <- add_futility(d$decisions)
  if (!is.null(fut)) {
    variants[[paste0(nm, " | futility")]] <-
      list(n1_seq = d$n1_seq, n2_seq = d$n2_seq, decisions = fut,
           Base = nm, Looks = looks, Futility = "futility", K = length(d$decisions))
  }
}
cat("Evaluating", length(variants), "design variants at n_sim =", n_sim_eval, "\n\n")


# ---------------------------------------------------------------------------
# Simulate every variant on the delta grid
# ---------------------------------------------------------------------------
res_fut <- lapply(names(variants), function(nm) {
  v <- variants[[nm]]
  avgoc2_seq_mc.normMix(
    prior_1        = p_vague,
    prior_2        = p_MAP,
    n1_seq         = v$n1_seq,
    n2_seq         = v$n2_seq,
    decisions_list = v$decisions,
    delta          = delta_values,
    design_prior_c = p_MAP,
    n_sim          = n_sim_eval,
    seed           = cal$mc_seed
  )
})
names(res_fut) <- names(variants)


# ---------------------------------------------------------------------------
# EUII across delta, with the prior_H1 band
# ---------------------------------------------------------------------------
df_euii <- lapply(names(variants), function(nm) {
  v <- variants[[nm]]
  e <- compute_euii(res_fut[[nm]], prior_H1 = prior_H1)
  bind_rows(lapply(names(e), function(g) {
    e[[g]] |> mutate(prior_H1 = as.numeric(g),
                     Design   = nm,
                     Base     = v$Base,
                     Futility = v$Futility,
                     Looks    = v$Looks)
  }))
}) |> bind_rows()

band <- value_band(df_euii, "EUII", c("Design", "Base", "Futility", "Looks", "Delta"))


# ---------------------------------------------------------------------------
# Anchor table: what futility actually costs and buys at the calibration point
# ---------------------------------------------------------------------------
anchor <- lapply(names(variants), function(nm) {
  v  <- variants[[nm]]
  r  <- res_fut[[nm]]
  ov <- r[[paste0("delta.", cal$delta_MCID)]]$Overall
  nl <- r[["delta.0"]]$Overall
  eu <- band |> filter(Design == nm, Delta == cal$delta_MCID) |> pull(mid)
  data.frame(
    Design    = nm,
    Looks     = v$Looks,
    Futility  = v$Futility,
    n_max     = max(v$n1_seq) + max(v$n2_seq),
    avgT1E    = nl[["Power"]],
    avgPower  = ov[["Power"]],
    EN_null   = nl[["EN_t"]] + nl[["EN_c"]],
    EN_MCID   = ov[["EN_t"]] + ov[["EN_c"]],
    EUII      = eu
  )
}) |> bind_rows()

cat("Futility bolted onto the calibrated SC designs (nothing recalibrated).\n")
cat("Targets were avgT1E =", cal$target_t1e, "and avgPower =", cal$target_power,
    "at delta =", cal$delta_MCID, ".\n")
cat("EUII quoted at prior_H1 =", ph1_ref, ", delta =", cal$delta_MCID, ".\n\n")
print(anchor |> mutate(across(c(avgT1E, avgPower, EUII), \(x) round(x, 4)),
                       across(c(EN_null, EN_MCID), \(x) round(x, 1))),
      row.names = FALSE)

saveRDS(list(res_fut = res_fut, df_euii = df_euii, band = band, anchor = anchor,
             n_sim_eval = n_sim_eval, prior_H1 = prior_H1,
             target_t1e = cal$target_t1e, target_power = cal$target_power,
             delta_MCID = cal$delta_MCID, mc_seed = cal$mc_seed),
        file = file.path(out_dir, "futility.RDS"))
cat("\nResults saved to", file.path(out_dir, "futility.RDS"), "\n")
