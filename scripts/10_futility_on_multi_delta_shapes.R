# =============================================================================
# Multi delta calibration - futility added on the calibrated designs 
# =============================================================================
# Add futility onto the boundary shapes calibrated in
# 08_calibrate_shapes_multi_delta.R. 
#
# n and the efficacy threshold are left exactly as calibrated in script 08;
# only the futility rule is added at the two interims (never the final
# look). Not a comparison at equal error rates: futility can only remove
# trials, so it lowers avgT1E and avgPower slightly below the calibrated
# targets, and those shifts are reported next to the EUII rather than tuned
# away.
#
# Requires 08_calibrate_shapes_multi_delta.R to have already been run.
# =============================================================================

library(dplyr)
source("scripts/00_shared_setup.R")          # p_MAP, p_vague, fut.crit
source("scripts/00_functions.R")             # avgoc2_seq_mc.normMix, compute_euii

out_dir <- "Output/10_futility_on_multi_delta_shapes"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)


# ---------------------------------------------------------------------------
# Load the calibrated shapes from script 08 (no futility, one per
# shape x delta target, K = 3 throughout)
# ---------------------------------------------------------------------------
cal <- readRDS("Output/08_calibrate_shapes_multi_delta/results.RDS")
combos      <- cal$combos
designs_cal <- cal$designs_cal
n_sim_eval  <- cal$n_sim_check
mc_seed     <- cal$mc_seed
prior_H1    <- cal$prior_H1
shape_names <- cal$shape_names

cat("Loaded", nrow(combos), "calibrated shapes from script 08 (delta_targets:",
    paste(cal$delta_targets, collapse = ", "), ")\n\n")


# ---------------------------------------------------------------------------
# Futility twin: the futility rule is added at both interims
# (K = 3 throughout here), never at the final look. 
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
for (i in seq_len(nrow(combos))) {
  d  <- designs_cal[[i]]
  nm <- combos$Label[i]
  variants[[paste0(nm, " | no futility")]] <-
    list(n1_seq = d$n1_seq, n2_seq = d$n2_seq, decisions = d$decisions,
         Shape = combos$Shape[i], DeltaTarget = combos$DeltaTarget[i],
         Futility = "no futility")
  fut <- add_futility(d$decisions)
  variants[[paste0(nm, " | futility")]] <-
    list(n1_seq = d$n1_seq, n2_seq = d$n2_seq, decisions = fut,
         Shape = combos$Shape[i], DeltaTarget = combos$DeltaTarget[i],
         Futility = "futility")
}
cat("Evaluating", length(variants), "design variants at n_sim =", n_sim_eval, "\n\n")

cores <- min(32L, length(variants))
if (.Platform$OS.type != "unix") cores <- 1L


# ---------------------------------------------------------------------------
# Simulate every variant at its own delta_target
# ---------------------------------------------------------------------------
results <- parallel::mclapply(names(variants), function(nm) {
  v  <- variants[[nm]]
  dt <- v$DeltaTarget

  r <- avgoc2_seq_mc.normMix(
    prior_1        = p_vague,
    prior_2        = p_MAP,
    n1_seq         = v$n1_seq,
    n2_seq         = v$n2_seq,
    decisions_list = v$decisions,
    delta          = c(0, dt),
    design_prior_c = p_MAP,
    n_sim          = n_sim_eval,
    seed           = mc_seed + 1
  )

  ov <- r[[paste0("delta.", dt)]]$Overall
  nl <- r[["delta.0"]]$Overall
  # EUII at every prior_H1 in the grid, so the plot can show the same
  # Pr(H1)-sensitivity band as the other EUII figures: EUII is the point
  # value at Pr(H1) = 0.5 (used throughout the prose and tables), EUII_lo/hi
  # the min/max across the full prior_H1 grid (used only for the ribbon).
  eu_all  <- compute_euii(r, prior_H1 = prior_H1)
  eu_vals <- vapply(eu_all, function(e) e$EUII[e$Delta == dt], numeric(1))
  # Effective sample size on each branch, across the full prior_H1 grid, for
  # the third panel of the EUII figures: (E[1/N | success])^-1 and
  # (E[1/N | failure])^-1. EN_success/EN_failure are the point values at
  # Pr(H1) = 0.5 (used throughout the prose and tables), the _lo/_hi the
  # min/max across the grid (used only for the ribbon).
  en_success_vals <- vapply(eu_all, function(e) 1 / e$E_invN_sig[e$Delta == dt], numeric(1))
  en_failure_vals <- vapply(eu_all, function(e) 1 / e$E_invN_nonsig[e$Delta == dt], numeric(1))

  list(n_max    = max(v$n1_seq) + max(v$n2_seq),
       avgT1E   = nl[["Power"]],
       avgPower = ov[["Power"]],
       EN_null  = nl[["EN_t"]] + nl[["EN_c"]],
       EN_MCID  = ov[["EN_t"]] + ov[["EN_c"]],
       EN_success    = en_success_vals[["0.5"]],
       EN_success_lo = min(en_success_vals),
       EN_success_hi = max(en_success_vals),
       EN_failure    = en_failure_vals[["0.5"]],
       EN_failure_lo = min(en_failure_vals),
       EN_failure_hi = max(en_failure_vals),
       EUII     = eu_vals[["0.5"]],
       EUII_lo  = min(eu_vals),
       EUII_hi  = max(eu_vals))
}, mc.cores = cores)
names(results) <- names(variants)

bad <- vapply(results, inherits, logical(1), "try-error")
if (any(bad)) {
  cat("Evaluation failed for:\n"); cat(paste0("   ", names(results)[bad], "\n"), sep = "")
  stop("evaluation failed for ", sum(bad), " of ", length(variants), " variants.")
}


# ---------------------------------------------------------------------------
# Assemble
# ---------------------------------------------------------------------------
anchor <- lapply(names(variants), function(nm) {
  v <- variants[[nm]]; r <- results[[nm]]
  data.frame(DeltaTarget = v$DeltaTarget, Shape = v$Shape, Futility = v$Futility,
             n_max = r$n_max, avgT1E = r$avgT1E, avgPower = r$avgPower,
             EN_null = r$EN_null, EN_MCID = r$EN_MCID,
             EN_success = r$EN_success, EN_success_lo = r$EN_success_lo, EN_success_hi = r$EN_success_hi,
             EN_failure = r$EN_failure, EN_failure_lo = r$EN_failure_lo, EN_failure_hi = r$EN_failure_hi,
             EUII = r$EUII, EUII_lo = r$EUII_lo, EUII_hi = r$EUII_hi)
}) |> bind_rows() |>
  mutate(Shape = factor(Shape, levels = shape_names)) |>
  arrange(DeltaTarget, Shape, Futility)

cat("Futility bolted onto the multi delta calibrated shapes (nothing recalibrated).\n")
cat("EUII quoted at Pr(H1) = 0.5, at each design's own delta target, K = 3 throughout.\n\n")
print(as.data.frame(anchor |>
        mutate(across(c(avgT1E, avgPower, EUII), \(x) round(x, 4)),
               across(c(EN_null, EN_MCID), \(x) round(x, 1)))),
      row.names = FALSE)


# ---------------------------------------------------------------------------
# The actual question: does futility change the ranking of shapes, and does
# it help the strict boundaries (Haybittle-Peto) relatively more than the
# permissive ones (Pocock)?
# ---------------------------------------------------------------------------
cat("\n=== Ranking check: does futility change which shape wins? ===\n")
for (dt in unique(anchor$DeltaTarget)) {
  sub <- anchor |> filter(DeltaTarget == dt)
  best_no_fut <- sub |> filter(Futility == "no futility") |> slice_max(EUII, n = 1)
  best_fut    <- sub |> filter(Futility == "futility")    |> slice_max(EUII, n = 1)
  cat(" delta =", format(dt, width = 2),
      ": best w/o futility =", format(as.character(best_no_fut$Shape), width = 15),
      "(", formatC(best_no_fut$EUII, format = "f", digits = 4), ")",
      "| best w/ futility =", format(as.character(best_fut$Shape), width = 15),
      "(", formatC(best_fut$EUII, format = "f", digits = 4), ")\n")
}

cat("\n=== Relative EUII gain from futility, by shape (average across targets) ===\n")
gain <- anchor |>
  select(DeltaTarget, Shape, Futility, EUII) |>
  tidyr::pivot_wider(names_from = Futility, values_from = EUII) |>
  mutate(rel_gain = (futility - `no futility`) / `no futility`)

gain_summary <- gain |> group_by(Shape) |> summarise(mean_rel_gain = mean(rel_gain), .groups = "drop") |>
  arrange(desc(mean_rel_gain))
print(as.data.frame(gain_summary |> mutate(mean_rel_gain = formatC(100 * mean_rel_gain, format = "f", digits = 2))),
      row.names = FALSE)
cat("(values are % EUII gain from adding futility, averaged over delta_target =",
    paste(range(anchor$DeltaTarget), collapse = ".."), ")\n")

# gain_summary is printed above as a sanity check but not saved: Applications.Rnw
# reads `gain` directly and recomputes this same summary itself, so keeping a
# second copy in the RDS would just be dead weight.
saveRDS(list(anchor = anchor, gain = gain,
             n_sim_eval = n_sim_eval, prior_H1 = prior_H1, mc_seed = mc_seed),
        file = file.path(out_dir, "results.RDS"))
cat("\nResults saved to", file.path(out_dir, "results.RDS"), "\n")
