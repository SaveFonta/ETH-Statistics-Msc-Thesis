# =============================================================================
# Sequential Design with TWO interim analyses — MC evaluation of AVG OC and EUII
#
# Same analysis as 2.3_SEQUENTIAL.MC.R, but the trial has three looks, timed as in
# the original protocol (Hueber et al. 2012): total enrollment of 17, 41 and 60
# patients at the two interims and the final analysis. With the 2:1 allocation
# this gives cumulative arm sizes that are only approximately 2:1 at the interims
# (11:6 and 27:14).
# Futility can stop the trial at both interims, efficacy at every look.
# Results saved to Output/SEQ2_MC_avg.RDS and Output/SEQ2_MC_euii.RDS.
# =============================================================================

library(RBesT)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(patchwork)   # plot_layout(guides = "collect") for a single shared legend
library(parallel)
source("scripts/prior_on_control/00.functions.MC.R")


# ---------------------------------------------------------------------------
# Design parameters
# ---------------------------------------------------------------------------
sigma      <- 88
delta_MCID <- 60
DV         <- 50   # decision value of the dual criterion (see dual.crit below)

n1_seq <- c(11, 27, 40)   # treatment arm cumulative sizes (three looks)
n2_seq <- c(6, 14, 20)    # control arm cumulative sizes, stage totals 17, 41, 60

delta_values <- c(0, seq(10, 100, by = 10))   
mc_seed      <- 123                           # shared by both sweeps (paired comparison)

# Prior probability of H1 in the EUII. It does not enter the simulation, only the
# processing with compute_euii, so the whole grid comes from a single
# sweep. 
prior_H1 <- c(0.01, 0.1, 0.5)


# ---------------------------------------------------------------------------
# Analysis priors
# ---------------------------------------------------------------------------
p_MAP <- mixnorm(
  c(0.4848, -52.457, 21.154), c(0.4598, -47.465, 7.843), c(0.0554, -50.355, 48.164),
  sigma = sigma, param = "ms"
)
p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")

analysis_priors <- list(
  "MAP"            = p_MAP,
  "Robust (w=0.2)" = robustify(p_MAP, 0.2, mean = -50),
  "Robust (w=0.4)" = robustify(p_MAP, 0.4, mean = -50),
  "Robust (w=0.6)" = robustify(p_MAP, 0.6, mean = -50),
  "Robust (w=0.8)" = robustify(p_MAP, 0.8, mean = -50),
  "Vague"          = p_vague
)
prior_names <- names(analysis_priors)

omega_map <- c(
  "MAP"            = 0.0,
  "Robust (w=0.2)" = 0.2,
  "Robust (w=0.4)" = 0.4,
  "Robust (w=0.6)" = 0.6,
  "Robust (w=0.8)" = 0.8,
  "Vague"          = 1.0
)


# ---------------------------------------------------------------------------
# Design priors
# ---------------------------------------------------------------------------
design_priors <- list(
  "Dirac (-50)"        = mixnorm(c(1, -50, 1e-16), sigma = sigma, param = "ms"),
  "MAP"                = p_MAP,
  "Skeptical (-90)"    = mixnorm(c(1, -90, 17.6), sigma = sigma, param = "ms"),
  "Misspecified (-10)" = mixnorm(c(1, -10, 25),   sigma = sigma, param = "ms")
)
dprior_names <- names(design_priors)


# ---------------------------------------------------------------------------
# Decision criteria
# ---------------------------------------------------------------------------
# qc is on the scale of theta_T - theta_C = -delta, so delta thresholds flip sign.
#
# Efficacy (both stages), dual criterion:
#   Pr(delta > 0 | data) > 0.95  AND  Pr(delta > 50 | data) > 0.50
dual.crit <- decision2S(pc = c(0.95, 0.50), qc = c(0, -50), lower.tail = TRUE)

sign.crit <- decision2S(pc = 0.95, qc = 0,  lower.tail = TRUE)

# Futility (both interims): Pr(delta < 40 | data) >= 0.90
#   Pr(delta < 40) = Pr(theta_T - theta_C > -40), and decision2S fires on a ">="
#   condition, hence the upper tail at qc = -40 with pc = 0.90.
fut.crit <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)

decisions_fut <- list(
  list(success = dual.crit, futility = fut.crit),   # interim 1
  list(success = dual.crit, futility = fut.crit),   # interim 2
  list(success = dual.crit, futility = NULL)        # final look
)
decisions_no_fut <- list(
  list(success = dual.crit, futility = NULL),
  list(success = dual.crit, futility = NULL),
  list(success = dual.crit, futility = NULL)
)

decisions_fut.single <- list(
  list(success = sign.crit, futility = fut.crit),   # interim 1
  list(success = sign.crit, futility = fut.crit),   # interim 2
  list(success = sign.crit, futility = NULL)        # final look
)
decisions_no_fut.single <- list(
  list(success = sign.crit, futility = NULL),
  list(success = sign.crit, futility = NULL),
  list(success = sign.crit, futility = NULL)
)


# ---------------------------------------------------------------------------
# Shared theme
# ---------------------------------------------------------------------------
my_theme <- theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())


# ---------------------------------------------------------------------------
# One MC run over every (analysis prior x design prior) combination
# ---------------------------------------------------------------------------
combos <- expand.grid(Analysis_Prior = prior_names,
                      Design_Prior   = dprior_names,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)



run_sweep <- function(decisions_list,
                      cores = 6, n_sim_avg = 1e6) {

  if (.Platform$OS.type != "unix") cores <- 1L   # mclapply cannot fork on Windows

  cat("\n", nrow(combos), " prior combinations, n_sim = ", n_sim_avg,
      ", cores = ", cores, "\n", sep = "")

  res <- parallel::mclapply(seq_len(nrow(combos)), function(i) {
    avgoc2_seq_mc.normMix(
      prior_1        = p_vague,
      prior_2        = analysis_priors[[combos$Analysis_Prior[i]]],
      n1_seq         = n1_seq,
      n2_seq         = n2_seq,
      decisions_list = decisions_list,
      delta          = delta_values,
      design_prior_c = design_priors[[combos$Design_Prior[i]]],
      n_sim          = n_sim_avg,
      seed           = mc_seed
    )
  }, mc.cores = cores)

  # mclapply reports failures as try-error elements instead of stopping
  bad <- vapply(res, inherits, logical(1), "try-error")
  if (any(bad)) {
    cat("MC failed for the following prior combinations:\n")
    cat(paste0("   ", combos$Analysis_Prior[bad], " | ", combos$Design_Prior[bad], "\n"), sep = "")
    stop("MC failed for ", sum(bad), " of ", nrow(combos), " prior combinations.")
  }
  res
}


summarise_total <- function (sweeped) {

list_res <- lapply(seq_along(sweeped), function(i) {
  # first extract each of the 24 analysis x design
  res <- sweeped[[i]]
  #  merge all the deltas evaluated
lapply(names(res), function(d){
    ov  <- res[[d]]$Overall
    Pow <- ov[["Power"]] #extract Power for each delta value
    data.frame(
      delta = sub("delta\\.", "", d)  |> as.numeric(),
      Power = Pow,
      # expected total sample size. 
      EN    = ov[["EN_t"]] + ov[["EN_c"]]
    )
  })  |> bind_rows()
})

lapply(seq_len(nrow(combos)), function(combo){
  list_res[[combo]]  |>  
      mutate(
        Analysis_Prior = combos$Analysis_Prior[[combo]], 
        Design_Prior = combos$Design_Prior[[combo]]
      )
})  |>  bind_rows()

}


summarise_per.stage <- function (sweeped) {
list_res <- lapply(seq_along(sweeped), function(i) {
  # extract specific combination
  res <- sweeped[[i]]

  lapply(names(res), function(d){ 
    res[[d]]$Per_Stage  |> 
              select(Stage, P_Succ,  P_Fut, Cum_P_Succ, Cum_P_Fut)  |> 
              mutate(delta = sub("delta\\.", "", d))
  })  |>  bind_rows()
})


lapply (seq_len(nrow(combos)), function(combo){
  list_res[[combo]]  |>  mutate(
        Analysis_Prior = combos$Analysis_Prior[[combo]], 
        Design_Prior = combos$Design_Prior[[combo]]
  )
})  |>  bind_rows()
}



# EUII for every prior_H1, 
euii_table <- function(sweep) {
  bind_rows(lapply(seq_along(sweep), function(i) {
    e <- compute_euii(sweep[[i]], prior_H1 = prior_H1)
    bind_rows(lapply(names(e), function(g) {
      e[[g]] |>
        mutate(prior_H1       = as.numeric(g),
               Analysis_Prior = combos$Analysis_Prior[i],
               Design_Prior   = combos$Design_Prior[i],
               omega          = omega_map[[combos$Analysis_Prior[i]]])
    }))
  })) |>
    mutate(Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
           Design_Prior   = factor(Design_Prior,   levels = dprior_names))
}




# =============================================================================
# Run one MC sweep per decision criterion
# =============================================================================
# Four designs are compared. All use the SAME seed and the SAME n_sim, so theta_C
# and the stage 1 data are common to all four and the comparison is paired.
n_sim_avg <- 1e6

designs <- list(
  "SC"            = decisions_no_fut.single,
  "SC + Futility" = decisions_fut.single,
  "DC"            = decisions_no_fut,
  "DC + Futility" = decisions_fut
)
design_labels <- names(designs)

sweeps <- lapply(design_labels, function(nm) {
  cat("\n=== ", nm, " ===")
  run_sweep(designs[[nm]], n_sim_avg = n_sim_avg)
})
names(sweeps) <- design_labels


# ---------------------------------------------------------------------------
# Stack the four sweeps, tagging each with its Criterion
# ---------------------------------------------------------------------------
# Total OC. create one column for avgt1e
df_oc <- bind_rows(lapply(design_labels, function(nm) {
  summarise_total(sweeps[[nm]]) |> mutate(Criterion = nm)
})) |>
  group_by(Criterion, Analysis_Prior, Design_Prior) |>
  mutate(avgT1E = Power[delta == 0]) |>
  ungroup() 



# Per stage OC
df_stage <- bind_rows(lapply(design_labels, function(nm) {
  summarise_per.stage(sweeps[[nm]]) |> mutate(Criterion = nm)
})) |>
  mutate(delta = as.numeric(delta))


# Trasform to factors
df_oc <- df_oc  |>  mutate(
  Criterion = factor(Criterion, levels = design_labels),
    Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
    Design_Prior   = factor(Design_Prior,   levels = dprior_names)
  )


df_stage <- df_stage  |>  mutate(
  Criterion = factor(Criterion, levels = design_labels),
    Analysis_Prior = factor(Analysis_Prior, levels = prior_names),
    Design_Prior   = factor(Design_Prior,   levels = dprior_names)
  )




# =============================================================================
# PART 1 — Average T1E and average Power, one panel per criterion
# =============================================================================
# function to create the panel
oc_layout <- function(p) {
  p +
    geom_point(size = 2.2) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~Criterion, nrow = 1) +
    scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"), guide = "none") +
    scale_color_viridis_d(option = "turbo", direction = -1) +
    my_theme +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
}

t1e <- df_oc  |>  filter (delta == 0)
p_t1e <- oc_layout(
  ggplot(t1e,
         aes(x = Design_Prior, y = Power,
             color = Analysis_Prior, group = Analysis_Prior,
             linetype = Analysis_Prior == "Vague"))) +
  geom_hline(yintercept = 0.05, linetype = "dashed", color = "red", linewidth = 0.5) +
  labs(x = "Design Prior", y = "Average T1E")

print(p_t1e)

p_pow <- oc_layout(
  ggplot(filter(df_oc, delta == delta_MCID),
         aes(x = Design_Prior, y = Power,
             color = Analysis_Prior, group = Analysis_Prior,
             linetype = Analysis_Prior == "Vague"))) +
  labs(x = "Design Prior",
       y = bquote("Average Power"))

print(p_pow)

# Expected sample size. Futility leaves the average T1E and the average Power
# basically unchanged, because it only stops trials that were going to fail
# anyway. --> from here we don't see the advantage of futility analyis. Luckily WE 
# HAVE DEVELOPED THE EUII




# =============================================================================
# PART 2 — Per stage behaviour, one panel per criterion
# =============================================================================
# What the interims actually do: how often the trial stops early for efficacy when
# the treatment works, and how often it stops early for futility under the null.
# With two interims there is one row per interim and per metric, four rows in total.
df_stage1 <- bind_rows(
  df_stage |> filter(Stage <= 2, delta == delta_MCID) |>
    transmute(Criterion, Analysis_Prior, Design_Prior, Prob = P_Succ,
              Metric = paste0("Efficacy stop, interim ", Stage,
                              " (delta = ", delta_MCID, ")")),
  df_stage |> filter(Stage <= 2, delta == 0) |>
    transmute(Criterion, Analysis_Prior, Design_Prior, Prob = P_Fut,
              Metric = paste0("Futility stop, interim ", Stage, " (delta = 0)"))
)

p_stage <- ggplot(df_stage1,
                  aes(x = Design_Prior, y = Prob,
                      color = Analysis_Prior, group = Analysis_Prior,
                      linetype = Analysis_Prior == "Vague")) +
  geom_point(size = 2.2) +
  geom_line(linewidth = 0.9) +
  facet_grid(Metric ~ Criterion) +
  scale_linetype_manual(values = c("TRUE" = "dashed", "FALSE" = "solid"), guide = "none") +
  scale_color_viridis_d(option = "turbo", direction = -1) +
  labs(x = "Design Prior", y = "Probability of stopping at the interim") +
  my_theme +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p_stage)

saveRDS(list(df_oc = df_oc, df_stage = df_stage, sweeps = sweeps),
        file = "Output/SEQ2_MC_avg.RDS")
cat("\nMC average and per stage results saved to Output/SEQ2_MC_avg.RDS\n")


# =============================================================================
# PART 3 — Sequential EUII, 
# =============================================================================
# hERE, WE ONLY USE THE correctly specified analysis prior (so MAP for both analysis and design prior)
# Left: EUII of the four sequential designs, with the two fixed designs dashed.
# Right: ratio of each sequential design to ITS OWN fixed counterpart, so SC and
# SC + Futility are divided by the fixed SC, and DC and DC + Futility by the fixed DC.


df_euii <- bind_rows(lapply(design_labels, function(nm) {
  euii_table(sweeps[[nm]]) |> mutate(Criterion = nm)
})) |>
  filter(Analysis_Prior == "MAP", Design_Prior == "MAP") |>
  mutate(Criterion = factor(Criterion, levels = design_labels),
         Base      = if_else(grepl("^SC", as.character(Criterion)), "SC", "DC"))

## Fixed baselines
read_fixed <- function(path, base_label) {
  readRDS(path)$df_euii |>
    filter(Analysis_Prior == "MAP", Design_Prior == "MAP") |>
    transmute(Delta = delta, EUII_fixed = EUII, Base = base_label)
}
fixed_base <- bind_rows(
  read_fixed("Output/FIXED_euii.RDS", "SC"),   
  read_fixed("Output/DUAL_euii.RDS",  "DC")    
)

df_euii <- df_euii |>
  left_join(fixed_base, by = c("Delta", "Base")) |>
  mutate(EUII_ratio = EUII / EUII_fixed)



# --- colour encodes the DESIGN --------------------------------------------------
crit_levels <- c(design_labels, "Fixed SC", "Fixed DC")

crit_cols <- c("SC" = "#0072B2", "SC + Futility" = "#56B4E9",
               "DC" = "#D55E00", "DC + Futility" = "#E69F00",
               "Fixed SC" = "#0072B2", "Fixed DC" = "#D55E00")





# prior_H1 shifts the EUII by about one percent and never changes the ranking of
# the designs, so plotting one line per value only makes the plot more compilcated Instead the line is
# drawn at prior_H1 = 0.5 and use the ribbon for the others
ph1_ref <- 0.5

euii_band <- function(df, value_col) {
  df |>
    mutate(Criterion = factor(as.character(Criterion), levels = crit_levels)) |>
    group_by(Criterion, Delta) |>
    summarise(lo  = min(.data[[value_col]]),
              hi  = max(.data[[value_col]]),
              mid = .data[[value_col]][prior_H1 == ph1_ref],
              .groups = "drop")
}

df_band_abs   <- euii_band(df_euii, "EUII")
df_band_ratio <- euii_band(df_euii, "EUII_ratio")

# The two fixed designs carry no prior_H1: in a fixed design N is constant. They are single reference lines, with no band.
df_fixed_lines <- fixed_base |>
  transmute(Delta, EUII = EUII_fixed,
            Criterion = factor(paste("Fixed", Base), levels = crit_levels))

# A fixed design shares the colour of its sequential counterpart, so colour alone cannot
# tell "SC" from "Fixed SC". The linetype is therefore mapped to Criterion as well. Since
# colour and linetype then have the same name and the same breaks, ggplot merges them into
# a SINGLE legend whose keys show both the real colour and the real linetype.
crit_ltys <- c("SC" = "solid", "SC + Futility" = "solid",
               "DC" = "solid", "DC + Futility" = "solid",
               "Fixed SC" = "twodash", "Fixed DC" = "twodash")

# --- top panel: absolute EUII, with the two fixed designs as references ---------
# The panels are stacked and share the delta axis, so the top one drops its x axis.
p_euii_abs <- ggplot() +
  geom_ribbon(data = df_band_abs,
              aes(x = Delta, ymin = lo, ymax = hi, fill = Criterion), alpha = 0.30) +
  geom_line(data = df_band_abs,
            aes(x = Delta, y = mid, color = Criterion, linetype = Criterion),
            linewidth = 1.2) +
  geom_line(data = df_fixed_lines,
            aes(x = Delta, y = EUII, color = Criterion, linetype = Criterion),
            linewidth = 1.0) +
  geom_hline(yintercept = 1, linetype = "dotted", color = "black", linewidth = 0.4) +
  scale_color_manual(values = crit_cols, drop = FALSE) +
  scale_linetype_manual(values = crit_ltys, drop = FALSE) +
  scale_fill_manual(values = crit_cols, guide = "none") +
  guides(color = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  coord_cartesian(xlim = c(DV, max(delta_values))) +
  labs(x = NULL, y = "EUII") +
  my_theme +
  theme(axis.text.x = element_blank())

# --- bottom panel: ratio to the corresponding fixed design ---------------------

p_euii_ratio <- ggplot(df_band_ratio) +
  geom_ribbon(aes(x = Delta, ymin = lo, ymax = hi, fill = Criterion), alpha = 0.30) +
  geom_line(aes(x = Delta, y = mid, color = Criterion, linetype = Criterion),
            linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = crit_cols) +
  scale_linetype_manual(values = crit_ltys) +
  scale_fill_manual(values = crit_cols, guide = "none") +
  guides(color = "none", linetype = "none") +
  coord_cartesian(xlim = c(DV, max(delta_values))) +
  labs(x = expression(delta),
       y = "EUII ratio") +
  my_theme

# patchwork stacks the panels and collects the guides, so they share one legend
print((p_euii_abs / p_euii_ratio) +
        patchwork::plot_layout(guides = "collect") &
        theme(legend.position = "bottom"))

saveRDS(list(df_euii = df_euii, prior_H1 = prior_H1),
        file = "Output/SEQ2_MC_euii.RDS")
cat("MC EUII results saved to Output/SEQ2_MC_euii.RDS\n")




###########################
# -------- Effective sample size, 1 / E[1/N | outcome]
#
# As in the EUII panels, the line is drawn at prior_H1 = 0.5 and there is ribbon for the others

library(tidyr)

n_min <- n1_seq[1] + n2_seq[1]                              #  first interim total, 20
n_mid <- n1_seq[2] + n2_seq[2]                              #  second interim total, 40
n_max <- n1_seq[length(n1_seq)] + n2_seq[length(n2_seq)]    #  final total, 60


sample_size <- df_euii |>
  filter(grepl("Futility", Criterion)) |>
  select(Delta, prior_H1, Criterion, E_invN_sig, E_invN_nonsig) |>
  pivot_longer(cols = c(E_invN_sig, E_invN_nonsig),
               names_to = "E_inv_type", values_to = "E_inv_val") |>
  mutate(
    N_eff     = 1 / E_inv_val,     
    Outcome   = factor(if_else(E_inv_type == "E_invN_sig", "Success", "Non success"),
                       levels = c("Success", "Non success")),
    Criterion = factor(as.character(Criterion), levels = design_labels)
  )

# band to generate the ribbon
band_N <- sample_size |>
  group_by(Criterion, Outcome, Delta) |>
  summarise(lo  = min(N_eff),
            hi  = max(N_eff),
            mid = N_eff[prior_H1 == ph1_ref],
            .groups = "drop")

p_N_eff <- ggplot(band_N) +
  geom_ribbon(aes(x = Delta, ymin = lo, ymax = hi, fill = Criterion), alpha = 0.30) +
  geom_line(aes(x = Delta, y = mid, color = Criterion), linewidth = 1.2) +
  geom_hline(yintercept = n_max, linetype = "dashed", color = "black",  linewidth = 0.4) +
  geom_hline(yintercept = n_mid, linetype = "dotted", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = n_min, linetype = "dotted", color = "grey40", linewidth = 0.4) +
  annotate("text", x = DV + 2, y = n_max - 1.3, hjust = 0, size = 3, color = "grey30",
           label = paste0("without futility (= ", n_max, ")")) +
  annotate("text", x = DV + 2, y = n_mid + 1.1, hjust = 0, size = 3, color = "grey50",
           label = paste0("second interim (= ", n_mid, ")")) +
  annotate("text", x = DV + 2, y = n_min + 1.1, hjust = 0, size = 3, color = "grey50",
           label = paste0("first interim (= ", n_min, ")")) +
  facet_wrap(~Outcome) +
  scale_color_manual(values = crit_cols) +
  scale_fill_manual(values = crit_cols, guide = "none") +
  # same delta window as the EUII panels, plus the sample size bounds
  coord_cartesian(xlim = c(DV, max(delta_values)), ylim = c(n_min, n_max)) +
  labs(x = expression(delta),
       y = expression("Effective sample size  " * (E * "[1/N | outcome]")^-1)) +
  my_theme

print(p_N_eff)


###########################
# -------- Two interims against one: EUII ratio
#
# The natural question of this script: does the second interim add information per
# patient beyond the first? The one interim results of 2.3_SEQUENTIAL.MC.R are read
# from Output/SEQ_MC_euii.RDS and the ratio EUII(two interims) / EUII(one interim)
# is computed per design and prior_H1. Both runs share the same seed, so theta_C is
# common and the comparison is paired.

df_euii_1int <- readRDS("Output/SEQ_MC_euii.RDS")$df_euii |>
  select(Delta, prior_H1, Criterion, EUII_1int = EUII)

df_ratio_12 <- df_euii |>
  select(Delta, prior_H1, Criterion, EUII) |>
  mutate(Criterion = as.character(Criterion)) |>
  left_join(df_euii_1int |> mutate(Criterion = as.character(Criterion)),
            by = c("Delta", "prior_H1", "Criterion")) |>
  mutate(ratio_12 = EUII / EUII_1int)

band_12 <- df_ratio_12 |>
  mutate(Criterion = factor(Criterion, levels = design_labels)) |>
  group_by(Criterion, Delta) |>
  summarise(lo  = min(ratio_12),
            hi  = max(ratio_12),
            mid = ratio_12[prior_H1 == ph1_ref],
            .groups = "drop")

p_ratio_12 <- ggplot(band_12) +
  geom_ribbon(aes(x = Delta, ymin = lo, ymax = hi, fill = Criterion), alpha = 0.30) +
  geom_line(aes(x = Delta, y = mid, color = Criterion), linewidth = 1.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 0.5) +
  scale_color_manual(values = crit_cols) +
  scale_fill_manual(values = crit_cols, guide = "none") +
  coord_cartesian(xlim = c(DV, max(delta_values))) +
  labs(x = expression(delta),
       y = "EUII (two interims) / EUII (one interim)") +
  my_theme

print(p_ratio_12)
