# =============================================================================
# Sequential Design — Monte Carlo evaluation of AVERAGE OC and EUII
# =============================================================================

library(RBesT)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(parallel)
source("scripts/prior_on_control/00.functions.MC.R")


# ---------------------------------------------------------------------------
# Design parameters
# ---------------------------------------------------------------------------
sigma      <- 88
delta_MCID <- 60

n1_seq <- c(20, 40)   # treatment arm cumulative sizes
n2_seq <- c(10, 20)   # control arm cumulative sizes

delta_values <- c(0, seq(10, 100, by = 10))   
mc_seed      <- 123                           # shared by both sweeps (paired comparison)

# Prior probability of H1 in the EUII. It does not enter the simulation, only the
# post processing inside compute_euii, so the whole grid comes from a single
# sweep. The EUII is reported for each value.
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
sign.crit <- decision2S(pc = c(0.95, 0.50), qc = c(0, -50), lower.tail = TRUE)

# Futility (stage 1 only): Pr(delta < 40 | data) >= 0.90
#   Pr(delta < 40) = Pr(theta_T - theta_C > -40), and decision2S fires on a ">="
#   condition, hence the upper tail at qc = -40 with pc = 0.90.
fut.crit <- decision2S(pc = 0.90, qc = -40, lower.tail = FALSE)

decisions_fut <- list(
  list(success = sign.crit, futility = fut.crit),   # stage 1: efficacy + futility
  list(success = sign.crit, futility = NULL)        # stage 2: final, efficacy only
)
decisions_no_fut <- list(
  list(success = sign.crit, futility = NULL),       # stage 1: efficacy only
  list(success = sign.crit, futility = NULL)        # stage 2: final
)


# ---------------------------------------------------------------------------
# Shared theme
# ---------------------------------------------------------------------------
my_theme <- theme_minimal() +
  theme(legend.position = "bottom", legend.title = element_blank())


# ---------------------------------------------------------------------------
# Helpers: one MC run over every (analysis prior x design prior) combination
# ---------------------------------------------------------------------------
combos <- expand.grid(Analysis_Prior = prior_names,
                      Design_Prior   = dprior_names,
                      stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)

# The combinations are independent, so they run in parallel. mclapply forks, so
# the workers inherit the priors, the decisions and the sourced functions for

run_sweep <- function(decisions_list,
                      cores = 6, n_sim_avg = 10000) {

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
    stop("MC failed for ", sum(bad), " of ", nrow(combos), " prior combinations (listed above).")
  }
  res
}


summarise_total <- function (sweeped) {

list_res <- lapply(seq_along(sweeped), function(i) {
  # first extract each of the 24 analysis x design
  res <- sweeped[[i]]
  #  merge all the deltas evaluated
lapply(names(res), function(d){
    Pow <- res[[d]]$Overall[["Power"]] #extract Power for each delta value
    data.frame(
      delta = sub("delta\\.", "", d)  |> as.numeric(),
      Power = Pow
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
# Run both MC sweeps
# =============================================================================
mc_fut    <- run_sweep(decisions_fut)
mc_no_fut <- run_sweep(decisions_no_fut)


summarise_per.stage(mc_fut)



# =============================================================================
# PART 1 — Average OC
# =============================================================================
summarise_total(mc_fut)
summarise_total(mc_no_fut)

summarise_per.stage(mc_fut)


# Now we need to make plots about this:
# make the classic plots of avgT1E and avg Power thnx 

# The idea is to have two panel of 4. Each one for one criteria:
1) SC 
2) SC + futility 
3) DC 
4) DC + futility 


# Then the first panel is avgT1E for each of the 4 combinations
# The second is power at delta = 60 for each of the 4 combinations 

#Maybe add also a Power varying delta? Need to think about it. 

Assumiamo correctly specified prior e compute per ognuno dei 4:
avg Metrics nello stesso grafico a sx, solo con 4 colori diversi. Nel grafico a dx invece ci metti il rapporto con il loro corrispettivo Fixed. Quindi se è un SC gli metti 
SC, se è DC gli metti DC. 




# =============================================================================
# PART 2 — Sequential EUII, for the moment we inspect only the correctly specified one
# =============================================================================
df_euii_seq_mc <- euii_table(mc_fut)
df_euii_no_fut <- euii_table(mc_no_fut)


# So we will do a comparison of 
SC vs SC futility vs DC vs DC futility and then also add the Fixed SC and Dc
# Then the plot next to it is the EUII compare 
Assumiamo correctly specified prior e compute per ognuno dei 4:
avg Metrics nello stesso grafico a sx, solo con 4 colori diversi. Nel grafico a dx invece ci metti il rapporto con il loro corrispettivo Fixed. Quindi se è un SC gli metti 
SC, se è DC gli metti DC. 
