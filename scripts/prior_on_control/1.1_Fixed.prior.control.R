source("scripts/prior_on_control/00.functions.R")


# this part will use it for the Rnw
library(RBesT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(checkmate)  
library(assertthat) 

# Create data directory if it doesn't exist
dir.create("data", showWarnings = FALSE)

sigma <- 88


# define success criteria
sign.crit.95 <- decision2S(pc = 0.95, qc = 0, lower.tail = TRUE)
sign.crit.975 <- decision2S(pc = 0.975, qc = 0, lower.tail = TRUE)

dual.crit.95 <- decision2S(pc = c(0.95, 0.5), qc = c(0, -50), lower.tail = TRUE)
dual.crit.975 <- decision2S(pc = c(0.975, 0.5), qc = c(0, -50), lower.tail = TRUE)





# define possible priors







p_MAP  <- mixnorm(c(0.51, -51, 19.9), c(0.44, -46.8, 7.6), c(0.05, -54.1, 51.7),
                  sigma = sigma, param = "ms")


  # I use this in the evaluation of the classic T1E and don't want to change 
p_rob <- robustify(p_MAP, 0.2, mean = -50)
# for the evaluation of Assurance, i use the versions with varying weights


p_rob0.2  <- robustify(p_MAP, 0.2, mean = -50)

p_rob_0.5 <- robustify(p_MAP, 0.5, mean = -50)


p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")

prior.t <- p_vague # Ok I decided to make prior.t super vague

p_skep  <- mixnorm(c(1, -90, 25), sigma = sigma, param = "ms") #posterior distribution of a stand-alone analysis of the historical study with the “most extreme” placebo effect

Normal <- mixnorm(c(1, -49, 20), sigma = sigma, param = "mn")


n.act <- 40
n.pbo <- 20

# --------------------
# CREATE OC
# --------------------


oc_vague.sign.95 <- oc2S(prior1 = p_vague, prior2 = p_vague, n1 = n.act, n2 = n.pbo, decision = sign.crit.95)
oc_MAP.sign.95 <- oc2S(prior1 = p_vague, prior2 = p_MAP, n1 = n.act, n2 = n.pbo, decision = sign.crit.95)
oc_rob.sign.95 <- oc2S(prior1 = p_vague, prior2 = p_rob, n1 = n.act, n2 = n.pbo, decision = sign.crit.95)

oc_vague.sign.975 <- oc2S(prior1 = p_vague, prior2 = p_vague, n1 = n.act, n2 = n.pbo, decision = sign.crit.975)
oc_MAP.sign.975 <- oc2S(prior1 = p_vague, prior2 = p_MAP, n1 = n.act, n2 = n.pbo, decision = sign.crit.975)
oc_rob.sign.975 <- oc2S(prior1 = p_vague, prior2 = p_rob, n1 = n.act, n2 = n.pbo, decision = sign.crit.975)

oc_vague.dual.95 <- oc2S(prior1 = p_vague, prior2 = p_vague, n1 = n.act, n2 = n.pbo, decision = dual.crit.95)
oc_MAP.dual.95 <- oc2S(prior1 = p_vague, prior2 = p_MAP, n1 = n.act, n2 = n.pbo, decision = dual.crit.95)
oc_rob.dual.95 <- oc2S(prior1 = p_vague, prior2 = p_rob, n1 = n.act, n2 = n.pbo, decision = dual.crit.95)

oc_vague.dual.975 <- oc2S(prior1 = p_vague, prior2 = p_vague, n1 = n.act, n2 = n.pbo, decision = dual.crit.975)
oc_MAP.dual.975 <- oc2S(prior1 = p_vague, prior2 = p_MAP, n1 = n.act, n2 = n.pbo, decision = dual.crit.975)
oc_rob.dual.975 <- oc2S(prior1 = p_vague, prior2 = p_rob, n1 = n.act, n2 = n.pbo, decision = dual.crit.975)



# -------------------
# Evalute classic T1E
# -------------------
theta_c <- seq(-150, 50)
theta_a <- theta_c

# Significance 95
T1E_vague.sign.95 <- oc_vague.sign.95(theta_c, theta_a)
T1E_MAP.sign.95   <- oc_MAP.sign.95(theta_c, theta_a)
T1E_rob.sign.95   <- oc_rob.sign.95(theta_c, theta_a)

# Significance 97.5
T1E_vague.sign.975 <- oc_vague.sign.975(theta_c, theta_a)
T1E_MAP.sign.975   <- oc_MAP.sign.975(theta_c, theta_a)
T1E_rob.sign.975   <- oc_rob.sign.975(theta_c, theta_a)

# Dual 95
T1E_vague.dual.95 <- oc_vague.dual.95(theta_c, theta_a)
T1E_MAP.dual.95   <- oc_MAP.dual.95(theta_c, theta_a)
T1E_rob.dual.95   <- oc_rob.dual.95(theta_c, theta_a)

# Dual 97.5
T1E_vague.dual.975 <- oc_vague.dual.975(theta_c, theta_a)
T1E_MAP.dual.975   <- oc_MAP.dual.975(theta_c, theta_a)
T1E_rob.dual.975   <- oc_rob.dual.975(theta_c, theta_a)


# -------------------
# Create df
# -------------------
df <- data.frame(
  theta = theta_c,
  "Vague_Significance (95%)"   = T1E_vague.sign.95,
  "MAP_Significance (95%)"     = T1E_MAP.sign.95,
  "Robust_Significance (95%)"  = T1E_rob.sign.95,
  
  "Vague_Significance (97.5%)" = T1E_vague.sign.975,
  "MAP_Significance (97.5%)"   = T1E_MAP.sign.975,
  "Robust_Significance (97.5%)"= T1E_rob.sign.975,
  
  "Vague_Dual (95%)"       = T1E_vague.dual.95,
  "MAP_Dual (95%)"         = T1E_MAP.dual.95,
  "Robust_Dual (95%)"      = T1E_rob.dual.95,
  
  "Vague_Dual (97.5%)"     = T1E_vague.dual.975,
  "MAP_Dual (97.5%)"       = T1E_MAP.dual.975,
  "Robust_Dual (97.5%)"    = T1E_rob.dual.975,
  
  check.names = FALSE # Prevents R from changing spaces to dots
)

# Pivot longer and split the column names into "Prior" and "Decision"
df_plot <- df %>%
  pivot_longer(
    cols = -theta,
    names_to = c("Prior", "Decision"),
    names_sep = "_",
    values_to = "T1E"
  ) %>%
  # Optional: Convert Prior to a factor so it plots in a specific order in the legend
  mutate(Prior = factor(Prior, levels = c("Vague", "MAP", "Robust")))



saveRDS(df_plot, file = "data/T1E.fixed.rds")
cat("T1E fixed data saved")

# -------------------
# Plotting
# -------------------
ggplot(df_plot, aes(x = theta, y = T1E, color = Prior)) + 
  geom_line(linewidth = 1) + 
  geom_vline(xintercept = -50, linetype = "dashed", alpha = 0.5) + 
  geom_hline(yintercept = 0.035, linetype = "dotted") +
  
  # Facet by the Decision criteria we extracted
  facet_wrap(~ Decision, ncol = 2) + 
  
  labs(
    title = "Type I Error by Decision Criteria", 
    x = "True Mean CDAI Change from Baseline",
    y = "Type I Error",
    color = "Analysis Prior"
  ) +
  theme_bw() + # Adding a clean theme makes facets look much better
  ylim(0, 1)   # Consider changing this to e.g., c(0, 0.2) if you want to zoom in on the T1E inflation!





















# -----------------------------------
# EVALUATE Average T1E
# ----------------------------------

get_avg_t1e_table <- function(crit_name, succ.crit) {
  
  # define the 3 oc for the 3 analysis prior
  avgoc_vague <- avgoc2S.normMix(
    prior1 = p_vague, prior2 = p_vague, n1 = n.act, n2 = n.pbo, 
    decision = succ.crit, delta = 0, design_prior2 = p_vague, sigma1 = sigma, sigma2 = sigma
  ) 
  
  avgoc_MAP <- avgoc2S.normMix(
    prior1 = p_vague, prior2 = p_MAP, n1 = n.act, n2 = n.pbo, 
    decision = succ.crit, delta = 0, design_prior2 = p_vague, sigma1 = sigma, sigma2 = sigma
  ) 
  
  avgoc_rob <- avgoc2S.normMix(
    prior1 = p_vague, prior2 = p_rob, n1 = n.act, n2 = n.pbo, 
    decision = succ.crit, delta = 0, design_prior2 = p_vague, sigma1 = sigma, sigma2 = sigma
  ) 
  
  # Evaluate for the 4 design priors
  col_vague <- c(avgoc_vague(), 
                 avgoc_MAP(), 
                 avgoc_rob())
  
  col_skep  <- c(avgoc_vague(design_prior2_new = p_skep), 
                 avgoc_MAP(design_prior2_new = p_skep), 
                 avgoc_rob(design_prior2_new = p_skep))
  
  col_MAP   <- c(avgoc_vague(design_prior2_new = p_MAP), 
                 avgoc_MAP(design_prior2_new = p_MAP), 
                 avgoc_rob(design_prior2_new = p_MAP))
  
  col_rob   <- c(avgoc_vague(design_prior2_new = p_rob), 
                 avgoc_MAP(design_prior2_new = p_rob), 
                 avgoc_rob(design_prior2_new = p_rob))
  
  # Tale
  tab <- data.frame(
    Vague        = col_vague,
    Skeptical    = col_skep,
    MAP          = col_MAP,
    Robust_p_MAP = col_rob
  ) 
  
  
  
  tab$Analysis_Prior <- c("Vague", "MAP", "Robust")
  tab$Decision_Criteria <- crit_name
  
  # Reorder cols
  tab <- tab[, c("Decision_Criteria", "Analysis_Prior", "Vague", "Skeptical", "MAP", "Robust_p_MAP")]
  
  return(tab)
}

# -----------------------------------
# GENERATE THE 4 TABLES
# -----------------------------------
tab1 <- get_avg_t1e_table("Significance (95%)", sign.crit.95)
tab2 <- get_avg_t1e_table("Significance (97.5%)", sign.crit.975)
tab3 <- get_avg_t1e_table("Dual (95%)", dual.crit.95)
tab4 <- get_avg_t1e_table("Dual (97.5%)", dual.crit.975)

final_table <- bind_rows(tab1, tab2, tab3, tab4)








# -------------------------------------
# COMPUTE ASSURANCE 
# -------------------------------------


# just for chdebugging, we try a prior that is way off --> Normal equal to 20 individuals, with a mean of -10 
p_off  <- mixnorm(c(1, -10, 20), sigma = sigma, param = "mn") 


library(parallel)
library(pbmcapply)

# Prior lists,  SKEP excluded from analysis priors, included in design priors
analysis_priors <- list(
  MAP        = p_MAP,
  Robust_0.2 = p_rob0.2,
  Robust_0.5 = p_rob_0.5,
  Vague      = p_vague,
  Normal     = Normal,
  Off = p_off
)

design_priors <- list(
  MAP        = p_MAP,
  Robust_0.2 = p_rob0.2,
  Robust_0.5 = p_rob_0.5,
  Vague      = p_vague,
  Normal     = Normal,
  Skeptical  = p_skep,
  Off = p_off
)

deltas <- -seq(0, 80, by = 10)  

# Build jobs: one per analysis prior
jobs <- names(analysis_priors)

RNGkind("L'Ecuyer-CMRG")
set.seed(123)

results_list <- pbmclapply(jobs, function(ap_name) {

  ap <- analysis_priors[[ap_name]]

  # Build the avgoc object for this analysis prior
  # (delta = 0 as placeholder, will be overridden at evaluation)
  avgoc_obj <- avgoc2S.normMix(
    prior1         = prior.t,
    prior2         = ap,
    n1             = n.act,
    n2             = n.pbo,
    decision       = dual.crit.975,
    delta          = 0,
    design_prior2  = p_vague,   
    sigma1         = sigma,
    sigma2         = sigma
  )

  # Loop over all design priors and all deltas
  bind_rows(lapply(names(design_priors), function(dp_name) {
    dp <- design_priors[[dp_name]]
    bind_rows(lapply(deltas, function(d) {
      data.frame(
        Analysis_Prior    = ap_name,
        Design_Prior      = dp_name,
        Delta             = d,
        Power             = avgoc_obj(delta_new = d, design_prior2_new = dp)
      )
    }))
  }))

}, mc.cores = 5)

# Combine all results
df_power_full <- bind_rows(results_list) |>
  mutate(
    Analysis_Prior = factor(Analysis_Prior, levels = names(analysis_priors)),
    Design_Prior   = factor(Design_Prior,   levels = names(design_priors))
  )









saveRDS(final_table, file = "data/avgT1E.fixed.rds")
saveRDS(df_power_full, file = "data/Assurance.fixed.rds")
cat("Assurance data saved")

final_table <- readRDS("data/avgT1E.fixed.rds")








# PLOT FOR ASSURANCE for each design x analysis prios

df_power_full<- readRDS("data/Assurance.fixed.rds")
















# Plot
p_power_all <- df_power_full |>
  ggplot(aes(x = abs(Delta), y = Power,
             colour = Analysis_Prior,
             group  = Analysis_Prior)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  facet_wrap(~ Design_Prior, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 90, 20),
                     labels = seq(0, 90, 20)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    subtitle = "Facets = design prior  |  Colour = analysis prior",
    x        = "Delta",
    y        = "Power",
    colour   = "Analysis prior"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom"  )






# Plot
p_power_all <- df_power_full |>
  ggplot(aes(x = abs(Delta), y = Power,
             colour = Design_Prior,
             group  = Design_Prior)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.5) +
  geom_hline(yintercept = 0.80, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  facet_wrap(~ Analysis_Prior, ncol = 3) +
  scale_x_continuous(breaks = seq(0, 90, 20),
                     labels = seq(0, 90, 20)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(
    subtitle = "Facets = design prior  |  Colour = design prior",
    x        = "|Delta|",
    y        = "Power",
    colour   = "Design prior"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position  = "bottom"  )
