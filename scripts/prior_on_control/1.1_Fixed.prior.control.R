source("00.functions.R")


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
p_MAP <- mixnorm(c(0.51, -51, 19.9), c(0.44, -46.8, 7.6), c(0.05, -54.1, 51.7), sigma = sigma, param = "ms")

p_rob <- robustify(p_MAP, 0.2, mean = -50)

p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma)

p_skep <- mixnorm(c(1, -90, 25), sigma = sigma)






# --------------------
# CREATE OC
# --------------------
n.act <- 40
n.pbo <- 20

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


saveRDS(final_table, file = "data/avgT1E.fixed.rds")
cat("avgT1E fixed data saved")

final_table <- readRDS("data/avgT1E.fixed.rds")






# Plot experiments

colnames(final_table)[6] <- "Robust"

df_plot_avg <- final_table  |> 
  pivot_longer(
    cols = c("Vague", "Skeptical", "MAP", "Robust"),
    names_to = "Design_Prior",
    values_to = "Avg_T1E"
  )  |> 
  mutate(
    Analysis_Prior = factor(Analysis_Prior, levels = c("Vague", "MAP", "Robust")),
    Design_Prior = factor(Design_Prior, levels = c("Vague", "Skeptical", "MAP", "Robust")),
    Decision_Criteria = factor(Decision_Criteria, 
                               levels = c("Significance (95%)", "Significance (97.5%)", 
                                          "Dual (95%)", "Dual (97.5%)"))
  )
  # Faceted Bar Chart
  ggplot(df_plot_avg, aes(x = Design_Prior, y = Avg_T1E, fill = Analysis_Prior)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7, color = "black") +
    
    # Add a threshold line (Change 2.5 to 5.0 if you are targeting 5% alpha)
    geom_hline(yintercept = 0.25, linetype = "dashed", color = "red", linewidth = 0.8) +
    
    # Facet by the 4 Decision Criteria - 4 columns (2x2 layout: Sig95, Sig975, Dual95, Dual975)
    facet_wrap(~ Decision_Criteria, ncol = 4) +
    
    labs(
      title = "Average Type I Error Across Decision Criteria",
      subtitle = "Evaluating prior robustness against varying underlying true scenarios",
      x = "True Underlying Reality (Design Prior)",
      y = "Average Type I Error (%)",
      fill = "Analysis Prior"
    ) +
    scale_fill_manual(values = c("Vague" = "#E69F00", "MAP" = "#56B4E9", "Robust" = "#009E73")) +
    theme(
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "#f0f0f0"),
      figure.width = 14,
      figure.height = 4
    )














