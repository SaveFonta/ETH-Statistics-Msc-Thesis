# Script 2.1 
library(RBesT)
library(dplyr)
library(tidyr)
library(ggplot2)
library(checkmate)  
library(assertthat) 

# FIRST THING FIRST:
# Define Analysis priors and Design priors

# ---- Design parameters ---
sigma <- 88
n_T <- 40
n_C <- 20

# ----- Analysis priors -----

analysis_priors <- list (
    p_normal = mixnorm(c(1, -49, 20), sigma = sigma, param = "mn"),
    p_MAP  = mixnorm(c(0.51, -51, 19.9), c(0.44, -46.8, 7.6), c(0.05, -54.1, 51.7),
                  sigma = sigma, param = "ms"),
    p_rob.2  = robustify(p_MAP, 0.2, mean = -49), 
    p_rob.4  = robustify(p_MAP, 0.4, mean = -49),  # maybe once I will change this mean to adapt it to have the mean of data as Weru suggests, not today
    p_rob.6  = robustify(p_MAP, 0.6, mean = -49),
    p_rob.8  = robustify(p_MAP, 0.8, mean = -49),
    p_vague = mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms")
)
list2env(analysis_priors, envir = .GlobalEnv)

# ----- Design priors ----

design_priors <- list (
    p_dirac = mixnorm(c(1, -49, 1e-16), sigma = sigma, param = "ms"),
    p_normal = mixnorm(c(1, -49, 20), sigma = sigma, param = "mn"),
    p_MAP  = mixnorm(c(0.51, -51, 19.9), c(0.44, -46.8, 7.6), c(0.05, -54.1, 51.7),
                  sigma = sigma, param = "ms"),
    # p_vague = mixnorm(c(1, -50, 8800), sigma = sigma, param = "ms"),
    p_skep  = mixnorm(c(1, -90, 25), sigma = sigma, param = "ms"),  #posterior distribution of a stand-alone analysis of the historical study with the “most extreme” placebo effect
    p_off  = mixnorm(c(1, -10, 20), sigma = sigma, param = "mn") 
)
list2env(design_priors, envir = .GlobalEnv)



# --- Decision criteria ---
sign.crit <- decision2S(pc = 0.975, qc = 0, lower.tail = TRUE)
dual.crit <- decision2S(pc = c(0.975, 0.5), qc = c(0, -50), lower.tail = TRUE) #for the moment I will ignore this 






#################################################
#                 CLASSIC T1E
################################################ 
# It is P(success | theta_c, delta = 0) and it is a function of the real theta_c


# create functions 
oc_normal <- oc2S(prior1 = p_vague, prior2 = p_normal, n1 = n_T, n2 = n_C, decision = sign.crit)
oc_MAP <- oc2S(prior1 = p_vague, prior2 = p_MAP, n1 = n_T, n2 = n_C, decision = sign.crit)
oc_rob.2 <- oc2S(prior1 = p_vague, prior2 = p_rob.2, n1 = n_T, n2 = n_C, decision = sign.crit)
oc_rob.4 <- oc2S(prior1 = p_vague, prior2 = p_rob.4, n1 = n_T, n2 = n_C, decision = sign.crit)
oc_rob.6 <- oc2S(prior1 = p_vague, prior2 = p_rob.6, n1 = n_T, n2 = n_C, decision = sign.crit)
oc_rob.8 <- oc2S(prior1 = p_vague, prior2 = p_rob.8, n1 = n_T, n2 = n_C, decision = sign.crit)
oc_vague <- oc2S(prior1 = p_vague, prior2 = p_vague, n1 = n_T, n2 = n_C, decision = sign.crit)






# -------------------
# Evalute classic T1E
# -------------------
theta_c <- seq(-80, 0)
theta_t <- theta_c


T1E_normal <- oc_normal(theta_t, theta_c)
T1E_MAP <- oc_MAP(theta_t, theta_c)
T1E_rob.2 <- oc_rob.2(theta_t, theta_c)
T1E_rob.4 <- oc_rob.4(theta_t, theta_c)
T1E_rob.6 <- oc_rob.6(theta_t, theta_c)
T1E_rob.8 <- oc_rob.8(theta_t, theta_c)
T1E_vague <- oc_vague(theta_t, theta_c)


df_T1E <- data.frame(
      theta = theta_c,
      "Normal" = T1E_normal,
      "MAP" = T1E_MAP, 
      "Robust (w=0.2)" = T1E_rob.2, 
      "Robust (w=0.4)" = T1E_rob.4,
      "Robust (w=0.6)" = T1E_rob.6,
      "Robust (w=0.8)" = T1E_rob.8,
      "Vague" = T1E_vague,
        check.names = FALSE # Prevents R from changing spaces to dots
)

df_T1E.plot <- df_T1E  |> 
    pivot_longer(
        cols = -theta,
        names_to = "Prior", 
        values_to = "T1E"
    )

df_T1E.plot  |> 
ggplot(aes(x = theta, y = T1E, color = Prior))  +
geom_line(linewidth = 1) + 
geom_hline (yintercept = 0.025, linetype = "dashed", color = "red") + 
labs(
    x = "Theta",
y = "Type 1 Error",
    color = "Metodo Prior"
) + 
theme_minimal()


# COMMENTS:
# For the vague prior (prior 7), the curve is sitting exactly at the nominal α, this is the frequentist reference

# Historical data has a mean of around -49 (mean of the Normal prior one component, but very similar to the the MAP 3 components, since both come from 
#the same data) 
# For the MAP the curves gets further away from alfa with inflation of theta_c, so when there is conflict with historical data.

For the 1-component and 3-component MAP (priors 1-2), the curve departs substantially from α, with inflation where θ_c conflicts with the historical data
For robust MAP priors 3-6, the departure is moderated by ω, with higher ω pulling the curve back toward the vague baseline
The curves are not monotone in θ_c — there are regimes where borrowing helps (T1E below nominal) and regimes where it hurts (T1E above nominal)







# TRY TO PLOT POWER:
# sHOW THE 3d SURFACE, or can fix a delta and plot against theta_c
# Point is: neither T1E nor Power is a single number under Bayesian borrowing, so the compuations of classical EUII = DOR^(1/n) is not directly possible like the frequentist case.
# This motivates the avgT1E framework of Best et al. and the subsequent definition of a Bayesian EUII.













# So it is a mess to evaluate this blablabla 
# then define avgT1E, which is the formula with integral, note that the classic T1E is just an avergae T1E with dirac at a true value. 
# So now for ease of notation, we will always talk about avgT1E and avgPower (when delta fixed), since we don't want to put a prior on that since there is no knowdlede
# If we say avgT1E using dirac at something, then it is just the conditional T1E

# 