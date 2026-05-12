library(RBesT)



# deifne data
sigma <- 88
data <- crohn
data <- transform(data, y.se =sigma /sqrt(n))


# define possible priors
p_MAP <- mixnorm(c(0.51, -51, 19.9), c(0.44, -46.8, 7.6), c(0.05, -54.1, 51.7), sigma = sigma, param = "ms")
p_MAP



# --------------------------
# Can even do it automatically

map_mcmc <- gMAP(cbind(y, y.se) ~ 1 | study,
  weights = n, data = data,
  family = gaussian,
  beta.prior = cbind(0, sigma),
  tau.dist = "HalfNormal", tau.prior = cbind(0, sigma / 2)
)

map <- automixfit(map_mcmc)
print(map)
plot(map)$mix

round(ess(map)) ## default elir method


# Now lets fit only one Normal
one <- mixfit(map_mcmc, Nc = 1) 
round(ess(one)) #almost 20 patients, not bad, similar to the 20 from Gsponer .



# --------------------------------------

p_rob <- robustify(p_MAP, 0.2, mean = -50)

p_vague <- mixnorm(c(1, -50, 8800), sigma = sigma)



# define success criterion 
succ.crit <- decision2S(pc = 0.95, qc = 0, lower.tail = TRUE)





# --------------------
# CREATE OC
# --------------------
n.act <- 40
n.pbo <- 20

oc_vague <- oc2S(prior1 = p_vague, prior2 = p_vague, n1 = n.act, n2 = n.pbo, decision = succ.crit)
oc_MAP <- oc2S(prior1 = p_vague, prior2 = p_MAP, n1 = n.act, n2 = n.pbo, decision = succ.crit)
oc_rob <- oc2S(prior1 = p_vague, prior2 = p_rob, n1 = n.act, n2 = n.pbo, decision = succ.crit)



# -------------------
# Evalute classic T1E
# -------------------
theta_c <- seq(-150, 50)
theta_a <- theta_c

T1E_vague <- oc_vague(theta_c, theta_a)
T1E_MAP <- oc_MAP(theta_c, theta_a)
T1E_rob <- oc_rob(theta_c, theta_a)

df <- data.frame (theta = theta_c, T1E_vague = T1E_vague, T1E_MAP = T1E_MAP, T1E_rob= T1E_rob )

library(dplyr)
library(tidyr)

df_plot <- pivot_longer(df, cols = c(T1E_vague, T1E_MAP, T1E_rob), names_to = "prior")


library(ggplot2)
ggplot(df_plot, aes(x = theta, y = value, color = prior)) + 
    geom_line(linewidth = 1) + 
    geom_vline (xintercept = -50, linetype = "dashed", alpha = 0.5) + 
    geom_hline (yintercept = 0.035, linetype = "dotted") +
    labs (
        title = "Classic T1E", 
x = "True Mean CDAI Change from Baseline",
    y = "T1E",
    color = "Analysis Prior"
    ) +
    ylim(0,1)




# -----------------------------------
# EVALUATE Average T1E
# -----------------------------------

# WRONG APPROACH:
pos_vague <- pos2S(prior1 = p_vague, prior2 = p_vague, n1 = n.act, n2 = n.pbo, decision = succ.crit)
pos_MAP <- pos2S(prior1 = p_vague, prior2 = p_MAP, n1 = n.act, n2 = n.pbo, decision = succ.crit)
pos_rob <- pos2S(prior1 = p_vague, prior2 = p_rob, n1 = n.act, n2 = n.pbo, decision = succ.crit) 


avg.T1E_vague <- pos_vague(p_vague, p_vague)
avg.T1E_MAP <- pos_MAP(p_vague, p_MAP)
avg.T1E_rob <- pos_rob(p_vague, p_rob)




# no this doesnt work, since like this we are assuming the two priors are completely indepedent 
# we need to sample and then plug those values 


# SIMULATION APPROACH
samples_vague <- rmix(p_vague, n = 10000)
samples_MAP   <- rmix(p_MAP, n = 10000)
samples_rob   <- rmix(p_rob, n = 10000)
 



pointwise_T1E_vague <- oc_vague(samples_vague, samples_vague)
pointwise_T1E_MAP   <- oc_MAP(samples_MAP, samples_MAP)
pointwise_T1E_rob   <- oc_rob(samples_rob, samples_rob)

mean(pointwise_T1E_vague)
mean(pointwise_T1E_MAP)
mean(pointwise_T1E_rob)



# Ok that we can simulate, but we can also define the integral!
#our function is integral[P(reject | theta_c = theta_a) * p(theta_c)] d theta_c



# -----------------------------------
# EVALUATE Average T1E for real
# -----------------------------------
average_T1E <- function(analysis.criteria, design.prior, lower = NA, upper = NA, n_nodes = 128, width_max = 2000) {

    
    # Extract the extremes for integration 
    if (is.na(lower) && is.na(upper)) { 

        # find the x values that cover 99.99 prob mass
        bounds <- qmix(design.prior, p = c(1e-5, 1 - 1e-5))
        lower <- bounds[1]
        upper <- bounds[2]
    }

    width <- upper - lower


    #create the function to integrate 
    to.integrate <- function(x) analysis.criteria(x,x) * dmix(design.prior, x)

    # integrate 
    if (width > width_max) {  # Gauss-Legendre 
    gl <- statmod::gauss.quad(n = n_nodes, kind = "legendre")
    
    # Scale nodes to our window
    mid  <- (upper + lower) / 2
    half.width <- (upper - lower) / 2
    x_nodes <- mid + half.width * gl$nodes
    
    # Weighted Summation
    fvals <- to.integrate(x_nodes)
    res_val <- half.width * sum(gl$weights * fvals)
    
    cat("Gauss-Legendre used with", n_nodes, "nodes since width larger ", width_max)
  }

  else { # classic integrate 
    res <- integrate(to.integrate, lower = lower, upper = upper)
    res_val <- res$value
  }
  
  return(res_val)
}






# Lets recreate the table 2
# we have an additional design prior, which is skeptical 
p_skep <- mixnorm(c(1, -90, 25), sigma = sigma)

T1E_vague_row <- c(
  average_T1E(oc_vague, p_vague),
  average_T1E(oc_vague, p_skep),
  average_T1E(oc_vague, p_MAP),
  average_T1E(oc_vague, p_rob)
)


T1E_MAP_row <- c(
  average_T1E(oc_MAP, p_vague),
  average_T1E(oc_MAP, p_skep),
  average_T1E(oc_MAP, p_MAP),
  average_T1E(oc_MAP, p_rob)
)

T1E_rob_row <-  c(
  average_T1E(oc_rob, p_vague),
  average_T1E(oc_rob, p_skep),
  average_T1E(oc_rob, p_MAP),
  average_T1E(oc_rob, p_rob)
)

Table_2 <- rbind(
  T1E_vague_row,
  T1E_MAP_row,
  T1E_rob_row
)
Table_2 <- as.data.frame(Table_2)



colnames(Table_2) <- c("Design_Vague", "Design_Skeptical", "Design_MAP", "Design_Robust")
rownames(Table_2) <- c("Analysis_Vague", "Analysis_MAP", "Analysis_Robust")
Table_2_Formatted <- round(Table_2 * 100, 1)


print(Table_2_Formatted)










library(checkmate)   # assert_number, assert_scalar, assert_list, assert_function
library(assertthat)  # assert_that

avgoc_vague <- avgoc2S.normMix(
  prior1 = p_vague, # analysis prior
  prior2 = p_vague,
  n1 = n.act,
  n2 = n.pbo,
  decision = succ.crit,
  delta = 0,        # <-- theta1 = theta2 + delta
  mix2 = p_vague,         # <- ( design prior)
  sigma1 = sigma,
  sigma2 = sigma
) 


avgoc_MAP <- avgoc2S.normMix(
  prior1 = p_vague, 
  prior2 = p_MAP,       # analysis prior
  n1 = n.act,
  n2 = n.pbo,
  decision = succ.crit,
  delta = 0,        # <-- theta1 = theta2 + delta
  mix2 = p_vague,         # <- ( design prior)
  sigma1 = sigma,
  sigma2 = sigma
) 


avgoc_rob <- avgoc2S.normMix(
  prior1 = p_vague, 
  prior2 = p_rob,       # analysis prior
  n1 = n.act,
  n2 = n.pbo,
  decision = succ.crit,
  delta = 0,        # <-- theta1 = theta2 + delta
  mix2 = p_vague,         # <- ( design prior)
  sigma1 = sigma,
  sigma2 = sigma
) 


col1 <- c(avgoc_vague(), avgoc_MAP(), avgoc_rob())

col2 <- c(avgoc_vague(mix2_new = p_skep), avgoc_MAP(mix2_new = p_skep), avgoc_rob(mix2_new = p_skep))

col3 <- c(avgoc_vague(mix2_new = p_MAP), avgoc_MAP(mix2_new = p_MAP), avgoc_rob(mix2_new = p_MAP))

col4 <- c(avgoc_vague(mix2_new = p_rob), avgoc_MAP(mix2_new = p_rob), avgoc_rob(mix2_new = p_rob))



full <- cbind(col1, col2, col3, col4)

tab <- as.data.frame(round(full,3))
rownames(tab) <- c("Vague", "MAP", "Robust p_MAP")
colnames(tab) <- c("Vague", "Skeptical", "MAP", "Robust p_MAP")
final <- tab * 100
final



















































































