# --------------------------------------------
# In here I want to understand how does the EUII behaves in a frequentist dual criterion design.
# -----------------------------------------------

# Need to create a dual criterion design. Then see how the EUII behaves for varying underlying truths. 

# nMin is the treshold of sample size. If we have a saple size bigger than nMin, then the hardest treshorld to overcome is DV,
# if n smaller nmin it is harder to get significance. 

# EXPLANATION
# significance if theta > sign_tresh, where sign_tres = NV + z_(1-alpha) * sigma / sqrt(n)
# relevance if theta > clin_tresh, where clin_tresh = DV

#Set the two treshold equal --> clin_tresh = sign_tresh  --> DV = NV + z_(1-alpha) * sigma / sqrt(n)

# Simplify and the n required is n =  [(z_(1-alpha) * sigma) / (DV - NV)]^2. 
# here we assume NV = 0


nMin <- function(DV, sigma=1, alpha=0.025, ceil=FALSE) {
  res <- (qnorm(1 - alpha) * sigma / DV)^2
  if(ceil) return(ceiling(res))
  return(res)
}

# Professor's Functions
# oneMinus = FALSE ---> T1E
# oneMinus = TRUE --> 1-T1E = Specificity 


DC.T1Erate_ <- function(n, DV, alpha=0.025, sigma=1, log.p=FALSE, oneMinus=FALSE) {
  nmin <- nMin(DV=DV, sigma=sigma, alpha=alpha, ceil=FALSE)
  if(!oneMinus) {
    if(!log.p) {
      if(n <= nmin) res <- alpha
      else res <- pnorm(-sqrt(n)*DV/sigma)
    }
    if(log.p) {
      if(n <= nmin) res <- log(alpha) #if the significant tresh is harder, then t1e is classic alpha
      else res <- pnorm(-sqrt(n)*DV/sigma, log.p=TRUE) # otw t1e is P_0(theta > DV). But under H0 theta N (0, sigma^2 / n). So Z = DV - 0 / sigma / sqrt(n) = sqrtn DV / sigma. CDF(-Z). Finds the prob of having a value bigger, so rejecting  
    }
  }
  if(oneMinus) {
    if(!log.p) {
      if(n <= nmin) res <- 1 - alpha
      else res <- pnorm(-sqrt(n)*DV/sigma, lower.tail=FALSE)
    }
    if(log.p) {
      if(n <= nmin) res <- log(1 - alpha)
      else res <- pnorm(-sqrt(n)*DV/sigma, log.p=TRUE, lower.tail=FALSE)
    }
  }
  return(res)
}

DC.Power_ <- function(n, theta, DV, sigma=1, alpha=0.025, log.p=FALSE, oneMinus=FALSE) {
    nmin <- nMin(DV=DV, sigma=sigma, alpha=alpha, ceil=FALSE)

    #use formula from my sheet:
  if (n <= nmin) { # Significance is the threshold to pass
      k <- qnorm(1 - alpha) * sigma / sqrt(n)
    } else {         # DV is the threshold to pass
      k <- DV
    }
# k <- pmax (qnorm(1 - alpha) * sigma / sqrt(n), DV) would be the same


  term <- sqrt(n) * (theta - k) / sigma

  if(!oneMinus) { # power
    if(!log.p) res <- pnorm(term)
    if(log.p) res <- pnorm(term, log.p=TRUE) #P_1(theta_hat > DV)
  }
  if(oneMinus) { # 1-power
    if(!log.p) res <- pnorm(term, lower.tail=FALSE)
    if(log.p) res <- pnorm(term, log.p=TRUE, lower.tail=FALSE)
  }
  return(res)
}


DC.Power <- Vectorize(DC.Power_)
DC.T1Erate <- Vectorize(DC.T1Erate_)

# Setup Parameters for Visualization
myDV <- 0.3
mytheta <- c(-0.1, 0, 0.1, 0.3, 0.5, 0.8) # Different true effects (note power is 50% max when true effect = DV) 
# QUESTION 1: does it even makes sense to look at the behavior of EUII when theta is smaller than the DV? 

nmin <- nMin(DV=myDV)       # Calculates to approx 42.68

# Create a sequence of sample sizes from 10 up to 2000
myn <- 2^(seq(log2(35), log2(2000), length.out=30)) #This creates an exponentially increasing gap between numbers.


# Note: professor code was starting the sequence at nmin. Ofc, if n < n_min, the EUII behaves like it is in a classic test,



# Calculate EUII
logT1Erate <- DC.T1Erate(n=myn, DV=myDV, log.p=TRUE)
logOneMinusT1Erate <- DC.T1Erate(n=myn, DV=myDV, log.p=TRUE, oneMinus=TRUE)

EUII <- matrix(NA, nrow=length(myn), ncol=length(mytheta))

for(i in 1:length(mytheta)) {
  logPower <- DC.Power(n=myn, theta=mytheta[i], DV=myDV, log.p=TRUE)
  logOneMinusPower <- DC.Power(n=myn, theta=mytheta[i], DV=myDV, log.p=TRUE, oneMinus=TRUE)
  
  logDOR <- logPower - logOneMinusPower - logT1Erate + logOneMinusT1Erate
  
  EUII[,i] <- exp(logDOR / myn)
}

par(las=1, mfrow=c(1,1))
matplot(myn, EUII, type="b", pch=19, log="x", lty=1,
        col=c("orange", "purple", "black", "red", "blue", "forestgreen"),  # add colors here
        xlab="Sample Size (n)", ylab="EUII",
        main="EUII in Dual-Criterion Design (DV = 0.3)")
abline(v=nmin, lty=2, col="gray50")
text(nmin, max(EUII, na.rm=TRUE)*0.9, labels="n_min", pos=4, col="gray50")
legend("topright", lty=1,
       col=c("orange", "purple", "black", "red", "blue", "forestgreen"),   # add here too
       legend=paste("theta =", mytheta), lwd=2)


# WHAT I DISCOVERED: 
# denote theta = true value
# if 0 < theta < DV, then EUII converges to 1 
# if theta < 0, then EUII is negative. This happens cause Power<T1E, biased test, and DOR is below 1. 




#Add the asymptotic limits I am not completely sure of
calc_DC_limit <- function(theta, DV, sigma=1) {
  if (theta <= DV) {
    return( exp( (theta * (2*DV - theta)) / (2*sigma^2) ) )
  } else {
    return( exp( ((theta - DV)^2 + DV^2) / (2*sigma^2) ) )
  }
}

# Apply the function to all of chosen thetas
DC_limits <- sapply(mytheta, calc_DC_limit, DV=myDV)

mycolors <- c("orange", "purple", "black", "red", "blue", "forestgreen")

abline(h=DC_limits, col=mycolors, lty=3, lwd=2) #lty = 3 is dotted line











# ------------------------------------------------
# Dual criterion vs standard
# ------------------------------------------------

# For a fixed DV = 0.4, which design has an higher EUII for different the true effects? 


#  Standard Design Functions
Standard.T1Erate_ <- function(n, alpha=0.025, log.p=FALSE, oneMinus=FALSE) {
  if(!oneMinus) {
    if(!log.p) res <- alpha
    if(log.p) res <- log(alpha)
  } else {
    if(!log.p) res <- 1 - alpha
    if(log.p) res <- log(1 - alpha)
  }
  return(res)
}

Standard.Power_ <- function(n, theta, sigma=1, alpha=0.025, log.p=FALSE, oneMinus=FALSE) {
  z_alpha <- qnorm(1 - alpha)
  term <- sqrt(n) * theta / sigma - z_alpha
  
  if(!oneMinus) {
    if(!log.p) res <- pnorm(term)
    if(log.p) res <- pnorm(term, log.p=TRUE)
  } else {
    if(!log.p) res <- pnorm(term, lower.tail=FALSE)
    if(log.p) res <- pnorm(term, log.p=TRUE, lower.tail=FALSE)
  }
  return(res)
}

Standard.T1Erate <- Vectorize(Standard.T1Erate_)
Standard.Power <- Vectorize(Standard.Power_)


##
myDV <- 0.3
mytheta <- c(0.3, 0.45, 0.6) 

nmin <- (qnorm(1 - 0.025) * 1 / myDV)^2 # approx 42.68

myn <- 2^(seq(log2(10), log2(2000), length.out=50))

# Matrices to store results
EUII_DC <- matrix(NA, nrow=length(myn), ncol=length(mytheta))
EUII_Std <- matrix(NA, nrow=length(myn), ncol=length(mytheta))


#Lets also store DOR
DOR_DC <- matrix(NA, nrow=length(myn), ncol=length(mytheta))
DOR_Std <- matrix(NA, nrow=length(myn), ncol=length(mytheta))



# T1E for Standard is constant across all thetas
logT1E_Std <- Standard.T1Erate(n=myn, log.p=TRUE)
logOneMinusT1E_Std <- Standard.T1Erate(n=myn, log.p=TRUE, oneMinus=TRUE)

# T1E for DC
logT1E_DC <- DC.T1Erate(n=myn, DV=myDV, log.p=TRUE)
logOneMinusT1E_DC <- DC.T1Erate(n=myn, DV=myDV, log.p=TRUE, oneMinus=TRUE)

for(i in 1:length(mytheta)) {
  # Standard Power
  logPow_Std <- Standard.Power(n=myn, theta=mytheta[i], log.p=TRUE)
  logOneMinusPow_Std <- Standard.Power(n=myn, theta=mytheta[i], log.p=TRUE, oneMinus=TRUE)
  
  logDOR_Std <- logPow_Std - logOneMinusPow_Std - logT1E_Std + logOneMinusT1E_Std
  DOR_Std[,i] <- exp(logDOR_Std)                 
  EUII_Std[,i] <- exp(logDOR_Std / myn)
  
  # DC Power
  logPow_DC <- DC.Power(n=myn, theta=mytheta[i], DV=myDV, log.p=TRUE)
  logOneMinusPow_DC <- DC.Power(n=myn, theta=mytheta[i], DV=myDV, log.p=TRUE, oneMinus=TRUE)
  
  logDOR_DC <- logPow_DC - logOneMinusPow_DC - logT1E_DC + logOneMinusT1E_DC
  DOR_DC[,i] <- exp(logDOR_DC)                   
  EUII_DC[,i] <- exp(logDOR_DC / myn)
}

# Plotting the Comparison
colors <- c("blue", "forestgreen", "orange")


par(las=1, mfrow=c(1,1), mar=c(5, 4, 4, 2) + 0.1) 

matplot(myn, EUII_Std, type="l", lty=2, lwd=2, col=colors, log="x",
        ylim=range(c(EUII_Std, EUII_DC), na.rm=TRUE),
        xlab="Sample Size (n)", ylab="EUII", 
        main="EUII: Dual-Criterion vs Standard Design (DV = 0.4)")

# Add the Dual-Criterion lines and n_min
matlines(myn, EUII_DC, type="l", lty=1, lwd=2, col=colors)
abline(v=nmin, col="gray50", lty=3)

legend("topright", title="True Effect",
       legend=paste("theta =", mytheta), col=colors, lwd=2, lty=1)

legend("topleft", title="Design",
       legend=c("Dual-Criterion", "Standard"), col="black", lwd=2, lty=c(1, 2))



DC_limits <- sapply(mytheta, calc_DC_limit, DV=myDV)
abline(h = DC_limits,lty=3, lwd=2)






calc_Std_limits <- function (theta){ 
  exp(theta^2 / 2)
}
Std_limits <- sapply(mytheta, calc_Std_limits)
abline(h = Std_limits,lty=3, lwd=1)


#Great! makes sense.
# when n is smaller than n_min, we are just in the standard design since the harder treshold is the significance. So the EUIIs overlaps perfectu

# INTERESTINGLY, for some specific thetas that are not too larger than the DV, at some fixed n the DC has
# larger EUII than standard!!!!!!








# --- EUII ratio easier to spot!!!!!
EUII_ratio <- EUII_DC / EUII_Std

matplot(myn, EUII_ratio, type="l", lty=1, lwd=2, col=colors, log="x",
        xlab="Sample Size (n)", ylab="DC / Standard",
        main="EUII Ratio (DC / Standard)")
abline(h=1, lty=2, col="gray50")
abline(v=nmin, lty=3, col="gray50")
text(nmin, max(EUII_ratio)*0.95, "n_min", pos=4, col="gray50")
legend("topright", legend=paste("theta =", mytheta),
       col=colors, lwd=2, lty=1)



# It is also clear that increasing n the ratio of EUII_DC / EUII_STnd converges to something since it has asyntotes. 
# Useful obtaining them ? Maybe not or maybe yes

# I think more interestingly I could try to obtain what is the condition s.t. 
# EUII_DC > EUII_STnd. Since it is fixed I could more easily prove it for DOR





# -----------------------------------------
# Comparison at Fixed n 
# -------------------------------------------


n_fixed <- 100
myDV <- 0.4
sigma <- 1
theta_test <- c(0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9)

#  Standard Design 
logPow_Std <- Standard.Power(n=n_fixed, theta=theta_test, log.p=TRUE)
log1mPow_Std <- Standard.Power(n=n_fixed, theta=theta_test, log.p=TRUE, oneMinus=TRUE)
logT1E_Std <- Standard.T1Erate(n=n_fixed, log.p=TRUE)
log1mT1E_Std <- Standard.T1Erate(n=n_fixed, log.p=TRUE, oneMinus=TRUE)

logDOR_Std <- logPow_Std - log1mPow_Std - logT1E_Std + log1mT1E_Std
EUII_Std <- exp(logDOR_Std / n_fixed)

# Convert log Power and T1E back to normal probabilities for the table
Pow_Std <- exp(logPow_Std)
T1E_Std <- exp(logT1E_Std)[1] # Same for all thetas

# DC Design
logPow_DC <- DC.Power(n=n_fixed, theta=theta_test, DV=myDV, log.p=TRUE)
log1mPow_DC <- DC.Power(n=n_fixed, theta=theta_test, DV=myDV, log.p=TRUE, oneMinus=TRUE)
logT1E_DC <- DC.T1Erate(n=n_fixed, DV=myDV, log.p=TRUE)
log1mT1E_DC <- DC.T1Erate(n=n_fixed, DV=myDV, log.p=TRUE, oneMinus=TRUE)

logDOR_DC <- logPow_DC - log1mPow_DC - logT1E_DC + log1mT1E_DC
EUII_DC <- exp(logDOR_DC / n_fixed)

Pow_DC <- exp(logPow_DC)
T1E_DC <- exp(logT1E_DC)[1]

# table
OC_Table <- data.frame(
  True_Theta = theta_test,
  
  # Standard Columns
  Std_Power = round(Pow_Std, 4),
  Std_EUII = round(EUII_Std, 4),
  
  # DC Columns
  DC_Power = round(Pow_DC, 4),
  DC_EUII = round(EUII_DC, 4)
)

cat("=========================================================\n")
cat("   OPERATING CHARACTERISTICS AT n =", n_fixed, "(DV =", myDV, ")\n")
cat("   Standard Type-I Error :", T1E_Std, "\n")
cat("   Dual-Crit Type-I Error:", formatC(T1E_DC, format="e", digits=2), "\n")
cat("=========================================================\n")
print(OC_Table)


# Good that power is 0.5 at theta = DV.
# Interestinlgy, for moderate effects, the T1E inflation lead to a higher EUII for DC, 
# but for large effect the standard is better 













































# -------------------------------------------------------------------------
# EUII Analysis: Fixed True Effect, Varying (DV)
# -------------------------------------------------------------------------

target_theta <- 0.5   # The true underlying effect is fixed
sigma <- 1
alpha <- 0.025

# We test varying levels of strictness for the Dual-Criterion
# Note: we keep DV < target_theta, otherwise power < 50% asymptotically
myDVs <- c(0.2, 0.3, 0.4, 0.45) 
colors <- c("dodgerblue", "purple", "firebrick", "darkorange")

myn <- 2^(seq(log2(10), log2(5000), length.out=100))

# 1. Standard Design Trajectory (Only need to calculate once)
logPow_Std <- Standard.Power(n=myn, theta=target_theta, log.p=TRUE)
log1mPow_Std <- Standard.Power(n=myn, theta=target_theta, log.p=TRUE, oneMinus=TRUE)
logT1E_Std <- Standard.T1Erate(n=myn, log.p=TRUE)
log1mT1E_Std <- Standard.T1Erate(n=myn, log.p=TRUE, oneMinus=TRUE)

EUII_Std <- exp((logPow_Std - log1mPow_Std - logT1E_Std + log1mT1E_Std) / myn)

# 2. Dual-Criterion Trajectories
EUII_DC <- matrix(NA, nrow=length(myn), ncol=length(myDVs))
nmin_vals <- numeric(length(myDVs))

for(j in 1:length(myDVs)) {
  current_DV <- myDVs[j]
  nmin_vals[j] <- nMin(DV=current_DV, sigma=sigma, alpha=alpha, ceil=FALSE)
  
  logPow_DC <- DC.Power(n=myn, theta=target_theta, DV=current_DV, log.p=TRUE)
  log1mPow_DC <- DC.Power(n=myn, theta=target_theta, DV=current_DV, log.p=TRUE, oneMinus=TRUE)
  logT1E_DC <- DC.T1Erate(n=myn, DV=current_DV, log.p=TRUE)
  log1mT1E_DC <- DC.T1Erate(n=myn, DV=current_DV, log.p=TRUE, oneMinus=TRUE)
  
  EUII_DC[,j] <- exp((logPow_DC - log1mPow_DC - logT1E_DC + log1mT1E_DC) / myn)
}

# -------------------------------------------------------------------------
# Plotting
# -------------------------------------------------------------------------
par(mfrow=c(1,2), las=1, mar=c(5, 4, 4, 2) + 0.1)

# Plot 1: Absolute EUII
matplot(myn, EUII_DC, type="l", lty=1, lwd=2, col=colors, log="x",
        ylim=range(c(EUII_Std, EUII_DC), na.rm=TRUE),
        xlab="Sample Size (n)", ylab="EUII",
        main=paste("EUII (Fixed Theta =", target_theta, ")"))

# Add Standard Design
lines(myn, EUII_Std, col="black", lty=2, lwd=3)

# Add n_min indicators
for(j in 1:length(myDVs)) {
  abline(v=nmin_vals[j], col=colors[j], lty=3, lwd=1.5)
}

legend("topright", title="Design",
       legend=c("Standard", paste("DC (DV =", myDVs, ")")),
       col=c("black", colors), lty=c(2, 1, 1, 1, 1), lwd=2)

# Plot 2: EUII Ratio (DC / Standard)
EUII_ratio <- sweep(EUII_DC, 1, EUII_Std, "/") # Fast way to divide columns by vector

matplot(myn, EUII_ratio, type="l", lty=1, lwd=2, col=colors, log="x",
        xlab="Sample Size (n)", ylab="Ratio (EUII_DC / EUII_Std)",
        main="Efficiency Ratio: DC vs Standard")
abline(h=1, col="gray50", lty=2, lwd=2)

for(j in 1:length(myDVs)) {
  abline(v=nmin_vals[j], col=colors[j], lty=3, lwd=1.5)
}

legend("topright", title="DV Threshold",
       legend=paste("DV =", myDVs), col=colors, lwd=2, lty=1)
































































# -----------------------------------------
# Comparison at Fixed Power 
# -------------------------------------------
target_power <- 0.80
alpha <- 0.025
sigma <- 1

#function for required sample size at a specific Power
required_n_std <- function(sigma = 1, theta, alpha, power) { 
  term1 <- qnorm(1-alpha)
  term2 <- qnorm(power)
  n <- ((term1 + term2) *sigma / theta)^2
  n
}

required_n_DV  <- function(sigma = 1, theta, alpha, power, DV) { 
  n_std <- required_n_std(theta = theta, alpha = alpha, power = power)
  n_relevance <- ((qnorm(power) * sigma ) / (theta - DV))^2
  n <- max(n_std,  n_relevance)
  n
}









theta_seq <- seq(0.21, 0.8, length.out = 150)

# The different designs we want to compare
my_DVs <- c(0.2, 0.3, 0.4)
colors <- c("dodgerblue", "firebrick", "forestgreen")

n_std_vec <- rep(NA, length(theta_seq))
EUII_std_vec <- rep(NA, length(theta_seq))

n_dc_mat <- matrix(NA, nrow=length(theta_seq), ncol=length(my_DVs))
EUII_dc_mat <- matrix(NA, nrow=length(theta_seq), ncol=length(my_DVs))

for(i in 1:length(theta_seq)) {
  current_theta <- theta_seq[i]
  
  # Standard
  n_std <- required_n_std(sigma=sigma, theta=current_theta, alpha=alpha, power=target_power)
  n_std_vec[i] <- n_std
  
  logPow_Std <- Standard.Power(n=n_std, theta=current_theta, log.p=TRUE)
  log1mPow_Std <- Standard.Power(n=n_std, theta=current_theta, log.p=TRUE, oneMinus=TRUE)
  logT1E_Std <- Standard.T1Erate(n=n_std, log.p=TRUE)
  log1mT1E_Std <- Standard.T1Erate(n=n_std, log.p=TRUE, oneMinus=TRUE)
  
  EUII_std_vec[i] <- exp((logPow_Std - log1mPow_Std - logT1E_Std + log1mT1E_Std) / n_std)
  
  # DC
  for(j in 1:length(my_DVs)) {
    current_DV <- my_DVs[j]
    
    # ONLY calculate if the true effect is strictly greater than the tresghold
    if(current_theta > current_DV) {
      n_dc <- required_n_DV(sigma=sigma, theta=current_theta, alpha=alpha, power=target_power, DV=current_DV)
      n_dc_mat[i, j] <- n_dc
      
      logPow_DC <- DC.Power(n=n_dc, theta=current_theta, DV=current_DV, log.p=TRUE)
      log1mPow_DC <- DC.Power(n=n_dc, theta=current_theta, DV=current_DV, log.p=TRUE, oneMinus=TRUE)
      logT1E_DC <- DC.T1Erate(n=n_dc, DV=current_DV, log.p=TRUE)
      log1mT1E_DC <- DC.T1Erate(n=n_dc, DV=current_DV, log.p=TRUE, oneMinus=TRUE)
      
      EUII_dc_mat[i, j] <- exp((logPow_DC - log1mPow_DC - logT1E_DC + log1mT1E_DC) / n_dc)
    }
  }
}


par(mfrow=c(1,2), las=1, mar=c(5, 4, 4, 2) + 0.1)

# Plot 1: n
plot(theta_seq, n_std_vec, type="l", col="black", lty=2, lwd=3,
     ylim=c(0, 300), xlim =c(0.19, 0.8),
     xlab="Target True Effect (Theta)", ylab="Required Sample Size (n)",
     main="Sample Size for 80% Power")

for(j in 1:length(my_DVs)) {
  lines(theta_seq, n_dc_mat[, j], col=colors[j], lty=1, lwd=2)
}

legend("topright", title="Design",
       legend=c("Standard", paste("DC (DV =", my_DVs, ")")),
       col=c("black", colors), lty=c(2, 1, 1, 1), lwd=c(3, 2, 2, 2))


# Plot 2: EUII
plot(theta_seq, EUII_std_vec, type="l", col="black", lty=2, lwd=3,
     ylim=c(1, 1.6), xlim =c(0.19, 0.8),
     xlab="Target True Effect (Theta)", ylab="EUII",
     main="Efficiency at 80% Power")

for(j in 1:length(my_DVs)) {
  lines(theta_seq, EUII_dc_mat[, j], col=colors[j], lty=1, lwd=2)
}

legend("bottomright", title="Design",
       legend=c("Standard", paste("DC (DV =", my_DVs, ")")),
       col=c("black", colors), lty=c(2, 1, 1, 1), lwd=c(3, 2, 2, 2))



# Insights: 
# Keeping the Power cosntant, as we increase the DV the EUII of DC is lower ffor a fixed theta, we get a smaller EUII, 
# Note that quite fast the DV fall back to standard design

# and note also that for a fixed Power, the standard EUII is always larger


















# DONT SHOW THIS TO LEO. La proof la faremo piu avanti


# -------------------------------------------------------------------------
# Simulation: Diagnostic Odds Ratio (DOR) of Standard vs. Dual-Criterion 
# Design as a function of the true effect size (theta).
# -------------------------------------------------------------------------

# 1. Define Design Parameters
se <- 1.0                              # Standard error (normalized to 1)
alpha <- 0.025                         # One-sided Type-I error rate
c_alpha <- qnorm(1 - alpha) * se       # Standard statistical threshold (~1.96)
dv <- 3.0                              # Clinical Decision Value (Must be > c_alpha)

# 2. Define a range of true effect sizes (theta)
# We need to go high enough to see the "blockbuster" effect where the proof applies
theta_seq <- seq(0.1, 7, by = 0.05)

# 3. Function to calculate the Diagnostic Odds Ratio (DOR)
calculate_dor <- function(threshold, theta, se) {
  # Power: Probability of passing the threshold given true effect theta
  power <- 1 - pnorm(threshold, mean = theta, sd = se)
  
  # Type-I Error (T1E): Probability of passing given no effect (theta = 0)
  t1e <- 1 - pnorm(threshold, mean = 0, sd = se)
  
  # Calculate odds
  power_odds <- power / (1 - power)
  t1e_odds <- t1e / (1 - t1e)
  
  # DOR is the ratio of Power Odds to T1E Odds
  return(power_odds / t1e_odds)
}

# 4. Calculate DOR across all theta values for both designs
dor_stnd <- calculate_dor(threshold = c_alpha, theta = theta_seq, se = se)
dor_dc   <- calculate_dor(threshold = dv, theta = theta_seq, se = se)

# 5. Plotting the results
# Using a log scale for the Y-axis (DOR) because the values explode exponentially
plot(theta_seq, dor_dc, type = "l", col = "darkgreen", lwd = 2, log = "y",
     xlab = expression(paste("True Effect Size (", theta, ")")),
     ylab = "Diagnostic Odds Ratio (Log Scale)",
     main = "DOR: Dual-Criterion vs Standard Design",
     ylim = range(c(dor_stnd, dor_dc), na.rm = TRUE))

lines(theta_seq, dor_stnd, col = "blue", lwd = 2, lty = 2)

# 6. Add the theoretical crossing point derived in our proof
crossing_point <- dv + c_alpha
abline(v = crossing_point, col = "red", lwd = 2, lty = 3)

# Add annotations to the plot
legend("topleft", 
       legend = c("Dual-Criterion (Threshold = DV)", 
                  "Standard Design (Threshold = c_alpha)",
                  paste0("Proof Crossing Point: theta = ", round(crossing_point, 2))),
       col = c("darkgreen", "blue", "red"), 
       lty = c(1, 2, 3), lwd = 2, cex = 0.9)

# 7. Print a summary check to the console
cat("--- Verification of the Proof ---\n")
cat("Standard Threshold (c_alpha):", round(c_alpha, 3), "\n")
cat("Clinical Threshold (DV):     ", round(dv, 3), "\n")
cat("Theoretical Crossing Point (DV + c_alpha):", round(crossing_point, 3), "\n\n")

# Let's check a theta value BEFORE the crossing point (e.g., theta = 4)
cat("At theta = 4.0 (Before crossing point):\n")
cat("  DOR Standard:", round(calculate_dor(c_alpha, 4.0, se), 2), "\n")
cat("  DOR DC:      ", round(calculate_dor(dv, 4.0, se), 2), " <- HIGHER\n\n")

# Let's check a theta value AFTER the crossing point (e.g., theta = 6)
cat("At theta = 6.0 (After crossing point):\n")
cat("  DOR Standard:", round(calculate_dor(c_alpha, 6.0, se), 2), " <- HIGHER\n")
cat("  DOR DC:      ", round(calculate_dor(dv, 6.0, se), 2), "\n")

































