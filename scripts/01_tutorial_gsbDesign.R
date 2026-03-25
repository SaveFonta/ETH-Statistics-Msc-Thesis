library(gsbDesign)




# Ma quindi assumiamo che il true value del placebo sia proprio quel 49 che abbiamo ottenuto come prior? 
# Pero facile cosi no? Significa che assumiamo con la prior di averci beccato nella mean

# NOTE: if the prior is provided per arm, you need simulations cannot use numerical integration




#FROM the paper on gsbDesing, they allow also the control to vary!! Bravo
desPoC5 <- gsbDesign(nr.stages = 2, patients = c(10, 20),
                         sigma = c(88, 88), criteria.success = c(0, 0.975, 50, 0.5),
                         criteria.futility = c(40, 0.9), prior.control = c(49, 20))
simPoC5 <- gsbSimulation(truth = cbind(rep(c(30, 50, 70), each = 5),
                                            c(30, 70, 80, 90, 100, 50, 90, 100, 110, 120, 70, 110, 120, 130, 140)),
                              nr.sim = 20000, type.update = "per arm", method = "simulation",
                         grid.type = "manually")
ocPoC5 <- gsb(desPoC5, simPoC5)





#NOTE 
#FROM THE HELP
#If type.update = "treatment effect", the Bayesian update from prior to posterior is calculated on 
#treatment effect delta. If type.update = "per arm", the update is calculated separately in the 
#placebo and the treatment arm. In this case it is possible to enter prior information in only one
#arm. For type.update = "per arm" only a simulation method is implemented.


# --> THIS doest mean that only one of the two priors can be define, but both are allowed, but also only one is allowed





####################

# Let' s do some practice
mio_design <- gsbDesign(nr.stages = 2, patients = c(50,60), sigma = c(88,88), criteria.success = c(0, 0.95, 50, 0.5), criteria.futility = c(40, 0.9), prior.control = c(49, 10), prior.treatment = c(52, 10))

#Curiouly --> prior.control and prior.treatment arre defined not as normals, but as c(mean, # patients), so  prior.control = c(0,2) means an informative prior equivalent to 2 patients in control arm with mean = 0 and I GUESS the same variance sigma.



# gsbSimualation is a bit more complicated. 
# Brief Eplanation 
# Core params are: 
# 1) type_update: tells how to perform the update 
#{"treatment effect"} means the bayesian update is computed on the treatment difference (not our case usually)
# if {"per arm"}: The update is calculated separately for the placebo and the treatment arms.

# 2) method: tells how to compute the probabilities 
# {"numerical integration"}: Uses exact mathematical integration (fastest and default when type.update = "treatment effect").

# {"simulation"}: Generates thousands of virtual trials. This is mandatory if you use type.update = "per arm".

#{"both"}

# Themn we have truth parameter, which changes depending on type.update and grid.type 

# Scenario A: type.update = "treatment effect"
# Don't need grid.type. THen truth myst be a vector of length 3:
# Element 1: The minimum true treatment difference to evaluate (e.g., -10).
#Element 2: The maximum true treatment difference to evaluate (e.g., 60).
#Element 3: The number of distinct values/steps to evaluate between the minimum and maximum

#otwise can even use grid.type = "manually" and giving truth = numeric containing true delta values 


#Scenario B: type.update = "per arm".
#Here need grid.type to tell the package how to interpret truth parameter 

# grid.type = "table": truth must be a list with two elements. The first element is a sequence of true values for the control arm, 
# and the second is a sequence for the treatment arm. The package evaluates every possible combination of these two sequences.

# grid.type = "sliced": truth must be a list with two elements. The first element is a vector of true control values, and the 
# second element is a vector of true deltas. Useful when you want to fix a specific control rate (like a historical placebo rate of 49) 
# and evaluate a specific sequence of relative treatment differences (e.g., 0, 40, 50, 60, 70).


# If grid.type = "manually": truth must be a matrix nx2. Each row of the matrix corresponds to a specific, hand-picked 
# pair of true values for the control and treatment arms.

# If grid.type = "plot": truth must be a vector of length 5. The first two numbers define the range for the control arm, the next two 
#define the range for the treatment arm, and the fifth number indicates the number of grid points to use to generate a contour plot.
 # truth = c(min.placebo, max.placebo,min.treatment, max.treatment, n), here n is the number of grid points.

mia_simulation <- gsbSimulation(type.update = "per arm", method = "simulation", grid.type = "table", truth = list (seq(from = 20, to = 60, by = 10), seq(30, 70, 10)))





# LAST FUNCTION: gsb

result <- gsb (simulation = mia_simulation, design = mio_design)

#two main elemtns:
result$delta.grid

#to see all the evaluated delta used


#best way to see at results:
table_ss <- tab(result, what = "sample size")
# max sample size is (50 + 60) + (50 + 60)
View(table_ss)

tab(result, what="cumulative all")

# you can specify:
# "what",  choices:

#"success", "futility", "indeterminate": Shows the probability of stopping at exactly that specific stage.
#"cumulative success", "cumulative futility", "cumulative all": Shows the running total probability of stopping by that stage or earlier.
#"sample size": Outputs the expected number of patients for each scenario.


# Note: This handy trick only works if your design used type.update = "treatment effect" 
#By default (atDelta = "default"), the table will only show the operating characteristics for the exact "true" effect sizes you specified in the truth argument of gsbSimulation().
#However, if you want to know the probability of success for a true effect size that you didn't explicitly simulate, you can provide a numeric vector here (e.g., atDelta = c(15, 25)). The function will automatically calculate the values for those specific deltas using linear interpolation between your closest simulated points.

#then can specify wide and digits, can even easily export

# ==============================================================================
# INSPECTION OF of the 'tab()' BUG
# ==============================================================================

cumall <- tab(result_new, what = "cumulative all")[1,] # The MATHEMATICALLY CORRECT table
normal <- tab(result_new)[1,]                          # The DEFAULT table (Contains the bug!)

# Specific extractions of marginal (absolute) probabilities
# Note: "success" and "futility" extract correct absolute probabilities.
# "indeterminate" extracts the buggy remainder for Stage 2.
tab(result_new, what = "success")[1,]
tab(result_new, what = "futility")[1,]
tab(result_new, what = "indeterminate")[1,] 

# ------------------------------------------------------------------------------
# 2. THE BUG IN THE DEFAULT TABLE (normal)
# ------------------------------------------------------------------------------
# At first glance, Stage 2 appears to be normalized (conditional probabilities) 
# because it sums to 1. 

sum(normal[c("stage1.suc", "stage1.fut", "stage1.ind")]) # Sums to 1 (Correct)
sum(normal[c("stage2.suc", "stage2.fut", "stage2.ind")]) # Sums to 1 (THE ILLUSION)

# THE TRUTH: stage2.suc and stage2.fut are NOT conditional probabilities. 
# They are ABSOLUTE marginal probabilities (out of the original 100% of trials).
# The package forces Stage 2 to sum to 1 by using a flawed formula for stage2.ind:
# (1.000 - stage2.suc - stage2.fut), ignoring that 64% of trials stopped in Stage 1!

# ------------------------------------------------------------------------------
# 3. WHY THE CONDITIONAL HYPOTHESIS FAILS
# ------------------------------------------------------------------------------
# Because stage2.suc is ALREADY an absolute probability (out of 100% of trials),
# you DO NOT multiply it by the trials that made it to Stage 2 (stage1.ind).
# Doing so applies a conditional penalty to a number that is already absolute.

# WRONG (Conditional logic applied to absolute numbers):
normal["stage1.suc"] + (normal["stage1.ind"] * normal["stage2.suc"]) # Does not match cumall

# ------------------------------------------------------------------------------
# 4. THE CORRECT MATH (ABSOLUTE ADDITION) vs. THE BUG
# ------------------------------------------------------------------------------
# Since success and futility are absorbing states, and the probabilities are absolute,
# you find the cumulative total by simply adding them directly.

# CORRECT: Cumulative Success
normal["stage1.suc"] + normal["stage2.suc"] # Matches cumall["stage2.suc"] exactly

# CORRECT: Cumulative Futility
normal["stage1.fut"] + normal["stage2.fut"] # Matches cumall["stage2.fut"] exactly

# THE BUG REVEALED: Indeterminate Trials
# "Indeterminate" trials are a transient state; they DO NOT accumulate. 
# The true number of trials finishing inconclusive is the trials that ENTERED Stage 2 
# minus the trials that STOPPED in Stage 2. Let's calculate the true remainder:

true_stage2_ind <- normal["stage1.ind"] - normal["stage2.suc"] - normal["stage2.fut"]

# Compare our true math with the package's outputs:
true_stage2_ind                   # Correct remainder (0.145)
cumall["stage2.ind"]              # Correct remainder (0.145) -> 'cumulative all' is right!
normal["stage2.ind"]              # BUGGY remainder   (0.785) -> Default table is wrong!

# ==============================================================================
# FINAL ADVICE & BEST PRACTICES FOR THESIS
# ==============================================================================
# 1. IGNORE THE DEFAULT TABLE: Never use tab(result) or tab(result, what="all") or 
#    tab(result, what="indeterminate"). The 'indeterminate' columns for any stage > 1 
#    are mathematically impossible and create an optical illusion.
# 2. USE CUMULATIVE: ALWAYS use tab(result, what="cumulative all") to see the 
#    true, mathematically sound progression of the simulated trial.






# NOW we have PLOT
plot(result, what="sample size")

plot(result, what="cumulative all")

#these graphs are pretty rettible to look, we should use the grid.type = "sliced" to have nicer plots. 



#fast esperimetno con slide (stesso di MING)

design <- gsbDesign(nr.stages = 2,
                    patients = c(10,20),
                    sigma=88,
                    criteria.success = c(0,0.95,50,0.5),
                    criteria.futility = c(40,0.9),
                    prior.control = c(49,20),
                    prior.treatment = c(49,0))


simulation <- gsbSimulation(
  truth=list(49,c(0,40,50,60,70)),
  grid.type="sliced",
  method = "both",
  type.update="per arm",
  nr.sim=10000,
  warnings.sensitivity = 500,
  seed=1)

result_new <- gsb(design=design,simulation=simulation)


plot(result_new, what="sample size", sliced=TRUE)
plot(result_new, what="cumulative all", sliced=TRUE)

#AHHH beautiful plots 4real 



