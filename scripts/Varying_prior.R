library(gsbDesign)

design <- gsbDesign(nr.stages = 2, 
                        patients = c(10, 20), 
                        sigma = c(88, 88), 
                        criteria.success = c(0, 0.95, 50, 0.5), 
                        criteria.futility = c(40, 0.9), 
                        prior.control = c(49, 20), 
                        prior.treatment = c(49, 0))


sim <- gsbSimulation(truth = cbind(rep(c(30, 50, 70), each = 5),  
                      c(30, 70, 80, 90, 100, 50, 90, 100, 110, 120, 70, 110, 120, 130, 140)), 
                      nr.sim = 100000, type.update = "per arm", method = "simulation", grid.type = "manually")

res <- gsb(design, sim)

tab(res, what = "cumulative all")


library(ggplot2)
df <- as.data.frame(tab(res, what = "cumulative success", export = FALSE))

ggplot(df, aes(x = delta, y = `stage2`, color = factor(control))) +
  geom_line(size = 1) + 
  geom_point(size = 3) + 
  labs(y = "Probability of Success", color = "True Control")

#dev.off()









#########################################################################################################


# Try to prove it is better not to include the virtual sample size with an extreme case
design <- gsbDesign(nr.stages = 1,
                    patients = c(1,1),
                    sigma=88,
                    criteria.success = c(0,0.95,50,0.5),
                    criteria.futility = c(40,0.9),
                    prior.control = c(49,200), #
                    prior.treatment = c(49,0))


simulation <- gsbSimulation(
  truth=list(49,c(0,40,50,60,70)),
  grid.type="sliced",
  method = "both",
  type.update="per arm",
  nr.sim=100000, 
  warnings.sensitivity = 500,
  seed=1)


result_new <- gsb(design=design,simulation=simulation)



table <- tab(result_new, what = "cumulative all")

T1E <- table[1, "stage1.suc"]
Power <- table[-1, "stage1.suc"]
LR_pos <- Power / T1E
LR_neg <- (1 - Power) / (1 - T1E)
DOR <- LR_pos / LR_neg

EUII_virtual <- DOR^(1/202) #if we use the total sample size 
EUII_virtual 

EUII_pis <- DOR^(1/2) #if we use the total sample size 
EUII_pis 





#build the same, no prior
design <- gsbDesign(nr.stages = 1,
                    patients = c(1,1),
                    sigma=88,
                    criteria.success = c(0,0.95,50,0.5),
                    criteria.futility = c(40,0.9),
                    prior.control = c(49,0), #
                    prior.treatment = c(49,0))


simulation <- gsbSimulation(
  truth=list(49,c(0,40,50,60,70)),
  grid.type="sliced",
  method = "both",
  type.update="per arm",
  nr.sim=100000, 
  warnings.sensitivity = 500,
  seed=1)


result_new <- gsb(design=design,simulation=simulation)
# 


table <- tab(result_new, what = "cumulative all")

T1E <- table[1, "stage1.suc"]
Power <- table[-1, "stage1.suc"]
LR_pos <- Power / T1E
LR_neg <- (1 - Power) / (1 - T1E)
DOR <- LR_pos / LR_neg


EUII <- DOR^(1/2) #if we use the total sample size 
EUII


