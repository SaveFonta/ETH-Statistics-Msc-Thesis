
## data <- read.table("table1.txt", header=TRUE)
## define matrices to store results
succMat <- futMat <- ExpNMat <- matrix(NA, ncol=5, nrow=5)
EUIIMat <- EUIIFutMat <- matrix(NA, ncol=5, nrow=4)


## Design 1
## GSD design with informative prior as in Gsponer et al 
design <- gsbDesign(nr.stages = 2, patients = c(10,20), sigma=mysigma,
                    criteria.success = c(0,0.95,50,0.5),
                    criteria.futility = c(40,0.9), prior.control = c(49,20),
                    prior.treatment = c(49,0))
## simulation scenarios
simulation <- gsbSimulation(truth=list(49,mydiff), grid.type="sliced",
                            type.update="per arm", nr.sim=100000,
                            warnings.sensitivity = 500, method="numerical integration")
## operating characteristics


## same results as in Gsponer paper, Table 1
result <- gsb(design=design,simulation=simulation)
## tables and graphics
tresult <- tab(result, what="cumulative all")
ssresult <- tab(result, what="sample size")
succ <- tresult[,7]
fut <- tresult[,8]
futMat[,1] <- fut
ExpN <- ssresult[,5]
succMat[,1] <- succ
ExpNMat[,1] <- ExpN

N0 <- ExpN[1]
N1 <- ExpN[-1]
SuccessOdds <- succ/(1-succ)
T1EOdds <- SuccessOdds[1]
PowerOdds <- SuccessOdds[-1]
EUII1 <- (PowerOdds^(1/N1))/(T1EOdds^(1/N0))
EUIIMat[,1] <- EUII1
SuccessFutOdds <- succ/fut
T1EOdds <- SuccessFutOdds[1]
PowerOdds <- SuccessFutOdds[-1]
EUIIFut1 <- (PowerOdds^(1/N1))/(T1EOdds^(1/N0))
EUIIFutMat[,1] <- EUIIFut1

## Design 2
## include 20 virtual patients from the prior in expected sample size
priorN <- 20
EUII0 <- (PowerOdds^(1/(priorN+N1)))/(T1EOdds^(1/(priorN+N0)))
succMat[,2] <- succ
futMat[,2] <- fut
ExpNMat[,2] <- ExpN + priorN
EUIIMat[,2] <- EUII0
SuccessFutOdds <- succ/fut
T1EOdds <- SuccessFutOdds[1]
PowerOdds <- SuccessFutOdds[-1]
EUIIFut0 <- (PowerOdds^(1/(priorN+N1)))/(T1EOdds^(1/(priorN+N0)))
EUIIFutMat[,2] <- EUIIFut0

## Design 3
## GSD design as before, but with uninformative prior
design <- gsbDesign(nr.stages = 2, patients = c(10,20), sigma=mysigma,
                    criteria.success = c(0,0.95,50,0.5),
                    criteria.futility = c(40,0.9), prior.control = c(49,0),
                    prior.treatment = c(49,0))
## simulation scenarios
simulation <- gsbSimulation(truth=list(49,mydiff), grid.type="sliced",
                            type.update="per arm", nr.sim=100000,
                            warnings.sensitivity = 500, method="numerical integration")
## operating characteristics
result <- gsb(design=design,simulation=simulation)
## tables and graphics
tresult <- tab(result, what="cumulative all")
ssresult <- tab(result, what="sample size")
succ <- tresult[,7]
fut <- tresult[,8]
futMat[,3] <- fut
ExpN <- ssresult[,5]
succMat[,3] <- succ
ExpNMat[,3] <- ExpN
SuccessOdds <- succ/(1-succ)
T1EOdds <- SuccessOdds[1]
PowerOdds <- SuccessOdds[-1]
N0 <- ExpN[1]
N1 <- ExpN[-1]
EUII2 <- (PowerOdds^(1/N1))/(T1EOdds^(1/N0))
EUIIMat[,3] <- EUII2
SuccessFutOdds <- succ/fut
T1EOdds <- SuccessFutOdds[1]
PowerOdds <- SuccessFutOdds[-1]
EUIIFut2 <- (PowerOdds^(1/N1))/(T1EOdds^(1/N0))
EUIIFutMat[,3] <- EUIIFut1

## Design 4
## fixed design with 40 patients per group
design <- gsbDesign(nr.stages = 1, patients = c(40,40), sigma=mysigma,
                    criteria.success = c(0,0.95,50,0.5),
                    criteria.futility = c(40,0.9), prior.control = c(49,0),
                    prior.treatment = c(49,0))
## simulation scenarios
simulation <- gsbSimulation(truth=list(49,mydiff), grid.type="sliced",
                            type.update="per arm", nr.sim=100000,
                            warnings.sensitivity = 500, method="numerical integration")
## operating characteristics
result <- gsb(design=design,simulation=simulation)
## tables and graphics
tresult <- tab(result, what="cumulative all")
ssresult <- tab(result, what="sample size")
succ <- tresult[,4]
fut <- tresult[,5]
futMat[,4] <- fut
succMat[,4] <- succ
N0 <- N1 <- 80
ExpNMat[,4] <- N0
T1EOdds <- succ[1]/(1-succ[1])
PowerOdds <- succ[-1]/(1-succ[-1])
EUIIFix <- (PowerOdds^(1/N1))/(T1EOdds^(1/N0))
EUIIMat[,4] <- EUIIFix
SuccessFutOdds <- succ/fut
T1EOdds <- SuccessFutOdds[1]
PowerOdds <- SuccessFutOdds[-1]
EUIIFutFix <- (PowerOdds^(1/N1))/(T1EOdds^(1/N0))
EUIIFutMat[,4] <- EUIIFutFix

## Design 5
## fixed design with 20 virtual prior patients and 20 / 40 patients per group
design <- gsbDesign(nr.stages = 1, patients = c(20,40), sigma=mysigma,
                    criteria.success = c(0,0.95,50,0.5),
                    criteria.futility = c(40,0.9), prior.control = c(49,20),
                    prior.treatment = c(49,0))
## simulation scenarios
simulation <- gsbSimulation(truth=list(49,mydiff), grid.type="sliced",
                            type.update="per arm", nr.sim=100000,
                            warnings.sensitivity = 500, method="numerical integration")
## operating characteristics
result <- gsb(design=design,simulation=simulation)
## tables and graphics
tresult <- tab(result, what="cumulative all")
ssresult <- tab(result, what="sample size")
succ <- tresult[,4]
fut <- tresult[,5]
futMat[,5] <- fut
N0 <- 80
N1 <- 80
succMat[,5] <- succ
ExpNMat[,5] <- N0

T1EOdds <- succ[1]/(1-succ[1])
PowerOdds <- succ[-1]/(1-succ[-1])
EUIIFix2 <- (PowerOdds^(1/N1))/(T1EOdds^(1/N0))
EUIIMat[,5] <- EUIIFix2
SuccessFutOdds <- succ/fut
T1EOdds <- SuccessFutOdds[1]
PowerOdds <- SuccessFutOdds[-1]
EUIIFutFix2 <- (PowerOdds^(1/N1))/(T1EOdds^(1/N0))
EUIIFutMat[,5] <- EUIIFutFix2

## results are identical to my functions for n=40 and sigma = sqrt(2)*mysigma (!)
succ2 <- c(DC.T1Erate(n=40, DV=50, sigma=mysigma*sqrt(2), alpha=0.025),
           DC.Power(n=40, theta=c(40,50,60,70), DV=50, sigma=mysigma*sqrt(2)))
## print(round(succ, 3))
## print(round(succ2, 3))
## and also to my functions for n=20 and sigma = mysigma (!)
succ3 <- c(DC.T1Erate(n=40, DV=50, sigma=mysigma*sqrt(2), alpha=0.025),
           DC.Power(n=40, theta=c(40,50,60,70), DV=50, sigma=mysigma*sqrt(2)))

## print(round(succMat[,c(2,3,5,4)], 3))
## print(round(futMat[,c(2,3,5,4)], 3))
## print(round(ExpNMat[,c(2,3,5,4)], 1))
## print(round(EUIIMat[,c(2,3,5,4)], 2))
## print(round(EUIIFutMat[,c(2,3,5,4)], 2))


EUIIInf <- DC.EUII(n=20000000,  theta=c(40,50,60,70), DV=50,
                     sigma=mysigma*sqrt(2))

EUIIFixLim <- exp(0.5*(((mydiff[-1]-50)/(mysigma*sqrt(2)))^2+
                        (50/(mysigma*sqrt(2)))^2)/2)
EUIIFixLim[1] <- NA

## print(EUIIInf)
## print(EUIIFixLim)


par(las=1)
d <- mydiff[-1]/mysigma
matplot(d, EUIIMat[,c(2,3,5,4)], type="b", lty=1, pch=19,
        xlab="Standardized effect size", ylab="EUII", log="y", ylim=c(1, max(EUIIMat[,-1])), lwd=2)
## matplot(d, cbind(EUII0, EUII2, EUIIFix2, EUIIFix), type="b", lty=1, pch=19,
##         xlab="Standardized effect size", ylab="EUII", log="y", ylim=c(1, max(EUIIMat[,-1])), lwd=2)
## for(i in c(1:4))
##     lines(d, EUIIFutMat[,c(2,3,5,4)][,i], lty=2, type="b", pch=19, lwd=2, col=i)
abline(h=1, lty=2, col="grey")
abline(v=50/mysigma, lty=2, col=1)
text(50/mysigma, max(EUII0), "standardized decision value", pos=4)
## legend("bottomright",
##        c("Informative BGSD not taking prior patients into account", "Informative BGSD taking prior patients into account", "Uninformative BGSD", "Bayesian fixed design"),
##        lty=1, col=c(1,2,3,4,5), bg="white", lwd=2, cex=0.85)
legend("bottomright",
       c("Informative BGS design", "Uninformative BGS design", "Informative Bayesian fixed design", "Uninformative Bayesian fixed design"),
       lty=1, col=c(1,2,3,4,5), bg="white", lwd=2, cex=0.85)

@ 
