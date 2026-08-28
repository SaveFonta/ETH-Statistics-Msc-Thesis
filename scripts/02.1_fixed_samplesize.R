## ---------------------------------------------------------------------------
## Required control sample size for each analysis prior.
##
## The closed form of the sample size applies only to a single Normal analysis
## prior. So we can use numerical inversion: avgoc2S.normMix evaluates avgPower
## for any Normal mixture, and avgPower is monotone in n_C, so bisection finds
## the smallest adequate n_C.
##
## The design prior is fixed at pi_MAP for every analysis prior. pi_MAP is a
## statement about theta_C derived from the historical data; the robust
## component of a robustified pi_A is a device for discounting those data, not
## a claim about where theta_C lies, so it is not used as the design prior.
## ---------------------------------------------------------------------------

suppressMessages(library(RBesT))
library(parallel)
source("scripts/00_functions.R")
source("scripts/00_shared_setup.R")

TARGET <- 0.80     # target average power
DELTA  <- 60       # treatment effect at which power is evaluated
rT     <- 2        # allocation ratio, n_T = rT * n_C
OUT    <- "Output/02.1_fixed_samplesize"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## the single Normal used for the closed form results of the methodology
p_norm1 <- mixnorm(c(1, mu_A, sigma / sqrt(20)), sigma = sigma, param = "ms")
priors  <- c(list("1 component Normal" = p_norm1), analysis_priors)

avgpow <- function(pa, pd, n_C) {
  f <- avgoc2S.normMix(prior1 = p_vague, prior2 = pa, n1 = rT * n_C, n2 = n_C,
                       decision = sign.crit, delta = 0, design_prior2 = pd)
  as.numeric(f(-DELTA, pd))
}

avgt1e <- function(pa, pd, n_C) {
  f <- avgoc2S.normMix(prior1 = p_vague, prior2 = pa, n1 = rT * n_C, n2 = n_C,
                       decision = sign.crit, delta = 0, design_prior2 = pd)
  as.numeric(f(0, pd))
}

smallest_n <- function(pa, pd, lo = 2, hi = 400) {
  if (avgpow(pa, pd, hi) < TARGET) return(NA_integer_) #at 400 almost impossible
  while (lo < hi) {
    mid <- (lo + hi) %/% 2
    if (avgpow(pa, pd, mid) >= TARGET) hi <- mid else lo <- mid + 1
  }
  lo
}

res <- mclapply(names(priors), function(nm) {
  pa <- priors[[nm]]
  n  <- smallest_n(pa, p_MAP)
  e  <- as.numeric(ess(pa))
  data.frame(prior = nm, ess = e, n_C = n, n_T = rT * n,
             n_total = (1 + rT) * n, avg_power = avgpow(pa, p_MAP, n),
             avg_t1e = avgt1e(pa, p_MAP, n))
}, mc.cores = length(priors))
res <- do.call(rbind, res)

## closed form reference for the single Normal, Equation (n_c_formula)
K  <- sigma^2 * (qnorm(0.95) + qnorm(TARGET))^2 / DELTA^2
nA <- 20
n_closed <- (K * (rT + 1) - rT * nA +
             sqrt((K * (rT + 1) - rT * nA)^2 + 4 * rT * K * nA)) / (2 * rT)

saveRDS(list(res = res, n_closed = n_closed, target = TARGET,
             delta = DELTA, rT = rT),
        file.path(OUT, "samplesize.RDS"))
