# Sensitivity analysis using PRCC ----
# Mathematical modelling of immunization strategies

# Load function simulation SIRS
setwd("~/papers/modelling_csf/Modelling-imnunization-csf-ecuador/Sensitivity/")
source("SIVRP_se.R")

# Se analysis using PRCC (Partial Rank Correlation Coefficient) ----
# Supposing an uniform correlation
library(zoo)
library(deSolve)
library(ggplot2)
library(lubridate)
library(dplyr)
library(epiR)

set.seed(555)

start_time_t <- Sys.time()

# 1 All together ----
start_time <- Sys.time()
beta <- runif(n = 1000, min = 1.1664 * 0.85, max = 1.1664 * 1.15)
betab <- runif(n = 1000, min = 4.0706 * 0.38, max = 4.0706 * 0.42)
gama <- runif(n = 1000, min = 0.6937 * 0.85, max = 0.6937 * 1.15)
mu <- runif(n = 1000, min = (1 / (370 / 30)), max = (1 / (330 / 30)))
w <- runif(n = 1000, min =  (1 / (200 / 30)), max = (1 / (160 / 30)))
theta <- runif(n = 1000, min = 0.35 * 0.85, max = 0.35 * 1.15)
tau <- runif(n = 1000, min = 2.7941 / 33, max = 2.7941 / 25)
d <- runif(n = 1000, min = (1/(12/30)), max = (1/(8/30)))
r <- runif(n = 1000, min = (1/(17/30)), max = (1/(13/30)))

# Parameters, commented the one that ----
# beta <- 1.1658530 #best
# # beta <- 4.0706 # Calculated nls monthly
# betab <- beta * 0.4 # Assumption
# gama <- 0.6937 #best
# beta <- 4.0706 # Calculated nls monthly
# betab <- 4.0706 * 0.4 # Assumption
# gama <- 2.7941 # Calculated nls monthly                
# mu <- (1 / (350 / 30)) # natural death rate
# w <- 1 / (180 / 30) # omega loss of vaccine induced immunity
# theta <- 0.35 # percentage of animals persistent infected
# tau <- gama / 30 # death rate of persistent infected
# d <- (1/(10/30))  # Monthly number of days to natural decay of carcass to be used as swill
# r <- (1/(15/30))  # Monthly number of days to remove the carcass

# Defining the prevalence as the output of the simulation ----
prevalence <- numeric(length(beta))

for (i in 1:length(beta)) {
  mod <- SIRVP_se(beta[i],
                  betab[i],
                  gama[i],
                  mu[i],
                  w[i],
                  theta[i],
                  tau[i],
                  d[i],
                  r[i],
                  tf = tf)
  prevalence[i] <- sum(mod$In)
}

# Sensibility analysis using PRCC for prevalence ----

dat_1 <- data.frame(cbind(beta,
                          betab,
                          gama,
                          mu,
                          w,
                          theta,
                          tau,
                          d,
                          r,
                          Y = prevalence))

prcc_beta <- epi.prcc(dat_1, sided.test = 2)
prcc_beta


# Graphic
par(mfrow = c(3, 3))
colnames(dat_1)
plot(dat_1$beta, dat_1$Y, xlab = "beta", ylab = "Prevalence", cex = .1)
plot(dat_1$betab, dat_1$Y, xlab = "betab", ylab = "Prevalence", cex = .1)
plot(dat_1$gama, dat_1$Y, xlab = "gamma", ylab = "Prevalence", cex = .1)
plot(dat_1$mu, dat_1$Y, xlab = "mu", ylab = "Prevalence", cex = .1)
plot(dat_1$w, dat_1$Y, xlab = "w", ylab = "Prevalence", cex = .1)
plot(dat_1$theta, dat_1$Y, xlab = "theta", ylab = "Prevalence", cex = .1)
plot(dat_1$tau, dat_1$Y, xlab = "tau", ylab = "Prevalence", cex = .1)
plot(dat_1$d, dat_1$Y, xlab = "d", ylab = "Prevalence", cex = .1)
plot(dat_1$r, dat_1$Y, xlab = "r", ylab = "Prevalence", cex = .1)

par(mfrow = c(1, 1))

write.csv(prcc_beta, file = "prcc.csv")

prcc <- read.csv(file = "prcc.csv")

# Graphic PRCC ----
ggplot(prcc_beta) +
  geom_col(aes(x = reorder(var, est), y = est), fill = "gray60") +
  geom_text(aes(x=reorder(var, est), y=1.1*est, label=round(est,3))) +
  # geom_text(aes(x=reorder(var, est), y=-0.01, label=var)) +
  theme_minimal() +
  xlab(NULL) +
  ylab("PRCC") +
  theme(text = element_text(size = 14))

end_time <- Sys.time()
end_time - start_time
