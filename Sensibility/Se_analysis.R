# Se sensitivity analysis using PRCC ----
# Matemathical modelling of inmunization strategies

# Load function simulationSIRS
setwd("~/Dropbox/0.USP/7.Publicações/Mathematical modelling CSFV/sensibility/")
source("SIVRP_se.R")

# Se analysis using PRCC (Partial Rank Correlation Coefficient) ----
# Suposing a uniform correlation
library(zoo)
library(deSolve)
library(ggplot2)
library(lubridate)
library(dplyr)
library(epiR)

set.seed(555)

start_time_t <- Sys.time()

# 1 beta ----
start_time <- Sys.time()
beta <- runif(n = 100, min = 4.0706 * 0.85, max = 4.0706 * 1.15)
# Parameters, comented the one that ----

# beta <- 4.0706 #Calculada nls
gama <- 2.7941 # Calculada nls
tau <- gama / 30 #
mu <- (1 / (350 / 30))
w <- 1 / (180 / 30) # omega
theta <- 0.1

# Definindo o vetor para a prevalencia media da simulacao
prevalence <- numeric(length(beta))

for (i in 1:length(beta)) {
  mod <- SIRVP_se(beta[i], gama = gama, tau = tau, mu = mu, w = w, theta = theta, tf = tf)
  prevalence[i] <- sum(mod$p)
}

# Analise de sensibilidade utilizando PRCC para Infh
## Alfredo: o pacote epi.prcc entende que Y e os valores a calcular
dat_1 <- data.frame(cbind(X1 = beta, Y = prevalence))
prcc_beta <- epi.prcc(dat_1, sided.test = 2)
plot(dat_1$X1, dat_1$Y, xlab = "beta", ylab = "Mean of Infected", cex = .5)

end_time <- Sys.time()
end_time - start_time

# 2 gama ----
start_time <- Sys.time()
gama <- runif(n = 100, min = 2.7941 * 0.85, max = 2.7941 * 1.15)

# Parameters, comented the one that is tested----

beta <- 4.0706 # Calculada nls
# gama <- 2.7941 #Calculada nls
tau <- 2.7941 / 30 #
mu <- (1 / (350 / 30))
w <- 1 / (180 / 30) # omega
theta <- 0.1

# Definindo o vetor para a prevalencia media da simulacao
prevalence <- numeric(length(gama))

for (i in 1:length(gama)) {
  mod <- SIRVP_se(beta = beta, gama[i], tau = tau, mu = mu, w = w, theta = theta, tf = tf)
  prevalence[i] <- sum(mod$p)
}

# Sensibility analysis using PRCC para Infh
## Alfredo: o pacote epi.prcc entende que Y e os valores a calcular
dat_2 <- data.frame(cbind(X1 = gama, Y = prevalence))
prcc_gama <- epi.prcc(dat_2, sided.test = 2)
plot(dat_2$X1, dat_2$Y)
plot(dat_2$X1, dat_2$Y, xlab = "gama", ylab = "Mean of Infected", cex = .5)

# Time
end_time <- Sys.time()
start_time - end_time

# 3 tau ----
start_time <- Sys.time()
tau <- runif(n = 100, min = gama / 30 * 0.85, max = gama / 30 * 1.15)

# Parameters, comented the one that is tested----

beta <- 4.0706 # Calculada nls
gama <- 2.7941 # Calculada nls
# tau <- gama/30 #
mu <- (1 / (350 / 30))
w <- 1 / (180 / 30) # omega
theta <- 0.1

# Definindo o vetor para a prevalencia media da simulacao
prevalence <- numeric(length(tau))

for (i in 1:length(tau)) {
  mod <- SIRVP_se(beta = beta, gama = gama, tau[i], mu = mu, w = w, theta = theta, tf = tf)
  prevalence[i] <- sum(mod$p)
}

# Sensibility analysis using PRCC para Infh
## Alfredo: o pacote epi.prcc entende que Y e os valores a calcular
dat_3 <- data.frame(cbind(X1 = tau, Y = prevalence))
prcc_tau <- epi.prcc(dat_3, sided.test = 2)
plot(dat_2$X1, dat_2$Y)
plot(dat_2$X1, dat_2$Y, xlab = "tau", ylab = "Mean of Infected", cex = .5)

# Time
end_time <- Sys.time()
start_time - end_time

# 4 mu ----
start_time <- Sys.time()
mu <- runif(n = 100, min = (1 / (350 / 30)) * 0.85, max = (1 / (350 / 30)) * 1.15)

# Parameters, comented the one that is tested----

beta <- 4.0706 # Calculada nls
gama <- 2.7941 # Calculada nls
tau <- gama / 30 #
# mu <- (1/(350/30))
w <- 1 / (180 / 30) # omega
theta <- 0.1

# Definindo o vetor para a prevalencia media da simulacao
prevalence <- numeric(length(mu))

for (i in 1:length(mu)) {
  mod <- SIRVP_se(beta = beta, gama = gama, tau = tau, mu[i], w = w, theta = theta, tf = tf)
  prevalence[i] <- sum(mod$p)
}

# Sensibility analysis using PRCC para Infh
## Alfredo: o pacote epi.prcc entende que Y e os valores a calcular
dat_4 <- data.frame(cbind(X1 = mu, Y = prevalence))
prcc_mu <- epi.prcc(dat_4, sided.test = 2)
plot(dat_2$X1, dat_2$Y)
plot(dat_2$X1, dat_2$Y, xlab = "mu", ylab = "Mean of Infected", cex = .5)

# Time
end_time <- Sys.time()
start_time - end_time

# 5 w ----
start_time <- Sys.time()
w <- runif(n = 100, min = 1 / (180 / 30) * 0.85, max = 1 / (180 / 30) * 1.15)

# Parameters, comented the one that is tested----

beta <- 4.0706 # Calculada nls
gama <- 2.7941 # Calculada nls
tau <- gama / 30 #
mu <- (1 / (350 / 30))
# w <- 1/(180/30) #omega
theta <- 0.1

# Definindo o vetor para a prevalencia media da simulacao
prevalence <- numeric(length(w))

for (i in 1:length(w)) {
  mod <- SIRVP_se(beta = beta, gama = gama, tau = tau, mu = mu, w[i], theta = theta, tf = tf)
  prevalence[i] <- sum(mod$p)
}

# Sensibility analysis using PRCC para Infh
## Alfredo: o pacote epi.prcc entende que Y e os valores a calcular
dat_5 <- data.frame(cbind(X1 = w, Y = prevalence))
prcc_w <- epi.prcc(dat_5, sided.test = 2)
plot(dat_2$X1, dat_2$Y, xlab = "w", ylab = "Mean of Infected", cex = .5)

# Time
end_time <- Sys.time()
end_time - start_time

# 6 theta ----
start_time <- Sys.time()
theta <- runif(n = 100, min = 0.1 * 0.85, max = 0.1 * 1.15)

# Parameters, comented the one that is tested----
beta <- 4.0706 # Calculada nls
gama <- 2.7941 # Calculada nls
tau <- gama / 30 #
mu <- (1 / (350 / 30))
w <- 1 / (180 / 30) # omega
# theta <- 0.1

# Definindo o vetor para a prevalencia media da simulacao
prevalence <- numeric(length(theta))

for (i in 1:length(theta)) {
  mod <- SIRVP_se(beta = beta, gama = gama, tau = tau, mu = mu, w = w, theta[i], tf = tf)
  prevalence[i] <- sum(mod$p)
}

# 7 Sensibility analysis using PRCC para Infh ----
## Alfredo: o pacote epi.prcc entende que Y e os valores a calcular
dat_6 <- data.frame(cbind(X1 = theta, Y = prevalence))
prcc_theta <- epi.prcc(dat_6, sided.test = 2)
plot(dat_2$X1, dat_2$Y, xlab = "theta", ylab = "Mean of Infected", cex = .5)

prcc_beta
prcc_gama
prcc_tau
prcc_mu
prcc_w
prcc_theta

# Time
end_time <- Sys.time()
end_time - start_time

# Graphic
setwd("/home/alfredo/Dropbox/0.USP/7.Publicações/Mathematical modelling CSFV/Graphics python/sensibility/")

par(mfrow = c(2, 3))

plot(dat_1$X1, dat_1$Y, xlab = "beta", ylab = "Cummulative cases", cex = .1)
plot(dat_2$X1, dat_2$Y, xlab = "gama", ylab = "Cummulative cases", cex = .1)
plot(dat_3$X1, dat_3$Y, xlab = "tau", ylab = "Cummulative cases", cex = .1)
plot(dat_4$X1, dat_4$Y, xlab = "mu", ylab = "Cummulative cases", cex = .1)
plot(dat_5$X1, dat_5$Y, xlab = "w", ylab = "Cummulative cases", cex = .1)
plot(dat_6$X1, dat_6$Y, xlab = "theta", ylab = "Cummulative cases", cex = .1)

par(mfrow = c(1, 1))

# Total Time
end_time_t <- Sys.time()
end_time_t - start_time_t


pr1 <- data.frame("beta", c(prcc_beta[c(2:7)]))
pr2 <- data.frame("gama", c(prcc_gama[c(2:7)]))
pr3 <- data.frame("tau", c(prcc_tau[c(2:7)]))
pr4 <- data.frame("mu", c(prcc_mu[c(2:7)]))
pr5 <- data.frame("w", c(prcc_w[c(2:7)]))
pr6 <- data.frame("theta", c(prcc_theta[c(2:7)]))
colnames(pr1) <- c("Parameter", "est", "Lower-Bound", "Upper-Bound", "test-stat", "df", "p.value")
colnames(pr2) <- c("Parameter", "est", "Lower-Bound", "Upper-Bound", "test-stat", "df", "p.value")
colnames(pr3) <- c("Parameter", "est", "Lower-Bound", "Upper-Bound", "test-stat", "df", "p.value")
colnames(pr4) <- c("Parameter", "est", "Lower-Bound", "Upper-Bound", "test-stat", "df", "p.value")
colnames(pr5) <- c("Parameter", "est", "Lower-Bound", "Upper-Bound", "test-stat", "df", "p.value")
colnames(pr6) <- c("Parameter", "est", "Lower-Bound", "Upper-Bound", "test-stat", "df", "p.value")

df <- rbind(pr1, pr2, pr3, pr4, pr5, pr6)
df$UL_Difference <- df[4] - df[5]

# Graficando o PRCC ----
# tor <-

ggplot(df) +
  geom_col(aes(x = reorder(Parameter, est), y = est), fill = "gray60") +
  geom_text(aes(x=reorder(Parameter, est), y=0.9*est, label=round(est,5))) +
  theme_minimal() +
  xlab(NULL) +
  ylab("PRCC") +
  theme(text = element_text(size = 14))

# scale_y_continuous(trans='sqrt')

# tiff(
#   filename = "Tornado_Plot_1000_mumulative_infected_15_90%_coverage%.tiff",
#   width = 9, height = 6, units = "cm",
#   compression = "lzw", pointsize = 12, res = 600
# )
# tor
# dev.off()
#     