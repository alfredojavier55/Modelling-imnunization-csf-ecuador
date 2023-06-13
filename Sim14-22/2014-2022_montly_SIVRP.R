###### 2014-2022 S-I-PI-V Simulation ######----

library(deSolve)
library(ggplot2)
library(lubridate)
library(dplyr)
library(zoo)
library(readr)

# Diretorio de trabalho
setwd("~/Dropbox/0.USP/7.Publicações/Mathematical modelling CSFV/Analise/")

# Vaccine coverage by months
coverage <- read.csv(file = "cobertura_vac_2017_21.csv")
coverage <- coverage[, -1]

# 12 months # 5 years (2016-2022) # 12*7 # 84 months
# Insert new lacking months (june and july 2016)
coverage <- rbind(coverage[1:5, ], c("2016-06-01", 0, "0.0"), coverage[-(1:5), ])
coverage <- rbind(coverage[1:6, ], c("2016-07-01", 0, "0.0"), coverage[-(1:6), ])
str(coverage)

# Creating lacking years 2014-2015
mes <- as.character(seq(as.Date("2014-01-01"), as.Date("2015-12-30"), by = "months"))
por <- c(rep(0.27, 12), rep(0.53, 12))
vacunados <- por * 278321
coverage_14_15 <- cbind(mes, vacunados, por)
str(coverage_14_15)

# joining the two
coverage <- rbind(coverage_14_15, coverage)
coverage$vacunados <- as.numeric(coverage$vacunados)

# Replacing from 2016-01 to 2016-11 53% of coverage p 278321
s <- coverage$vacunados[25:36] + 147510.13

coverage$vacunados[25:36] <- s

# Calculating the roll mean k=3 ----
# There are 82 moving means I recalculate the ones that are missing considering the last data
moving <- rollmean(coverage$vacunados, k = 3, align = c("right"))
moving2 <- coverage$vacunados[c(105:108)]
moving3 <- rollmean(moving2, k = 3)
moving <- c(moving, moving3)

# 95% of coverage is the max historic data
coverage$mov_mean <- round(moving, 2)
coverage$vac_cov <- round(coverage$mov_mean / (278321 * 1.00), 3)

plot(ymd(coverage$mes), coverage$vacunados, type = "l", col = "grey")
lines(ymd(coverage$mes), moving, type = "l", col = "red")
plot(ymd(coverage$mes), coverage$vac_c, type = "l", col = "red")

# Montly means
coverage %>%
  group_by(year(ymd(mes))) %>%
  summarize(
    x_vaccinated_pigs = mean(vacunados),
    x_coverage = mean(vac_cov),
    x_mobile_mean = mean(mov_mean)
  )

# Graphic historic coverage
coverage %>%
  group_by(
    month_l = month(ymd(mes)),
    year = factor(year(mes))
  ) %>%
  filter(year == 2021 | year == 2020 | year == 2019 |
    year == 2018 | year == 2017) %>%
  summarize(
    x_vaccinated_pigs = vacunados,
    x_coverage = vac_cov,
    x_mobile_mean = mean(mov_mean)
  ) %>%
  ggplot(aes(month_l, x_coverage, color = year)) +
  geom_line(size = 1.5, alpha = 0.7) +
  xlab("") +
  xlim(1, 12) +
  ylab("Vaccine coverage (%)") +
  theme_minimal() +
  scale_x_discrete(
    limit = c(1:12),
    labels = c(month(1:12, label = T))
  )

# Graphic historic number of animals
coverage %>%
  group_by(month_l = month(ymd(mes)), year = factor(year(mes))) %>%
  filter(year == 2021 | year == 2020 | year == 2019 |
    year == 2018 | year == 2017) %>%
  dplyr::summarize(
    x_vaccinated_pigs = vacunados,
    x_coverage = vac_cov,
    x_mobile_mean = mean(mov_mean)
  ) %>%
  ggplot(aes(month_l, x_vaccinated_pigs, color = year)) +
  geom_line(size = 1.5, alpha = 0.7) +
  xlab("") +
  xlim(1, 12) +
  ylab("Number of applied doses") +
  theme_minimal() +
  scale_x_discrete(
    limit = c(1:12),
    labels = c(month(1:12, label = T))
  ) +
  guides(colour = guide_legend(reverse = T))

# Graphic historic number of animals
coverage %>%
  group_by(month_l = month(ymd(mes)), year = factor(year(mes))) %>%
  filter(year == 2021 | year == 2020 | year == 2019 |
    year == 2018 | year == 2017) %>%
  dplyr::summarize(
    x_vaccinated_pigs = vacunados,
    x_coverage = vac_cov,
    x_mobile_mean = mean(mov_mean)
  ) %>%
  ggplot(aes(month_l, x_vaccinated_pigs / 1.55, color = year)) +
  geom_line(size = 1.5, alpha = 0.7) +
  xlab("") +
  xlim(1, 12) +
  ylab("Number of applied doses") +
  theme_minimal() +
  scale_x_discrete(
    limit = c(1:12),
    labels = c(month(1:12, label = TRUE))
  ) +
  guides(colour = guide_legend(reverse = TRUE))

# Yearly means
coverage %>%
  group_by(year(ymd(mes))) %>%
  summarize(
    x_vaccinated_pigs = sum(vacunados),
    x_coverage = mean(vac_cov),
    x_mobile_mean = sum(mov_mean)
  )

# Graphics
coverage %>%
  group_by(m = floor_date(ymd(mes), unit = "month")) %>%
  filter(m > "2017-01-01") %>%
  dplyr::summarize(
    Applied_doses = vacunados,
    x_coverage = mean(vac_cov),
    Mobile_mean = mean(mov_mean)
  ) %>%
  ggplot() +
  geom_line(aes(m, Applied_doses), size = 0.7) +
  geom_line(aes(m, Mobile_mean), size = 2, alpha = 0.6, color = "red") +
  theme_minimal() +
  xlab("") +
  ylab("Montly applied doses against CSF")

# Saving the file
# write.csv(coverage, file = "coverage.csv")

# Fig to python ----
cov <- coverage %>%
  group_by(m = floor_date(ymd(mes), unit = "month")) %>%
  filter(m > "2017-01-01") %>%
  dplyr::summarize(
    Applied_doses = vacunados,
    x_coverage = mean(vac_cov),
    Mobile_mean = mean(mov_mean)
  )

reshape2::melt(cov, id = "m") %>%
  filter(variable != "x_coverage") %>%
  ggplot() +
  geom_line(aes(m, value, color = variable), size = 1, alpha = 0.7) +
  scale_color_manual(values = c("grey60", "#0072B166")) +
  theme_minimal() +
  xlab("") +
  ylab("Montly applied doses against CSF") +
  theme(legend.title = element_blank())

# Yearly means
coverage %>%
  group_by(year(ymd(mes))) %>%
  summarize(
    x_vaccinated_pigs = sum(vacunados),
    x_coverage = mean(vac_cov),
    x_mobile_mean = sum(mov_mean)
  )

# Looking for the infected over time ----
Fig25 <- read.csv2(file = "Beta_2014-2017_diario_Adjustment.csv", header = FALSE, sep = ";")
Fig25$month <- seq(as.Date("2014-01-01"), as.Date("2017-10-26"), by = "days")

afected2014 <- Fig25 %>%
  group_by(Month = floor_date(month, unit = "month")) %>%
  dplyr::summarise(afectados = (sum(V2, na.rm = TRUE)))

# Reading the affected
afec <- read_csv("afectados.2014-2023.csv")
afected <- rbind(afected2014[1:7, ], afec[, c(2, 3)])
plot(afected$Month, afected$afectados)

# Parametros - taxas  ###----
# Taxas em meses ----
beta <- 4.0706 # Calculated nls montly
gama <- 2.7941 # Calculated nls montly
mu <- (1 / (350 / 30)) # natural death rate
w <- 1 / (180 / 30) # omega loss of vaccine induced inmunity
theta <- 0.9 # percentage of animals persistent infected
tau <- gama / 30 # death rate of persistent infected
R0 <- beta / gama

1 / R0 # vaccine coverage needed

# Death rate in days of the p.i
1 / (gama / 30) * 30 # 21
11 * 30
1 / (gama) * 30 # number of days that the animals die in the infected

# About tau Persisten infected parameter
# 1/(gama/30)
# 45 semanas / 4 semanas por mes = 11.25 meses
# 1/(45/4) = 0.0888


# About the prevalence

# Coverage ----
c <- coverage$vac_cov
# cf coverage function
fp <- function(t) {
  {
    fp <- c[t]
  }
  return(fp)
}

fp(105) # testing

# Parametros ----
par.SIRVsd <- c(beta = beta, gama = gama, mu = mu, w = w, theta = theta, tau = tau)

# Define system ----
# Initial state vector SVIRPIR----
# calculating number of animals by S, I and V ----

# Initial conditions
mean(2582409, 3075957, 2890209)
2582409 # mean applied doses of the
# Mean Coverage
summary(coverage$vac_cov)

2582409 / 1.55 # Using the coverage adjustment to define the population (Acosta et al, 2022)

pop <- 2000000
infe <- round(((770 * 100) / 4) / 12, 0) # 770/12
p.i <- round(infe * 0.05, 0) # persistent infected
r <- infe
vac <- round((pop - infe - p.i - r) * 0.29, 0) # initial vaccinated pop
sus <- pop - vac - infe - p.i - r
infe + p.i + r + vac + sus - pop

# Apparente prevalence of CSF at initial conditions on the first month
12 * (infe + p.i) / pop # 0.098

#  Set intial vector
init.sd <- c(S = sus, I = infe, PI = p.i, V = vac, R = r)
(infe + p.i) / pop * 100

# Deterministic modell  ----
# Calling function

# Time interval 2014-2022  ----
Dt <- 1
tf <- 12 * 9 # Final time 2014-2022

# SIRV with population dynamics  ----
SIRVsd <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # rate of change
    dV <- (fp(t)) * mu * (S + I + PI + V + R) - w * V - mu * V
    dS <- (1 - (fp(t))) * mu * (S + I + PI + V + R) - (beta * S * I) / (S + I + PI + V) + w * V - mu * S
    dI <- (beta * S * I) / (S + I + PI + V + R) - gama * I * theta - gama * I * (1 - theta) - mu * I
    dPI <- gama * I * (1 - theta) - tau * PI - mu * PI
    dR <- gama * I * theta + tau * PI - mu * R
    # return the output of the model
    return(list(c(dS, dI, dPI, dV, dR)))
  })
}

times <- seq(1, tf, by = Dt)

# Simulation ----
modSIRVsd <- ode(y = init.sd, times = times, func = SIRVsd, parms = par.SIRVsd, method = "ode45")
modSIRVsd <- as.data.frame(modSIRVsd)
modSIRVsd$N <- (modSIRVsd$S + modSIRVsd$I + modSIRVsd$PI + modSIRVsd$V + modSIRVsd$R)
modSIRVsd$month <- coverage$mes

# Simulation Dynamics

# Add the cases to the simulation
modSIRVsd$cases <- afected$afectados[match(ymd(modSIRVsd$month), afected$Month)]

# FIG 3 ###ggplot#####
library(ggplot2)
library(lubridate)

cbPalette <- c("red", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

ggplot(modSIRVsd, aes(floor_date(ymd(month), unit = "month"))) +
  geom_line(aes(y = PI, colour = "P.Infected"), size = 1) +
  geom_line(aes(y = V, colour = "Vaccined"), size = 1) +
  geom_line(aes(y = S, colour = "susceptible"), size = 1) +
  geom_line(aes(y = I, colour = "Infected"), size = 1) +
  geom_line(aes(y = R, colour = "Removed"), size = 1) +
  geom_line(aes(y = N, colour = "Total"), size = 1) +
  ylab("Population") +
  xlab("Time (months)") +
  scale_color_manual(values = cbPalette) +
  labs(colour = "Simulation+
  dynamics") +
  theme_minimal()

# Writting the file to further use with python
# write.csv(modSIRVsd, file = "modsirv.csv")

ggplot(modSIRVsd, aes(floor_date(ymd(month), unit = "month"))) +
  geom_line(aes(y = PI, colour = "Persistent I"), size = 1.2) +
  geom_line(aes(y = V, colour = "Vaccined"), size = 1.2) +
  geom_line(aes(y = S, colour = "susceptible"), size = 1.2) +
  geom_line(aes(y = I, colour = "Infected"), size = 1.2) +
  ylim(0, 45000) +
  ylab("Population") +
  xlab("Time (months)") +
  scale_color_manual(values = cbPalette) +
  labs(colour = "Simulation
  dynamics") +
  theme_minimal()

# Grafico dos infectados sqrt scale on Y ----
# GGPLOT
ggplot(modSIRVsd) +
  geom_line(aes(ymd(month), I, colour = "Infected"), size = 1) +
  geom_line(aes(ymd(month), PI, colour = "P Infected"), size = 1) +
  geom_point(aes(ymd(month), cases, fill = "Observed"), size = 2, colour = "red", alpha = 0.4) +
  geom_point(aes(ymd(month), cases * 100 / 6, fill = "Predicted (Se)"), size = 1.5, shape = 17, colour = "#10A53DFF", alpha = 0.7) +
  ylim(0, 42000) +
  ylab("Number of cases") +
  xlab("Time (Months)") +
  labs(
    fill = "Cases",
    colour = "Simulation"
  ) +
  scale_y_continuous(trans = "sqrt") +
  theme_minimal()

# ----


geom_line(aes(ymd(month), PI, colour = "P Infected"), size = 1) +
  geom_point(aes(ymd(month), cases, colour = "observed"), size = 0.3, color = "red") +
  geom_point(aes(ymd(month), cases * 100 / 6), size = 0.05, colour = "#009E73") +
  ylim(0, 42000) +
  ylab("Number of cases") +
  xlab("Time (Months)") +
  labs(colour = "Compartment
    dynamics") +
  scale_y_continuous(trans = "sqrt") +
  theme_minimal()


ggplot(modSIRVsd) +
  geom_line(aes(ymd(month), I, colour = "Infected"), size = 1) +
  geom_line(aes(ymd(month), PI, colour = "P Infected"), size = 1) +
  geom_point(aes(ymd(month), cases), size = 0.5, color = "red") +
  geom_point(aes(ymd(month), cases * 100 / 6), size = 0.05, color = "black") +
  ylim(0, 4200) +
  ylab("Population") +
  xlab("Time (Months)") +
  labs(colour = "Compartment
    dynamics") +
  theme_minimal()


# Transfering
coverage$afectados <- afected$afectados[match(as.character(coverage$mes), as.character(afected$Month))]
coverage$afectados[is.na(coverage$afectados)] <- 0
modSIRVsd$afectados <- coverage$afectados

# Zoom infected ----
tiff(
  filename = "Fig4.1.tif", width = 420, height = 200, pointsize = 40,
  units = "mm", res = 400, compression = "lzw"
)

floor_date(ymd(month), unit = "month")

ggplot(modSIRVsd, aes(floor_date(ymd(month), unit = "month"))) +
  geom_line(aes(y = I, colour = "Infected"), size = 1.4) +
  geom_line(aes(y = PI, colour = "Persistent I"), size = 1.4) +
  geom_line(aes(y = (afectados * 95) / 6.87, colour = "Most probable"), size = 0.3) +
  geom_line(aes(y = afectados, colour = "Detected"), size = 0.2) +
  ylim(0, 10000) +
  ylab("Population") +
  xlab("Time (months)") +
  labs(colour = "Simulation
dynamics")

dev.off()


# Palette para os cores dos graficos

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbPalette <- c("#D55E00", "#009E73", "#56B4E9", "#999999")

cbPalette2 <- c("#D55E00", "#009E73", "#56B4E9", "#999999", "#CC79A7")


### Highest and lower values ----
# Infected
# Looking the higest and lowest values
max(modSIRVsd$I) # 67475.05
min(modSIRVsd$I) # 24.3

which.min(modSIRVsd$I) # 833 onde e o tempo t do valor menor
min(modSIRVsd$I) # 24.3 qual e o valor menor

which.max(modSIRVsd$I) # onde e o tempo t do valor menor
max(modSIRVsd$I) # qual e o valor menor

# Second epidemic curve
which.max(modSIRVsd$I) # onde e o tempo t do valor menor
max(modSIRVsd$I) # qual e o valor menor



# Persistent infected
which.max(modSIRVsd$PI) # onde e o tempo t do valor menor
max(modSIRVsd$PI) # qual e o valor menor

which.max(modSIRVsd$PI) # onde e o tempo t do valor menor
max(modSIRVsd$PI) # qual e o valor menor

# # Second epidemic curve
# I= 210 at 1294
# PI= 420 at 1396

summary(modSIRVsd$I)
# > summary(modSIRVsd$I)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 24.37    61.65   205.35  7648.41  4204.48 67475.05

summary(modSIRVsd$PI)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.0   397.5  2394.3 16585.7 22885.7 84589.7

##########
boxplot(modSIRVsd$I, modSIRVsd$PI)

###########################################
# Surveillance system Sensitivity analysis

# 2014#######################3
# Infected
summary(modSIRVsd$I[1:365])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 770    8553   21073   27853   47308   67475

# Persistent infected
summary(modSIRVsd$PI[1:365])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0    7140   53948   43959   74705   84590

27853 + 43959
71812 # suppose is the I+PI animals number

770 / 71812 #  1,07% System sensitivity

71812 / 1527114 # 4.7% prevalence

# 2017 #############################

365 * 3 <- 1095

summary(modSIRVsd$I[1095:1395])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 108.4   174.4   201.8   187.5   207.9   210.1

summary(modSIRVsd$PI[1095:1395])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 254.4   286.0   350.9   343.3   400.1   420.8

sum(Fig25$Casos[1095:1395], na.rm = TRUE)

187.5 + 343.3
530.8 # Total number of cases

102 / 530.8 # 19.2% system sentitivity

(530.8 / 2044000) * 100 # 0.025% prevalence (round)







############################################# 3

##################### GREY ####################

tiff(
  filename = "0_SIRV_SIRVip_2014-2017_OIE_B&W.tif", width = 200, height = 180, pointsize = 12,
  units = "mm", res = 300, compression = "lzw"
)

plot(modSIRVsd$time, modSIRVsd$S,
  col = "grey", lty = 1, lwd = 2.5, type = "l",
  ylim = c(0, 560000), xlab = "Time (days)", ylab = "Population"
)
lines(modSIRVsd$time, modSIRVsd$I, col = "black", lty = 1, lwd = 3.5)
lines(modSIRVsd$time, modSIRVsd$PI, col = "grey", lty = 1, lwd = 2.5)
lines(modSIRVsd$time, modSIRVsd$R, col = "black", lty = 3, lwd = 3.5)
lines(modSIRVsd$time, modSIRVsd$V, col = "grey", lty = 8, lwd = 3.5)
# lines(modSIRVsd$time, modSIRVsd$N, col="black")
title(main = "SIRV CSFV 2014 - 2017 in Ecuador", sub = "Vaccine coverage simulation")
grid(NULL, NULL, lty = 3, col = "cornsilk2")
# legend("topright", legend=c("S"," I","pI","R","V"), lty=1, cex=0.95, col=c("blue","red", "pink", "green","grey"))

dev.off()







tiff(
  filename = "1_Infectados_&_pI_2014-2017_2_OIE.tif", width = 200, height = 180, pointsize = 12,
  units = "mm", res = 300, compression = "lzw"
)


plot(afected$afectados,
  ylim = c(0, 85000), pch = 1, cex = 0.7,
  xlab = "Time (days)",
  ylab = "Number of Cases"
)
lines(modSIRVsd$time, modSIRVsd$I, lty = 1)
lines(modSIRVsd$time, modSIRVsd$PI, lty = 2)

title(main = "SIRV CSFV 2014 - 2017 in Ecuador", sub = ("Infected and persistent infection interraction"))
grid(NULL, NULL, lty = 3, col = "cornsilk2")
legend("topright", legend = c("I", "pI"), lty = c(1, 2))

dev.off()


# Zoom infected ----
plot(Fig25,
  ylim = c(0, 1000), pch = 1, cex = 0.7,
  xlab = "Time (days)",
  ylab = "Number of Cases"
)
lines(modSIRVsd$time, modSIRVsd$I, col = "red")
lines(modSIRVsd$time, modSIRVsd$PI, col = "pink")

Fig25 <- read.csv2(file = "Infected.csv", header = FALSE, sep = ";")
names(Fig25) <- c("Tempo", "Casos")



# Grafico dos infectados zoom GREY ----
tiff(
  filename = "2_Infectados_&_pI_2014-2017_2_OIE.tif", width = 200, height = 180, pointsize = 12,
  units = "mm", res = 300, compression = "lzw"
)

Fig25 <- read.csv2(file = "Beta_2014-2017_diario_Adjustment.csv", header = FALSE, sep = ";")
names(Fig25) <- c("Tempo", "Casos")

plot(Fig25,
  ylim = c(0, 800), pch = 1,
  xlab = "Time (days)",
  ylab = "Number of Cases"
)
lines(modSIRVsd$time, modSIRVsd$I, col = "black", lty = 1, lwd = 3.5)
lines(modSIRVsd$time, modSIRVsd$PI, col = "grey", lty = 2, lwd = 3.5)

# title(main = "SIRV CSFV 2014 - 2017 in Ecuador", sub = ("Infected and persistent infection interraction"))
grid(NULL, NULL, lty = 3, col = "cornsilk2")
# legend("topright", legend=c("I","pI"), col=c("red", "pink"), cex=1.2, lty=c(1, 2))

dev.off()

##############
# Recovered compartment----
tiff(
  filename = "3_Recovered_2014-2017.tif", width = 200, height = 180, pointsize = 12,
  units = "mm", res = 300, compression = "lzw"
)

Fig25 <- read.csv2(file = "Beta_2014-2017_diario_Adjustment.csv", header = FALSE, sep = ";")
names(Fig25) <- c("Tempo", "Casos")

plot(Fig25,
  ylim = c(0, 550000), pch = 1,
  xlab = "Time (days)",
  ylab = "Number of Cases"
)
lines(modSIRVsd$time, modSIRVsd$I, col = "red")
lines(modSIRVsd$time, modSIRVsd$PI, col = "pink")
lines(modSIRVsd$time, modSIRVsd$R, col = "green")

title(main = "Infected, persistent I & Recovered CSFV 2014 - 2017 in Ecuador", sub = ("Infected and persistent infection interraction"))
grid(NULL, NULL, lty = 3, col = "cornsilk2")
legend("topright", legend = c("R", "I", "pI"), lty = 1, cex = 0.95, col = c("green", "red", "pink"))


dev.off()
