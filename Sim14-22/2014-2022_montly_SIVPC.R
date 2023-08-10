###### 2014-2022 S-I-V-P-C Simulation ######----

#Working directory
setwd("~/papers/modelling_csf/Modelling-imnunization-csf-ecuador/Sim14-22")

# Dependencies
library(deSolve)
library(ggplot2)
library(lubridate)
library(dplyr)
library(zoo)
library(readr)
library(ggpubr)

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
afec <- read.csv(file="afectados.2014-2023.csv") %>% select(-1)
afected <- rbind(afected2014[1:7, ], afec)
plot(afected$Month, afected$afectados)

# Simulation
# Parameters ----

# Dependencies
library(ggplot2)
library(lubridate)

coverage2 <- coverage[,-3]
colnames(coverage2) <- c("date", "number_of_vaccinated", "rollmean", "coverage")
write.csv(coverage2 ,file = "coverage.csv")

# Palette 
# cbPalette <- c("red", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
# 
# cbPalette <- c("#0072B2", "red", "#E69F00", "#009E73", "#56B4E9", "#F0E442")

cbPalette <- c("#56B4E9","#440154FF", "#22A884FF","#414487FF", "#FDE725FF","#7AD151FF")

# Monthly parameters ----
beta <- 1.1664 #best
gama <- 0.6939 #best
betab <- beta * 0.6 # Assumption
mu <- (1 / (350 / 30)) # natural death rate
w <- 1 / (180 / 30) # omega loss of vaccine induced immunity
theta <- 0.15 # percentage of animals persistent infected
tau <- gama / 30 # death rate of persistent infected
d <- (1/(10/30))  # Monthly number of days to natural decay of carcass to be used as swill
r <- (1/(15/30))  # Monthly number of days to remove the carcass

R0 <- beta / gama

# Checking days
1/(beta / 30)
1/(betab / 30)
1/(gama / 30)
1/(mu / 30)
1/(w / 30)
tau <- gama / 34 # death rate of persistent infected
(1/(tau / 30))/30
1/(r / 30)  # 4
1/(d / 30)  # 4

# About the prevalence
# Coverage ----
c <- coverage$vac_cov * 0.9

# fc coverage function
fc <- function(t) {
      c[t]
}

fc(108) # testing

# Parametros ----
par.SIRVsd <- c(beta = beta,
                betab = betab,
                gama = gama,
                mu = mu,
                w = w,
                theta = theta,
                tau = tau,
                d = d,
                r = r)

# Define system ----
# Initial state vector SVIRPIR----

# Initial conditions
mean(2582409, 3075957, 2890209)
2582409 # mean applied doses of the
# Mean Coverage
summary(coverage$vac_cov)

2582409 / 1.55 # Using the coverage adjustment to define the population (Acosta et al, 2022)

pop <- 2000000
infe <- round(((770/12 * 100) / 4), 0) # 770/12
p.i <- round(infe * 0.2, 0) # persistent infected
vac <- round((pop - infe - p.i) * 0.27, 0) # initial vaccinated pop
car <- (infe + p.i) * 0.2
sus <- pop - vac - infe - p.i
infe + p.i + vac + sus - pop

# Apparent prevalence of CSF at initial conditions on the first month
12 * (infe + p.i) / pop # 0.098

#  Set initial vector
init.sd <- c(S = sus,
             I = infe,
             P = p.i,
             V = vac,
             C = car)
init.sd
# Deterministic model  ----
# Calling function

# Time interval 2014-2022  ----
Dt <- 1
tf <- 12 * 8 # Final time 2014-2021
times <- seq(1, tf, by = Dt)

# SIRV with population dynamics  ----
SIRVsd <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {

    dV <- fc(t)*mu*(V+S+I+P+C) - w*mu*V - mu*V -(betab*V*C)/(V+S+I+P+C)
    dS <- (1 - fc(t))*mu*(V+S+I+P+C) + w*mu*V - (beta*S*I)/(V+S+I+P+C) - mu*S - (betab*C*S)/(V+S+I+P+C)
    dI <- (beta*S*I)/(V+S+I+P+C) - (1-theta)*gama*I - theta*gama*I + (betab*C*S)/(V+S+I+P+C) - mu * I
    dP <- theta*gama*I - tau*P - mu*P + (betab*V*C)/(V+S+I+P+C)
    dC <- tau*P + (1-theta)*gama*I - d*C - r*C
    
    # return the output of the model
    return(list(c(dS, dI, dP, dV, dC)))
  })
}

# Simulation ----
modSIRVsd <- ode(y = init.sd, times = times, func = SIRVsd, parms = par.SIRVsd, method = "ode45")
# plot(modSIRVsd)
# dev.off()
# matplot.0D(modSIRVsd, lwd = 2)

modSIRVsd <- as.data.frame(modSIRVsd)
modSIRVsd$N <- (modSIRVsd$S + modSIRVsd$I + modSIRVsd$P + modSIRVsd$V)
modSIRVsd$In <- (modSIRVsd$I + modSIRVsd$P)
modSIRVsd$month <- coverage$mes[1:96]

# Coverage simulated
modSIRVsd$simCov <- modSIRVsd$V/modSIRVsd$N

# Simulation Dynamics
# Add the cases to the simulation
modSIRVsd$cases <- afected$afectados[match(ymd(modSIRVsd$month), afected$Month)]

afec_act <- read.csv(file = "afectados2014-2023-actualizacion.csv")
modSIRVsd$cases2 <- afec_act$afectados[match(modSIRVsd$month, afec_act$Month)]
afected2014 <- afected2014[1:9,]
modSIRVsd$cases2[1:9] <- afected2014$afectados

modSIRVsd$I <- round(modSIRVsd$I,0)
modSIRVsd$P <- round(modSIRVsd$P,0)
modSIRVsd$C <- round(modSIRVsd$C,0)
modSIRVsd$S <- round(modSIRVsd$S,0)

# Writting the file to further use with python
# modSIRVsd$cases <- NULL
# colnames(modSIRVsd)[11] <- "cases"
# 
# modSIRVsd$I[96] <-0
# write.csv(modSIRVsd, file = "modsirv.csv")

modSIRVsd$P[modSIRVsd$P <= 0] <- NA
modSIRVsd$I[modSIRVsd$I == 0] <- NA
modSIRVsd$C[modSIRVsd$C <= 0] <- NA
modSIRVsd$S[modSIRVsd$S <= 0] <- NA

modSIRVsd$In <- round(modSIRVsd$In,0)
modSIRVsd$In[modSIRVsd$In <= 0] <- NA

modSIRVsd$V <- round(modSIRVsd$V,0)
modSIRVsd$V[modSIRVsd$V <= 0] <- NA

# FIG 3 ###ggplot#####
ggplot(modSIRVsd, aes(floor_date(ymd(month), unit = "month"))) +
  geom_line(aes(y = P, colour = "P.Infected"), size = 2) +
  geom_line(aes(y = V, colour = "Vaccined"), size = 2) +
  geom_line(aes(y = S, colour = "susceptible"), size = 2) +
  geom_line(aes(y = I, colour = "Infected"), size = 2) +
  geom_line(aes(y = C, colour = "Carcases"), size = 0.4, color="black") +
  geom_line(aes(y = N, colour = "Total"), size = 2) +
  ylab("Population") +
  xlab("Time (months)") +
  scale_color_manual(values = cbPalette) +
  labs(colour = "Simulation 
  dynamics") +
  theme_minimal()


# Grafico dos infectados sqrt scale on Y ----
# GGPLOT
ggplot(modSIRVsd) +
  geom_line(aes(ymd(month), I, colour = "Infected"), size = 1.5) +
  geom_line(aes(ymd(month), P, colour = "Persistent"), size = 1.5) +
  geom_line(aes(ymd(month), C, colour = "Carcass"), size = 0.5, colour="black", lty=3) +
  geom_point(aes(ymd(month), cases, fill = "Observed"), size = 2, colour = "red", alpha = 0.4) +
  geom_point(aes(ymd(month), cases * 100 / 4, fill = "Predicted (Se)"), size = 1.5, shape = 17, colour = "#10A53DFF", alpha = 0.7) +
  ylab("Number of cases") +
  xlab("Time (Months)") +
  labs(
    fill = "Cases",
    colour = "Simulation"
  ) +
  scale_y_continuous(trans = "sqrt") +
  theme_minimal()

# ----
# Table simulation
summary(modSIRVsd$S)
summary(modSIRVsd$I)
summary(modSIRVsd$P)
summary(modSIRVsd$V)
summary(modSIRVsd$simCov)
summary(modSIRVsd$C)
summary(modSIRVsd$N)

summary(modSIRVsd$I[84:96])
summary(modSIRVsd$P[84:96])

mean(modSIRVsd$In[84:96])/mean(modSIRVsd$N[84:96]) 
mean(modSIRVsd$In[1:12])/mean(modSIRVsd$N[1:12])



### Highest and lower values ----
# Infected
# Looking the higest and lowest values
max(modSIRVsd$I, na.rm=TRUE) # 67475.05
min(modSIRVsd$I) # 24.3

which.min(modSIRVsd$I) # 833 onde e o tempo t do valor menor
min(modSIRVsd$I) # 24.3 qual e o valor menor

which.max(modSIRVsd$I) # onde e o tempo t do valor menor
max(modSIRVsd$I) # qual e o valor menor

# Second epidemic curve
which.max(modSIRVsd$I) # onde e o tempo t do valor menor
max(modSIRVsd$I) # qual e o valor menor



# Persistent infected
which.max(modSIRVsd$P) # onde e o tempo t do valor menor
max(modSIRVsd$P) # qual e o valor menor

which.max(modSIRVsd$P) # onde e o tempo t do valor menor
max(modSIRVsd$P) # qual e o valor menor

# # Second epidemic curve
# I= 210 at 1294
# PI= 420 at 1396

summary(modSIRVsd$I)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.00     1.75   145.00  6721.80  3317.00 77413.00

summary(modSIRVsd$P)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 2.0    110.5   3869.0  43244.4  45920.2 274257.0 

summary(modSIRVsd$C)
# 0.0     0.0    24.0   740.9   424.0  7468.0 

##########
boxplot(modSIRVsd$I, modSIRVsd$P)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0     1.0    31.0   829.3   537.2  7824.0 

###########################################
# Surveillance system Sensitivity analysis

# 2014#######################3
# Infected
summary(modSIRVsd$I[1:365])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 770    8553   21073   27853   47308   67475

# Persistent infected
summary(modSIRVsd$P[1:365])
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
  xlab = "Time (days)", ylab = "Population"
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
  ylim = c(0, 95000), pch = 1, cex = 0.7,
  xlab = "Time (days)",
  ylab = "Number of Cases"
)
lines(modSIRVsd$time, modSIRVsd$I, lty = 1)
lines(modSIRVsd$time, modSIRVsd$P, lty = 2)

title(main = "SIRV CSFV 2014 - 2017 in Ecuador", sub = ("Infected and persistent infection interraction"))
grid(NULL, NULL, lty = 3, col = "cornsilk2")
legend("topright", legend = c("I", "pI"), lty = c(1, 2))

dev.off()

