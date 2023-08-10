# Adjustment data start on 2015-01-01 until 2021-12-01
setwd("~/papers/modelling_csf/Modelling-imnunization-csf-ecuador/Contact-death/")

# Libraries ----
library(deSolve)
source("ImodSIsd.R")
library(ggplot2)

# Data ----
afec_act <- read.csv(file = "afectados2014-2023-actualizacion.csv")

# Select 2015-2023 data ----
case_sim <- afec_act[c(7:98),]
colnames(case_sim)[3] <- "cases"
case_sim <- case_sim[case_sim$cases != 0, ]
case_sim_date <-  case_sim[,2]
case_sim <- case_sim[c(1,3)]
case_sim$time <- 1:nrow(case_sim)
case_sim <- case_sim[, c(3,2)]

# Cases ploting ----
plot(case_sim, pch=1, ylim=c(0,350),
     xlab = "time (month)",
     ylab = "number of cases")

# Initial conditions ----
N <- 1686
I <- 39
state.SIsd <- c(s=round((N-I)/N,4), i=round(I/N,4))
state.SIsd

# Simulation time ----
tsim <- nrow(case_sim)
Dt <- 1

# Parameters adjustment ----
aj.SIsd <- nls(cases ~ ImodSIsd(time, state.SIsd, N, tsim, Dt, beta, gama),
                data = case_sim, algorithm="port",
                start = list(beta=1, gama=1), # start
                # lower = list(beta=0, gama=0),
                # upper = list(beta=100, gama =100),
                trace = TRUE)

# Adjustment summary ----
summary(aj.SIsd)

# Coef adjusted ----
coef(aj.SIsd)

# Infected estimated by the model ----
Dt <- 1
tempos.aju2 <- seq(0, tsim, by=Dt)
Imodelo2 <- ImodSIsd(tempos.aju2,
                     state.SIsd, 
                     N = 1686, 
                     tsim, 
                     Dt, 
                     beta = coef(aj.SIsd)[1],
                     gama = coef(aj.SIsd)[2])
tail(Imodelo2)

# Add model results to observed plot ----
lines(tempos.aju2, Imodelo2, lty=3, col="blue")


t <- data.frame(tempo=tempos.aju2, imode=Imodelo2)

ggplot()+ 
  geom_point(data=case_sim, aes(x=time, y=cases), col= "BLACK", size= 1.5)+ 
  geom_line(data=t, aes(x=tempo, y=imode, col ="Estimated")) + 
  ylab('Number of cases 2014-2022')+
  xlab('Time (month)') +
  theme_light()

case_sim$date <- case_sim_date
case_sim$imode <-Imodelo2[c(-1)]
  
# Save data
write.csv(t, file = "t.csv")
write.csv(case_sim, file = "d.csv")

d <- read.csv(file="d.csv")
t <- read.csv(file="t.csv")

library(ggplot2)
library(lubridate)
ggplot()+ 
  geom_point(data=d, aes(x=ymd(date), y=cases), col= "BLACK", size= 1.5)+ 
  geom_line(data=d, aes(x=ymd(date), y=imode, col ="Adjusted")) + 
  labs(colour = "") +
  ylab('Number of cases 2014-2022')+
  xlab('Time (month)') +
  theme_light()
