###### 2014-2022 S-I-PI-V Simulation ######----
###### 2014-20 S-I-P-V Simulation ######----

library(zoo); library(deSolve); library(ggplot2); library(lubridate); library(dplyr)
# Diretorio de trabalho
setwd("~/Dropbox/0.USP/7.Publicações/Mathematical modelling CSFV/Analise/")

# Parametros - taxas  ###----
#taxas em dias ----
beta <- 4.0706 #Calculada nls
gama <- 2.7941 #Calculada nls
tau <- gama/30 #
mu <- (1/(350/30))
R0 <- beta/gama
w <- 1/(180/30) #omega
theta <- 0.1

# About tau Persisten infected parameter 
# 1/(gama/30)
# 45 semanas / 4 semanas por mes = 11.25 meses 
# 1/(45/4) = 0.0888

# Coverage ----

# v <- .9
# vc <- c(rep(c(v*.95, v*1.05),60))
# runif(120, 0.9025, 7)

rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

#95%
vc <- rnorm2(120,.95,.00625)
vc[vc > .95,] <- .95

#90%
vc <- rnorm2(120,.9,.00625)
vc[vc > .9,] <- .9

#85%
vc <- rnorm2(120,.85,.00625)
vc[vc > .85,] <- .85

#80%
vc <- rnorm2(120,.80,.00625)
vc[vc > .8,] <- .8

#75%
vc <- rnorm2(120,.75,.00625)
vc[vc > .75,] <- .75

#70%
vc <- rnorm2(120,.70,.00625)
vc[vc > .7,] <- .7

mean(vc); summary(vc)

# vaccine coverage function by years  ----
c <- vc
# cf coverage function
fc <- function(t){
  c[t]
}

fc(120) #testing

# Parametros ----
par.SIRVsd <- c(beta = beta, gama = gama, mu = mu, w = w, theta = theta, tau = tau)

# Define system ----
# Initial state vector ----
# calculating number of animals by S, I and V ----
pop <- 2000000
infe <-208 #136 #final conditions  Dec-2021
p_infe <- 520 #413 #persistent infected
r <- infe+p_infe
sus <- pop-infe-p_infe

# 80%
vac <-(sus)*(0.8) #initial vaccinated pop as jan 2022

sus <- pop-vac-r
vac+sus+r-pop

#  Set intial vector
init.sd <- c(S=sus, I=infe, P=p_infe, V=vac, R=r)                      

# Deterministic modell  ----
# Calling function
# Time interval 2014-2022  ----
Dt <- 1
tf <- 12*10   # Final time 2022-2033

# SIRV with population dynamics  ----
SIRVsd <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    
    # rate of change
    dV = (fc(t))*mu*(S+I+P+V+R) - w*V - mu*V
    dS = (1-(fc(t)))*mu*(S+I+P+V+R) - (beta*S*I)/(S+I+P+V) + w*V - mu*S
    dI = (beta*S*I)/(S+I+P+V+R) - gama*I*theta - gama*I*(1-theta) - mu*I
    dP= gama*I*(1-theta) - tau*P - mu*P
    dR= gama*I*theta + tau*P - mu*R
    # return the output of the model
    return(list(c(dS, dI, dP, dV, dR)))
  })
}

times <- seq(1, tf, by=Dt)

# Simulation ----
modSIRVsd <- ode(y = init.sd, times = times, func = SIRVsd, parms = par.SIRVsd, method = "ode45")
modSIRVsd <- as.data.frame(modSIRVsd)
modSIRVsd$N <- (modSIRVsd$S + modSIRVsd$I + modSIRVsd$P + modSIRVsd$V + modSIRVsd$R)
modSIRVsd$month <- ymd("2022-01-01")+ months(1:120)
modSIRVsd$I[modSIRVsd$I <=1 ] <- 0
modSIRVsd$P[modSIRVsd$P <=1 ] <- 0
modSIRVsd$infected <- (modSIRVsd$I + modSIRVsd$P)
modSIRVsd$I[modSIRVsd$I <=1 ] <- NA
modSIRVsd$P[modSIRVsd$P <=1 ] <- NA


# Calculating the prevalence 
# Prevalence
modSIRVsd %>% 
  filter(!is.na(infected)) %>% 
  summarise(prev.infected=round(mean(infected/N)*100,3))

modSIRVsd %>% 
  filter(!is.na(infected)) %>% 
  summarise(total.infected=mean(infected))

table(!is.na(modSIRVsd$I))
table(!is.na(modSIRVsd$P))

modSIRVsd$month[which.min(modSIRVsd$I)]
modSIRVsd$month[which.min(modSIRVsd$P)]

# Creating a df with all simulation results
# modresults <- data.frame(modSIRVsd$month)
# modresults$I95 <- modSIRVsd$I
# modresults$P95 <- modSIRVsd$P
 
# modresults$I90 <- modSIRVsd$I
# modresults$P90 <- modSIRVsd$P

modresults$I85 <- modSIRVsd$I
modresults$P85 <- modSIRVsd$P

modresults$I80 <- modSIRVsd$I
modresults$P80 <- modSIRVsd$P

modresults$I75 <- modSIRVsd$I
modresults$P75 <- modSIRVsd$P

modresults$I70 <- modSIRVsd$I
modresults$P70 <- modSIRVsd$P


modresults$x <- modSIRVsd$I
modresults$Px <- modSIRVsd$P

# Montly mean of animals by status I or IP
mean(modresults$I95[1:4])
mean(modresults$P95[1:35])

mean(modresults$I85[1:5])
mean(modresults$P85[1:35])

mean(modresults$I80)
mean(modresults$P80)

mean(modresults$I75)
mean(modresults$P75)


# Simulation Dynamics
# FIG 3 ###ggplot#####
library(ggplot2); library(lubridate)

library(reshape2)
df <- melt(modresults, id="modSIRVsd.month")

# General simulations
ggplot(df, aes(floor_date(ymd(modSIRVsd.month), unit = "month")))+ 
  geom_line(aes(y=value, colour=variable), size=2)+
  scale_colour_viridis_d()+
  ylab('Population')+
  xlab('Time (months)')+
  ylim(0, 14000)+
  labs(colour = "Infected
Simulation 
dynamics")+
  theme_linedraw()


# Mutating names
df <- df %>% 
  mutate(simulation=ifelse(variable == "I95" | variable == "P95","90-95%",NA)) %>% 
  mutate(simulation=ifelse(variable == "I90" | variable == "P90","90-95%",simulation)) %>% 
  mutate(simulation=ifelse(variable == "I85" | variable == "P85","85%",simulation)) %>% 
  mutate(simulation=ifelse(variable == "I80" | variable == "P80","80%",simulation)) %>% 
  mutate(simulation=ifelse(variable == "I75" | variable == "P75","75%",simulation)) %>% 
  mutate(simulation=ifelse(variable == "I70" | variable == "P70","70%",simulation))

table(df$simulation)

# Grouping
df$group <- rep(c(rep("Infected",120), rep("Persistent I",120)), 6)

# last plot  ----
df %>% 
  ggplot(aes(floor_date(ymd(modSIRVsd.month), unit = "month")))+ 
  geom_line(aes(y=value, colour=variable), size=1, alpha=0.2)+
  geom_point(aes(y=value, colour=variable, shape=factor(group)), size=1.5, alpha=0.7)+
  scale_colour_viridis_d(direction=-1)+
  scale_shape_manual(values = c(46, 20)) +
  facet_grid(rows = "simulation", scales="free")+
  ylab('Simulated infected population')+
  xlab(" ") +
  labs(shape = "Infected
Simulation 
dynamics")+
  theme_linedraw()+
  guides(colour="none")

write.csv(df, file = "df.csv")

df$per <- c(rep(95,240), rep(90,240), rep(85,240), rep(80,240), rep(75,240), rep(70,240))

# Montly prevalence (%)
df %>% 
  group_by(per) %>% 
  filter(!is.na(value)) %>% 
  summarise(prev.infected=round(mean(value/2000549)*100,4))

df %>% 
  group_by(year(df$modSIRVsd.month), per) %>% 
  filter(!is.na(value)) %>% 
  summarise(prev.infected=round(mean(value/2000549)*100,4))


#montly mean number of infected animals
df %>% 
  group_by(per) %>% 
  filter(!is.na(value)) %>% 
  summarise(infected=mean(value))

df %>% 
  group_by(variable) %>% 
  filter(!is.na(value)) %>% 
  summarise(infected=mean(value))



df %>% 
  group_by(variable) %>% 
  filter(!is.na(value)) %>% 
  summarise(prev.infected=round(mean(value/2000549)*100,3))


  filter(!is.na(infected)) %>% 
  summarise(prev.infected=round(mean(infected/N)*100,3))

modSIRVsd %>% 
  filter(!is.na(infected)) %>% 
  summarise(total.infected=mean(infected))




 # Grafico dos infectados ----
# GGPLOT

tiff(filename="Fig_infectados.1.tif", width=260, height=180, pointsize=12,
     units="mm", res=400, compression="lzw")

ggplot(modSIRVsd)+ 
  geom_line(aes(ymd(month), I, colour="Infected"), size=1)+ 
  geom_line(aes(ymd(month), P, colour="P Infected"), size=1)+
  geom_line(aes(ymd(month), afectados))+
  ylim(0, 10000)+
  ylab('Population')+
  xlab('Time (Months)')+
  labs(colour = "Compartment 
    dynamics")

# Transfering 
coverage$afectados <- afected$afectados[match(as.character(coverage$mes), as.character(afected$Month))]
coverage$afectados[is.na(coverage$afectados)] <- 0
modSIRVsd$afectados <- coverage$afectados

# Zoom infected ----
tiff(filename="Fig4.1.tif", width=420, height=200, pointsize=40,
     units="mm", res=400, compression="lzw")

floor_date(ymd(month), unit = "month")

ggplot(modSIRVsd, aes(floor_date(ymd(month), unit = "month")))+ 
  geom_line(aes(y=I, colour="Infected"), size=1.4)+ 
  geom_line(aes(y=PI, colour="Persistent I"), size=1.4)+
  geom_line(aes(y=(afectados*95)/6.87, colour="Most probable"), size=0.3)+
  geom_line(aes(y=afectados, colour="Detected"), size=0.2)+
  ylim(0, 10000)+
  ylab('Population')+
  xlab('Time (months)')+
    labs(colour = "Simulation 
dynamics")

dev.off()


# Palette para os cores dos graficos 

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cbPalette <- c("#D55E00", "#009E73","#56B4E9", "#999999")

cbPalette2 <- c("#D55E00", "#009E73","#56B4E9", "#999999", "#CC79A7")


###Highest and lower values ----
# Infected
#Looking the higest and lowest values 
max(modSIRVsd$I) #67475.05
min(modSIRVsd$I) #24.3

which.min(modSIRVsd$I) # 833 onde e o tempo t do valor menor
min(modSIRVsd$I) # 24.3 qual e o valor menor

which.max(modSIRVsd$I) #onde e o tempo t do valor menor
max(modSIRVsd$I) #qual e o valor menor

#Second epidemic curve
which.max(modSIRVsd$I) #onde e o tempo t do valor menor
max(modSIRVsd$I) #qual e o valor menor



# Persistent infected
which.max(modSIRVsd$P) #onde e o tempo t do valor menor
max(modSIRVsd$P) #qual e o valor menor

which.max(modSIRVsd$P) #onde e o tempo t do valor menor
max(modSIRVsd$P) #qual e o valor menor

# # Second epidemic curve
# I= 210 at 1294
# PI= 420 at 1396

summary(modSIRVsd$I)
# > summary(modSIRVsd$I)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 24.37    61.65   205.35  7648.41  4204.48 67475.05 

summary(modSIRVsd$P)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.0   397.5  2394.3 16585.7 22885.7 84589.7

##########
boxplot(modSIRVsd$I, modSIRVsd$P)

###########################################
#Surveillance system Sensitivity analysis 

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
71812 # suppose is the I+P animals number

770/71812 #  1,07% System sensitivity

71812/1527114  # 4.7% prevalence

# 2017 #############################

365*3=1095

summary(modSIRVsd$I[1095:1395])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 108.4   174.4   201.8   187.5   207.9   210.1 

summary(modSIRVsd$P[1095:1395])
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 254.4   286.0   350.9   343.3   400.1   420.8 

sum(Fig25$Casos[1095:1395], na.rm = TRUE)

187.5+343.3
530.8 #Total number of cases

102/530.8 # 19.2% system sentitivity

(530.8/2044000)*100 # 0.025% prevalence (round)







#############################################3

##################### GREY ####################

tiff(filename="0_SIRV_SIRVip_2014-2017_OIE_B&W.tif", width=200, height=180, pointsize=12,
     units="mm", res=300, compression="lzw")

plot(modSIRVsd$time, modSIRVsd$S, col="grey", lty= 1, lwd=2.5, type="l", 
     ylim=c(0,560000), xlab="Time (days)", ylab = "Population")
lines(modSIRVsd$time, modSIRVsd$I, col="black", lty= 1, lwd=3.5)
lines(modSIRVsd$time, modSIRVsd$P, col="grey", lty= 1, lwd=2.5)
lines(modSIRVsd$time, modSIRVsd$R, col="black", lty= 3, lwd=3.5)
lines(modSIRVsd$time, modSIRVsd$V, col="grey", lty= 8, lwd=3.5)
# lines(modSIRVsd$time, modSIRVsd$N, col="black")
title(main = "SIRV CSFV 2014 - 2017 in Ecuador", sub = "Vaccine coverage simulation")
grid (NULL,NULL, lty = 3, col = "cornsilk2")
# legend("topright", legend=c("S"," I","P","R","V"), lty=1, cex=0.95, col=c("blue","red", "pink", "green","grey"))

dev.off()

# Looking for the infected over time ----
Fig25 <- read.csv2(file="Beta_2014-2017_diario_Adjustment.csv", header=FALSE, sep=";")
View(Fig25)

Fig25$month <- seq(as.Date("2014-01-01"), as.Date("2017-10-26"), by = "days")

afected2014 <- Fig25 %>% 
  group_by(Month=floor_date(month, unit="month" )) %>% 
  summarise(afectados=(sum(V2, na.rm=TRUE)))

# Reading the affected
afec <- read_csv("afectados.2014-2021.csv")

afected <- rbind(afected2014[1:7,], afec[,c(2,3)])

plot(afectados$Month, afectados$afectados)

tiff(filename="1_Infectados_&_pI_2014-2017_2_OIE.tif", width=200, height=180, pointsize=12,
     units="mm", res=300, compression="lzw")


plot(afected$afectados, ylim=c(0,85000), pch=1, cex=0.7,
     xlab = "Time (days)",
     ylab = "Number of Cases")
lines(modSIRVsd$time, modSIRVsd$I, lty=1)
lines(modSIRVsd$time, modSIRVsd$PI, lty=2)

title(main = "SIRV CSFV 2014 - 2017 in Ecuador", sub = ("Infected and persistent infection interraction"))
grid (NULL,NULL, lty = 3, col = "cornsilk2")
legend("topright", legend=c("I","pI"), lty=c(1, 2))

dev.off()


# Zoom infected ----
plot(Fig25, ylim=c(0,1000), pch=1, cex=0.7,
     xlab = "Time (days)",
     ylab = "Number of Cases")
lines(modSIRVsd$time, modSIRVsd$I, col="red")
lines(modSIRVsd$time, modSIRVsd$PI, col="pink")

Fig25 <- read.csv2(file="Infected.csv", header=FALSE, sep=";")
names(Fig25) <- c("Tempo","Casos")



# Grafico dos infectados zoom GREY ----
tiff(filename="2_Infectados_&_pI_2014-2017_2_OIE.tif", width=200, height=180, pointsize=12,
     units="mm", res=300, compression="lzw")

Fig25 <- read.csv2(file="Beta_2014-2017_diario_Adjustment.csv", header=FALSE, sep=";")
names(Fig25) <- c("Tempo","Casos")

plot(Fig25, ylim=c(0,800), pch=1,
     xlab = "Time (days)",
     ylab = "Number of Cases")
lines(modSIRVsd$time, modSIRVsd$I, col="black", lty=1, lwd=3.5)
lines(modSIRVsd$time, modSIRVsd$PI, col="grey", lty=2, lwd=3.5)

# title(main = "SIRV CSFV 2014 - 2017 in Ecuador", sub = ("Infected and persistent infection interraction"))
grid (NULL,NULL, lty = 3, col = "cornsilk2")
# legend("topright", legend=c("I","pI"), col=c("red", "pink"), cex=1.2, lty=c(1, 2))

dev.off()

##############
#Recovered compartment----
tiff(filename="3_Recovered_2014-2017.tif", width=200, height=180, pointsize=12,
     units="mm", res=300, compression="lzw")

Fig25 <- read.csv2(file="Beta_2014-2017_diario_Adjustment.csv", header=FALSE, sep=";")
names(Fig25) <- c("Tempo","Casos")

plot(Fig25, ylim=c(0,550000), pch=1,
     xlab = "Time (days)",
     ylab = "Number of Cases")
lines(modSIRVsd$time, modSIRVsd$I, col="red")
lines(modSIRVsd$time, modSIRVsd$PI, col="pink")
lines(modSIRVsd$time, modSIRVsd$R, col="green")

title(main = "Infected, persistent I & Recovered CSFV 2014 - 2017 in Ecuador", sub = ("Infected and persistent infection interraction"))
grid (NULL,NULL, lty = 3, col = "cornsilk2")
legend("topright", legend=c("R","I","pI"), lty=1, cex=0.95, col=c("green","red", "pink"))


dev.off()
