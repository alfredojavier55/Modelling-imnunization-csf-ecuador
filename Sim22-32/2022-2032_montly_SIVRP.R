###### 2014-2022 S-I-PI-V Simulation ######----

library(zoo);
library(deSolve);
library(ggplot2);
library(lubridate);
library(dplyr)
setwd("~/papers/modelling_csf/Modelling-imnunization-csf-ecuador/Sim22-32/")

# Parameters ----
beta <- 1.1664 # Calculated nls monthly
betab <- 4.0706 * 0.5 # Assumption
gama <- 0.6939 # Calculated nls monthly                
mu <- (1 / (350 / 30)) # natural death rate
w <- 1 / (180 / 30) # omega loss of vaccine induced immunity
theta <- 0.35 # percentage of animals persistent infected
tau <- gama / 30 # death rate of persistent infected
# d <- (-30*log(0.001))/7  # removal of carcass 7 days
# r <- (-30*log(0.001))/1  # removal of carcass 7 days

d <- (1/(5/30))  # Monthly number of days to remove the carcass
r <- (1/(1/30))  # Monthly number of days to natural decay of carcass to be used as swill

R0 <- beta/gama

# Parameters ----
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
# Initial state vector ----
# calculating number of animals by S, I and V ----
pop <- 2000000

# Actualizar con la informacion de casos de PPC ----
infe <- round(((551/12 * 100) / 6), 0) # 770/12
p.i <- round(infe * 0.5, 0) # persistent infected
car <- (infe + p.i) * 0.5
sus <- pop-infe-p.i-car

# 80%
mean(coverage$vac_cov[64:70]) # Coverage of the last semester of 2021
vac <-(sus)*(0.85)*0.95 #initial vaccinated pop as jan 2022
sus <- pop-vac-infe-p.i-car
infe + p.i + vac + sus + car - pop
vac/(pop)

#  Set initial vector
init.sd <- c(S=sus, 
             I=infe, 
             P=p.i, 
             V=vac,
             C=car)                      
init.sd

# Deterministic model  ----
# Calling function
# Time interval 2014-2022  ----
Dt <- 1
tf <- 12*10   # Final time 2022-2033
times <- seq(1, tf, by=Dt)

# SIRV with population dynamics  ----
SIRVsd <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    # rate of change

    dV <- fc(t)*mu*(V+S+I+P+C) - w*mu*V - mu*V -(betab*V*C)/(V+S+I+P+C)
    dS <- (1 - fc(t))*mu*(V+S+I+P+C) + w*mu*V - (beta*S*I)/(V+S+I+P+C) - mu*S - (betab*C*S)/(V+S+I+P+C)
    dI <- (beta*S*I)/(V+S+I+P+C) - (1-theta)*gama*I - theta*gama*I + (betab*C*S)/(V+S+I+P+C) - mu * I
    dP <- theta*gama*I - tau*P - mu*P + (betab*V*C)/(V+S+I+P+C)
    dC <- tau*P + (1-theta)*gama*I - d*C - r*C
    
        # return the output of the model
    return(list(c(dS, dI, dP, dV, dC)))
  })
}

# Coverage ----
rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

#50%
vc <- rnorm2(120,.50,.00625)
vc[vc > .5,] <- .5

#60%
vc <- rnorm2(120,.60,.00625)
vc[vc > .6,] <- .6

# # #70%
vc <- rnorm2(120,.70,.00625)
vc[vc > .7,] <- .7

# # #75%
vc <- rnorm2(120,.75,.00625)
vc[vc > .75,] <- .75

# #80%
vc <- rnorm2(120,.80,.00625)
vc[vc > .8,] <- .8

# #85%
vc <- rnorm2(120,.85,.00625)
vc[vc > .85,] <- .85

#90%
vc <- rnorm2(120,.9,.00625)
vc[vc > .9,] <- .9

# #95%
vc <- rnorm2(120,.95,.00625)
vc[vc > .95,] <- .95


mean(vc); summary(vc)
# vaccine coverage function by years  ----
c <- vc
# cf coverage function
fc <- function(t){
  c[t]
}

fc(120) #testing
# Simulation ----
modSIRVsd <- ode(y = init.sd, times = times, func = SIRVsd, parms = par.SIRVsd, method = "ode45")
modSIRVsd <- as.data.frame(modSIRVsd)
modSIRVsd$N <- (modSIRVsd$S + modSIRVsd$I + modSIRVsd$P + modSIRVsd$V)
modSIRVsd$coverage <- modSIRVsd$V/modSIRVsd$N
modSIRVsd$month <- ymd("2022-01-01")+ months(1:120)
modSIRVsd$I[modSIRVsd$I <=0.5 ] <- 0
modSIRVsd$P[modSIRVsd$P <=0.5 ] <- 0
modSIRVsd$C[modSIRVsd$C <=0.5 ] <- 0
modSIRVsd$infected <- (modSIRVsd$I + modSIRVsd$P)
modSIRVsd$I[modSIRVsd$I == 0 ] <- NA
modSIRVsd$P[modSIRVsd$P == 0 ] <- NA
modSIRVsd$C[modSIRVsd$C == 0 ] <- NA

compare <- modSIRVsd

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
modresults <- data.frame(modSIRVsd$month)
modresults$I95 <- modSIRVsd$I
modresults$P95 <- modSIRVsd$P
#
modresults$I90 <- modSIRVsd$I
modresults$P90 <- modSIRVsd$P

modresults$I85 <- modSIRVsd$I
modresults$P85 <- modSIRVsd$P

modresults$I80 <- modSIRVsd$I
modresults$P80 <- modSIRVsd$P

modresults$I75 <- modSIRVsd$I
modresults$P75 <- modSIRVsd$P

modresults$I70 <- modSIRVsd$I
modresults$P70 <- modSIRVsd$P

modresults$I60 <- modSIRVsd$I
modresults$P60 <- modSIRVsd$P

modresults$I50 <- modSIRVsd$I
modresults$P50 <- modSIRVsd$P

# Monthly mean of animals by status I or IP
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

# Mutating names
df$variable <- paste(df$variable, "%")
table(df$variable)

# General simulations
ggplot(df, aes(floor_date(ymd(modSIRVsd.month), unit = "month")))+ 
  geom_line(aes(y=value, colour=variable), alpha=0.5, size=1,
            position=position_dodge2(width=2,
                                    padding=1.5,
                                    preserve = "single"))+
  scale_colour_viridis_d()+
  scale_y_continuous(trans = "sqrt") +
  ylab('Population')+
  xlab('Time (months)')+
  # ylim(0, 100)+
  labs(colour = "Infected
Simulation 
dynamics")+
  theme_linedraw()

# Grouping
df$group <- rep(c(rep("Infected",120), rep("Persistent I",120)), 8)
df$group <- rep(c(rep("Infected",120), rep("Persistent I",120)), 5)
df$perc <- substr(df$variable,2,5)


# last plot  ----
df %>% 
  filter(variable != "I70 %") %>% 
  filter(variable != "P70 %") %>% 
  filter(variable != "I85 %") %>% 
  filter(variable != "P85 %") %>% 
  filter(variable != "I75 %") %>% 
  filter(variable != "P75 %") %>% 
  filter(variable != "I60 %") %>%
  filter(variable != "P60 %") %>%
  ggplot(aes(floor_date(ymd(modSIRVsd.month), unit = "month"))) + 
  geom_line(aes(y=value, colour=variable), size=2, alpha=1) +
  geom_point(aes(y=value, colour=variable, shape=factor(group)), size=2, alpha=1) +
  scale_colour_viridis_d(direction=-1) +
  scale_shape_manual(values = c(46, 20)) +
  facet_grid(rows = "variable", scales = "free") +
  ylab('Simulated infected population') +
  xlab(NULL)+
  labs(shape = "Infected
Simulation 
dynamics") +
  theme_linedraw() +
  guides(colour="none") +
  guides(shape="none")


write.csv(df, file = "df.csv")
class(df$modSIRVsd.month)

df <- read.csv(file="df.csv")
df$modSIRVsd.month <- date(df$modSIRVsd.month)

df %>% 
  filter(variable != "I70 %") %>% 
  filter(variable != "P70 %") %>% 
  filter(variable != "I85 %") %>% 
  filter(variable != "P85 %") %>% 
  filter(variable != "I75 %") %>% 
  filter(variable != "P75 %") %>% 
  filter(variable != "I60 %") %>% 
  filter(variable != "P60 %") %>% 
  ggplot(aes(floor_date(ymd(modSIRVsd.month), unit = "month"))) + 
  geom_line(aes(y=value, colour=variable), size=2, alpha=1) +
  scale_colour_viridis_d(direction=-1) +
  facet_grid(rows = "variable", scales = "free") +
  ylab('Simulated infected population') +
  xlab(NULL)+
  theme_linedraw() +
  guides(colour="none") +
  guides(shape="none")



# Montly prevalence (%)
df %>% 
  group_by(perc) %>% 
  filter(!is.na(value)) %>% 
  summarise(prev.infected=round(mean(value/2000000)*100,4))

df %>% 
  group_by(year(df$modSIRVsd.month),group) %>% 
  filter(!is.na(value)) %>% 
  summarise(inf=sum(value))

df %>% 
  group_by(year(df$modSIRVsd.month), variable) %>% 
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


# modresults_before <-modresults

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
