cbPalette <- c("red", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")
setwd("~/papers/modelling_csf/Modelling-imnunization-csf-ecuador/Sim14-22")

anterior <- modSIRVsd

# Monthly parameters ----
beta <- 1.1664 #best
betab <- beta * 0.6 # Assumption
gama <- 0.6939 #best
mu <- (1 / (350 / 30)) # natural death rate
w <- 1 / (180 / 30) # omega loss of vaccine induced immunity
theta <- 0.15 # percentage of animals persistent infected
tau <- gama / 30 # death rate of persistent infected
# d <- (-30*log(0.001))/10  #natural decay rate 5 days
# r <- (-30*log(0.001))/15  # removal of carcass 15 days
beta
gama
# carcass
d <- (1/(10/30))  # Monthly number of days to natural decay of carcass to be used as swill
r <- (1/(15/30))  # Monthly number of days to remove the carcass

R0 <- beta / gama

# Checking days
1/(beta / 30)
1/(betab / 30)
1/(gama / 30)
1/(mu / 30)
1/(w / 30)
1/(tau / 30)
322/30

1/(d / 30) # 5 dias original
1/(r / 30)  # 1 dia original

1 / R0 # vaccine coverage needed

# Death rate in days of the p.i
1 / (gama / 19.5) * 30 # 21
7 * 30
1 / (gama) * 30 # number of days that the animals die in the infected

# About the prevalence
# Coverage ----
c <- coverage$vac_cov * 0.90

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
                r = r
                )

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

infe <- round(((770/12 * 100) / 4), 0) # 770/12
p.i <- round(infe * 0.2, 0) # persistent infected
vac <- round((pop - infe - p.i) * 0.27, 0) # initial vaccinated pop
car <- (infe + p.i) * 0.2
# sus <- pop - vac - infe - p.i - r
sus <- pop - vac - infe - p.i
infe + p.i + vac + sus - pop
vac/pop
# Apparente prevalence of CSF at initial conditions on the first month
12 * (infe + p.i) / pop # 0.098

#  Set intial vector
# init.sd <- c(S = sus, I = infe, P = p.i, V = vac, R = r)
init.sd <- c(S = sus,
             I = infe,
             P = p.i,
             V = vac,
             C = car)
init.sd
(infe + p.i) / pop * 100

# Deterministic model  ----
# Calling function

# Time interval 2014-2022  ----
Dt <- 1
tf <- 12 * 8 # Final time 2014-2022
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
modSIRVsd <- ode(y = init.sd, 
                 times = times, 
                 func = SIRVsd, 
                 parms = par.SIRVsd, 
                 method = "ode45")
# plot(modSIRVsd)
# dev.off()
# matplot.0D(modSIRVsd, lwd = 2)

modSIRVsd <- as.data.frame(modSIRVsd)
modSIRVsd$N <- (modSIRVsd$S + modSIRVsd$I + modSIRVsd$P + modSIRVsd$V)
modSIRVsd$In <- (modSIRVsd$I + modSIRVsd$P)
modSIRVsd$month <-  coverage$mes[1:96]

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
modSIRVsd$I[modSIRVsd$I == 0] <- NA

modSIRVsd$P <- round(modSIRVsd$P,0)
modSIRVsd$P[modSIRVsd$P <= 0] <- NA

modSIRVsd$C <- round(modSIRVsd$C,0)
modSIRVsd$C[modSIRVsd$C <= 0] <- NA

modSIRVsd$S <- round(modSIRVsd$S,0)
modSIRVsd$S[modSIRVsd$S <= 0] <- NA

modSIRVsd$In <- round(modSIRVsd$In,0)
modSIRVsd$In[modSIRVsd$In <= 0] <- NA

modSIRVsd$V <- round(modSIRVsd$V,0)
modSIRVsd$V[modSIRVsd$V <= 0] <- NA

# FIG 3 ###ggplot#####
ga <- ggplot(modSIRVsd, aes(floor_date(ymd(month), unit = "month"))) +
  geom_line(aes(y = P, colour = "Persistent"), size = 1) +
  geom_line(aes(y = V, colour = "Vaccined"), size = 1) +
  geom_line(aes(y = S, colour = "Susceptible"), size = 1) +
  geom_line(aes(y = I, colour = "Infected"), size = 1) +
  geom_line(aes(y = C, colour = "Carcases"), size = 0.4) +
  geom_line(aes(y = N, colour = "Total"), size = 1) +
  ylab("Population") +
  xlab("Time (months)") +
  scale_color_manual(values = cbPalette) +
  labs(colour = "Simulation 
  dynamics") +
  theme_minimal()

# Writting the file to further use with python
# write.csv(modSIRVsd, file = "modsirv.csv")


# Grafico dos infectados sqrt scale on Y ----
# GGPLOT
gb <- ggplot(modSIRVsd) +
  geom_line(aes(ymd(month), I, colour = "Infected"), size = 1) +
  geom_line(aes(ymd(month), P, colour = "Persistent"), size = 1) +
  geom_line(aes(ymd(month), C, colour = "Infected"), size = 0.2, color="black") +
  geom_point(aes(ymd(month), cases2, fill = "Observed"), size = 2, colour = "orange", alpha = 0.4) +
  geom_point(aes(ymd(month), cases2 * 100 / 4, fill = "Predicted (Se)"), size = 1.5, shape = 17, colour = "#10A53DFF", alpha = 0.7) +
  stat_smooth(aes(ymd(month), cases2 * 100 / 4, method = "lm", col = "lm(predicted)")) +
  ylim(0, 42000) +
  ylab("Number of cases") +
  xlab("Time (Months)") +
  labs(
    fill = "Cases",
    colour = "Simulation"
  ) +
  scale_y_continuous(trans = "sqrt") +
  theme_minimal()

ggarrange(ga, gb, ncol=1, heights = c(1,3))

tail(modSIRVsd)

# gb
# library(ggpubr)
ggplotly(gb)

gc <- ggplot(modSIRVsd) +
  geom_line(aes(ymd(month), In, colour = "cases"), size = 1) +
  geom_line(aes(ymd(month), C, colour = "Infected"), size = 0.2, color="black") +
  geom_point(aes(ymd(month), cases2, fill = "Observed"), size = 2, colour = "orange", alpha = 0.4) +
  geom_point(aes(ymd(month), cases2 * 100 / 4, fill = "Predicted (Se)"), size = 1.5, shape = 17, colour = "#10A53DFF", alpha = 0.7) +
  stat_smooth(aes(ymd(month), cases2 * 100 / 4, method = "lm", col = "lm(predicted)")) +
  ylim(0, 42000) +
  ylab("Number of cases") +
  xlab("Time (Months)") +
  labs(
    fill = "Cases",
    colour = "Simulation"
  ) +
  scale_y_continuous(trans = "sqrt") +
  theme_minimal()


modSIRVsd %>% group_by(year(month)) %>% 
  summarise(cases=sum(cases2, na.rm = TRUE))

# ----

summary(modSIRVsd$I)
summary(modSIRVsd$P)
summary(modSIRVsd$C)
summary(modSIRVsd$N)

mean(modSIRVsd$In[84:96])/mean(modSIRVsd$N[84:96])
mean(modSIRVsd$In[1:12])/mean(modSIRVsd$N[1:12])

summary(modSIRVsd$I[84:96])
summary(modSIRVsd$P[84:96])

summary(anterior$N)
summary(modSIRVsd$N)

summary(anterior$In)
summary(modSIRVsd$In)
