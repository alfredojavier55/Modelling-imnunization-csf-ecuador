SIRVP_se <- function (beta, betab, gama, mu, w, theta, tau, d, r, tf)
  {
# Parameters----
#Monthly parameters ----

beta <- 1.1664 #best
# beta <- 4.0706 # Calculated nls monthly
betab <- beta * 0.4 # Assumption
gama <- 0.6937 #best
# gama <- 2.7941 # Calculated nls monthly                
mu <- (1 / (350 / 30)) # natural death rate
w <- 1 / (180 / 30) # omega loss of vaccine induced immunity
theta <- 0.15 # percentage of animals persistent infected
tau <- gama / 30 # death rate of persistent infected
# d <- (-30*log(0.001))/7  # removal of carcass 7 days
# r <- (-30*log(0.001))/1  # removal of carcass 7 days
d <- (1/(10/30))  # Monthly number of days to remove the carcass
r <- (1/(15/30))  # Monthly number of days to natural decay of carcass to be used as swill


# Coverage ----
rnorm2 <- function(n,mean,sd) { mean+sd*scale(rnorm(n)) }

# 90%
# vc <- rnorm2(120,.9,.00625)
# vc[vc > .9,] <- .9

# #75%
vc <- rnorm2(120,.75,.00625)
vc[vc > .75,] <- .75

# vaccine coverage function by years  ----
c <- vc
# cf coverage function
fc <- function(t){
  c[t]
}

fc(10)

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
infe <- round(((551/12 * 100) / 6), 0) # 770/12
p.i <- round(infe * 0.5, 0) # persistent infected
# infe <-10 #136 #final conditions  Dec-2021
# p_infe <- 20 #413 #persistent infected
sus <- pop-infe-p.i
car <- (infe + p.i) * 0.5

# 80%
vac <-(sus)*(0.75)*0.95 #initial vaccinated pop as jan 2022
sus <- pop-vac-infe-p.i

infe + p.i + vac + sus -pop

#  Set initial vector
init.sd <- c(S=sus,
             I=infe,
             P=p.i, 
             V=vac, 
             C=car)                      

# Deterministic model  ----
# Calling function
# Time interval 2014-2022  ----
Dt <- 1
tf <- 12*1   # Final time 2022-2033

# SIRV with population dynamics  ----
SIRVsd <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    
    # Rate of change
    dV <- fc(t)*mu*(V+S+I+P+C) - w*mu*V - mu*V -(betab*V*C)/(V+S+I+P+C)
    dS <- (1 - fc(t))*mu*(V+S+I+P+C) + w*mu*V - (beta*S*I)/(V+S+I+P+C) - mu*S - (betab*C*S)/(V+S+I+P+C)
    dI <- (beta*S*I)/(V+S+I+P+C) - (1-theta)*gama*I - theta*gama*I + (betab*C*S)/(V+S+I+P+C) - mu * I
    dP <- theta*gama*I - tau*P - mu*P + (betab*V*C)/(V+S+I+P+C)
    dC <- tau*P + (1-theta)*gama*I - d*C - r*C
    # Return the output of the model
    return(list(c(dS, dI, dP, dV, dC)))
  })
}

times <- seq(1, tf, by=Dt)

# Simulation ----
modSIRVsd <- ode(y = init.sd, times = times, func = SIRVsd, parms = par.SIRVsd, method = "ode45")
# matplot.0D(modSIRVsd, which = c("I", "P", "C"), lwd = 2)
# plot(modSIRVsd)
modSIRVsd <- as.data.frame(modSIRVsd)
modSIRVsd$N <- (modSIRVsd$S + modSIRVsd$I + modSIRVsd$P + modSIRVsd$V + modSIRVsd$C)
modSIRVsd$In <- modSIRVsd$I+modSIRVsd$P
modSIRVsd$pre <- modSIRVsd$In/modSIRVsd$N
modSIRVsd <- as.data.frame(modSIRVsd)
}
