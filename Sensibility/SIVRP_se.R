SIRVP_se <- function (beta, gama, tau, mu, w, theta, tf)
  {
# Parameters----
#Montly parameters ----
# beta <- 4.0706 #Calculada nls
# gama <- 2.7941 #Calculada nls
# tau <- gama/30 #
# mu <- (1/(350/30))
# w <- 1/(180/30) #omega
# theta <- 0.1

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
fp <- function(t){
  {fp <- c[t]}
  return(fp)
}

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
vac <-(sus)*(0.75) #initial vaccinated pop as jan 2022

sus <- pop-vac-r
vac+sus+r-pop

#  Set intial vector
init.sd <- c(S=sus, I=infe, PI=p_infe, V=vac, R=r)                      

# Deterministic modell  ----
# Calling function
# Time interval 2014-2022  ----
Dt <- 1
tf <- 12*2   # Final time 2022-2033

# SIRV with population dynamics  ----
SIRVsd <- function(t, state, parameters){
  with(as.list(c(state, parameters)),{
    
    # Rate of change
    dV = (fp(t))*mu*(S+I+PI+V+R) - w*V - mu*V
    dS = (1-(fp(t)))*mu*(S+I+PI+V+R) - (beta*S*I)/(S+I+PI+V) + w*V - mu*S
    dI = (beta*S*I)/(S+I+PI+V+R) - gama*I*(1-theta) - gama*I*theta - mu*I
    dPI= gama*I*theta - tau*PI - mu*PI
    dR= gama*I*(1-theta) + tau*PI - mu*R
    # Return the output of the model
    return(list(c(dS, dI, dPI, dV, dR)))
  })
}

times <- seq(1, tf, by=Dt)

# Simulation ----
modSIRVsd <- ode(y = init.sd, times = times, func = SIRVsd, parms = par.SIRVsd, method = "ode45")
modSIRVsd <- as.data.frame(modSIRVsd)
modSIRVsd$N <- (modSIRVsd$S + modSIRVsd$I + modSIRVsd$PI + modSIRVsd$V + modSIRVsd$R)
modSIRVsd$p <- modSIRVsd$I+modSIRVsd$PI
modSIRVsd <- as.data.frame(modSIRVsd)
}