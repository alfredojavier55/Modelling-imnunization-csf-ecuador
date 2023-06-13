################ AJUSTE MENSUAL #####################################

# dados de ajuste deslocados a jun 2014
# setwd("~/Dropbox/0.USP/4.Publicações/Modelling CSFV/Analise")
setwd("~/Dropbox/0.USP/7.Publicações/Mathematical modelling CSFV/Analise")
getwd()

# Carregando deSolve
library(deSolve)

# Carregando os dados
Fig24 <- read.csv2(file="Dados_2014_betaMJUN.csv", header=FALSE, sep=";")

names(Fig24) <- c("Tempo","Casos")

plot(Fig24, ylim=c(0,310), pch=1,
     xlab = "time (month)",
     ylab = "number of cases")

#Number of cases in 2014-2017
plot(Fig24, ylim=c(0,310), pch=1,
     xlab = "time (month)",
     ylab = "number of cases")

source("ImodSIsd.R")

# Condicao inicial das variaveis

N <- 3718
state.SIsd <- c(s=3709/3718, i=9/3718)

3718-9

# Tempo de simulacao
tsim <- 12
Dt <- 1

# Ajuste dos parametros
aj.SIsd <- nls(Casos ~ ImodSIsd(Tempo, state.SIsd, N, tsim, Dt, beta, gama),
                data = Fig24, algorithm="port",
                start = list(beta=1, gama=1), # chute inicial
                # lower = list(beta=0, gama=0),
                # upper = list(beta=100, gama =100),
                trace = TRUE)

# Resumo do ajuste
summary(aj.SIsd)

# Parametros ajustados
coef(aj.SIsd)

# Infectados estimados pelo modelo
Dt <- 0.01 
tempos.aju2 <- seq(0, tsim, by=Dt)
Imodelo2 <- ImodSIsd(tempos.aju2, state.SIsd, N = 3718, tsim, Dt, beta = 4.0706, gama = 2.7941)

Imodelo2

# Acrescentando os resultados do modelo no grafico dos dados observados
lines(tempos.aju2, Imodelo2, lty=3)
legend("topright", c("Observed prevalence", "Estimated"), lwd=1, pch=1, lty=c(0,3)) 
title(main = "Beta adjustment CSF Ecuador 2014 Ecuador",
      sub="Beta= 4.0706 d gamma= 2.7941")
par(ps = 10, cex = 1, cex.main = 1)


#Exportando a TIFF

tiff(filename="beta-gamma2014_paper1.tif", width=100, height=90, pointsize=10,
     units="mm", res=800, compression="lzw")

plot(Fig24, ylim=c(0,310), pch=1,
     xlab = "Time (month)",
     ylab = "Number of cases 2014")
lines(tempos.aju2, Imodelo2, lty=3)

dev.off()

lines(tempos.aju2, Imodelo2, lty=3)
legend("topright", c("Observed prevalence", "Estimated"), lwd=1, pch=1, lty=c(0,3)) 
title(main = "Beta adjustment CSF Ecuador 2014 Ecuador",
      sub="Beta= 4.0706 d gamma= 2.7941")
par(ps = 10, cex = 1, cex.main = 1)
grid(NULL,NULL)

dev.off()

library(ggplot2)

t <- data.frame(tempo=tempos.aju2, imode=Imodelo2)


tiff(filename="beta_gamma2014_adjust.tif", width=200, height=180, pointsize=12,
     units="mm", res=400, compression="lzw")

ggplot()+ 
  geom_point(data=Fig24, aes(x=Tempo, y=Casos), col= "BLACK", size= 1.5)+ 
  geom_line(data=t, aes(x=tempo, y=imode, col ="Stimated")) + 
  ylab('Number of cases 2014')+
  xlab('Time (month)')

dev.off()

write.csv(t, file = "t.csv")