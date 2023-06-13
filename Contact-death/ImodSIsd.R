# Alfred E bom sempre separar os codigos para nao ficar confuso

# Funcao ImodSIRsd.R
# Solucao numerica do modelo SIR sem demografia
# Retorna o numero de infectados

ImodSIsd <- function(Tempo, state.SIsd, N, tsim, Dt, beta, gama){
  
  par.SIsd <- c(beta=beta, gama=gama)
  
  # Modelo SIR sem demografia
  SIsd <- function(t, state, parameters){
    with(as.list(c(state, parameters)),{
      # rate of change
      ds <- -beta*s*i
      di <- beta*s*i - gama*i
      
      # return the output of the model
      return(list(c(ds, di)))
    })
  }
  
  times <- seq(0, tsim, by=Dt)
  
  modSIsd <- ode(y = state.SIsd, times = times, func = SIsd, parms = par.SIsd, method = "ode45")
  
  modSIsd <- as.data.frame(modSIsd)
  
  # O modelo retorna a fracao de infectados (i), mas queremos fazer o ajuste para o numero
  # Multiplicamos i pelo tamanho da populacao
  Ih = N*modSIsd$i
  
  # Fazendo a interpolacao do numero de infectados para um dado Tempo
  # No caso deste exemplo, nao seria necessario interpolar, pois os valores ja estão saindo
  # com os valores de tempo correspondentes aos dados
  I = approx(x = modSIsd$t, y = Ih, xout = Tempo)$y;
  
  return(I);
}

# Se coloca ao final o cifrão ($) ipsylon ao final para que ele retorne so o valor de y que seriam os infectados
?approx 
