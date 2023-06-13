

library(ggplot2); library(lubridate)

cbPalette <- c("red", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")


tiff(filename="Fig3.1.tif", width=260, height=180, pointsize=12,
     units="mm", res=400, compression="lzw")

ggplot(modSIRVsd, aes(floor_date(ymd(month), unit = "month")))+ 
  geom_line(aes(y=PI, colour="P.Infected"), size=1)+ 
  geom_line(aes(y=V, colour="Vaccined"), size=1)+
  geom_line(aes(y=S, colour="susceptible"), size=1)+
  geom_line(aes(y=I, colour="Infected"), size=1)+ 
  geom_line(aes(y=R, colour="Removed"), size=1)+ 
  geom_line(aes(y=N, colour="Total"), size=1)+ 
  ylab('Population')+
  xlab('Time (months)')+
  scale_color_manual(values=cbPalette)+
  labs(colour = "Simulation 
  dynamics")+
  theme_minimal()



ggplot(modSIRVsd)+ 
  geom_line(aes(ymd(month), I, colour="Infected"), size=1)+ 
  geom_line(aes(ymd(month), PI, colour="P Infected"), size=1)+
  geom_point(aes(ymd(month), cases), size=0.5, color="red")+
  geom_point(aes(ymd(month), cases*100/6), size=0.05, colour="black")+
  ylim(0, 42000)+
  ylab('Population')+
  xlab('Time (Months)')+
  labs(colour = "Compartment 
    dynamics")+
  theme_minimal()
