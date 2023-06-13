setwd("~/Dropbox/0.USP/7.Publicações/Mathematical modelling CSFV/Graphics python/1-beta-gamma-adjustment/")

# Reading CSV
# Observed
Fig24 <- read.csv2(file="Dados_2014_betaMJUN.csv", header=FALSE, sep=";")
names(Fig24) <- c("Tempo","Casos")

# Adjusted
adjusted <- read.csv(file = "t.csv")


# Ploting and saving the file
tiff(filename="beta_gamma14-17.tif", width=200, height=180, pointsize=12,
     units="mm", res=400, compression="lzw")

ggplot()+ 
  geom_point(data=Fig24, aes(x=Tempo, y=Casos), col= "BLACK", size= 1.5)+ 
  geom_line(data=adjusted, aes(x=tempo, y=imode, col ="Adjusted")) + 
  ylab('Number of cases 2014') +
  xlab('Time (month)') +
  theme(legend.title=element_blank())

dev.off()

