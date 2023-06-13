setwd("/home/alfredo/Dropbox/0.USP/7.Publicações/Mathematical modelling CSFV/Graphics python/sensibility/")

# Monotonicity

write.csv(dat_1, file = "dat1.csv")
write.csv(dat_2, file = "dat2.csv")
write.csv(dat_3, file = "dat3.csv")
write.csv(dat_4, file = "dat4.csv")
write.csv(dat_5, file = "dat5.csv")
write.csv(dat_6, file = "dat6.csv")

dat_1 <- read_csv(file = "dat1.csv")
dat_2 <- read_csv(file = "dat2.csv")
dat_3 <- read_csv(file = "dat3.csv")
dat_4 <- read_csv(file = "dat4.csv")
dat_5 <- read_csv(file = "dat5.csv")
dat_6 <- read_csv(file = "dat6.csv")

par(mfrow = c(2, 3))

plot(dat_1$X1, dat_1$Y, xlab = "beta", ylab = "Cummulative cases", cex = .1)
plot(dat_2$X1, dat_2$Y, xlab = "gama", ylab = "Cummulative cases", cex = .1)
plot(dat_3$X1, dat_3$Y, xlab = "tau", ylab = "Cummulative cases", cex = .1)
plot(dat_4$X1, dat_4$Y, xlab = "mu", ylab = "Cummulative cases", cex = .1)
plot(dat_5$X1, dat_5$Y, xlab = "w", ylab = "Cummulative cases", cex = .1)
plot(dat_6$X1, dat_6$Y, xlab = "theta", ylab = "Cummulative cases", cex = .1)

par(mfrow = c(1, 1))

# PRCC
write.csv(df, file = "df.csv")
df <- read.csv(file = "df.csv")
df <- df[c(1,2,5), ]

ggplot(df) +
  geom_col(aes(x = reorder(Parameter, est), y = est), fill = "gray60") +
  geom_text(aes(x=reorder(Parameter, est), y=0.9*est, label=round(est,5))) +
  theme_minimal() +
  xlab(NULL) +
  ylab("PRCC") +
  theme(text = element_text(size = 14))
