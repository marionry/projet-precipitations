library(evd) 
library(ggpubr)
install.packages("ggpubr")
seuil=0
sigma=0.5
xi=1
x <- seq(0,4,by=0.0001)
gp <- dgpd(x,loc=seuil,scale=sigma,shape=xi)

plot(x,gp,type='l',col="red",lty=1,ylab=expression("densité de "~h[xi](x)),main=expression("Densité de la GP selon le signe de"~xi ))
lines(x,dgpd(x,loc=seuil,scale=sigma,shape=-0.5),lty=2,col="purple")
lines(x,dgpd(x,loc=seuil,scale=sigma,shape=0),lty=3,linetype = "dotted",col="black")

legend("topright", legend = c(expression(xi > 0), expression(xi == 0), expression(xi < 0)), col = c("red", "black","purple"),lwd =2,lty= c(1,2,3))
