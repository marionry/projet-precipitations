#pour juste tracer la generalised pareto density

install.packages("evd")
library(evd)  

x<-seq(0,7, by=0.0001)
par(mfrow=c(1, 1))
?pgpd
sigma=1
xi=0.5
seuil=0

gp<-dgpd(x, loc = seuil, scale = sigma, shape = xi)

plot(x, gp, type = 'l', col = "black", lwd = 1, lty=2,
     main = "Densité modèle 1 comparé à la GP et la Gamma", 
     xlab = "x (mm)", ylab = "Densité de probabilité")


# add la gamma

gd <- dgamma(x, shape = 1.4, scale = 1.4)
lines(x, gd, col = "blue", lwd = 2, lty = 3) # Densité Gamma (bleu, continu)

#pour ajouter le modele i
G1<-function(k){
  return((k/sigma*(xi*x/sigma+1)^(-1/xi/sigma-1))*(1-(xi*x/sigma+1)^(-1/xi))^(k-1))
}

lines(x,G1(2),col="red",lwd=2,lty=4)
lines(x,G1(5),col="purple",lwd=2,lty=1)

legend("topright", legend = c("GP (k=1)","Gamma (1.4,1.4)", "Modèle 1 (k=2)", "Modèle 1 (k=5)"),
       col = c("black","blue", "red", "purple"), lwd = 2, lty = c(2, 3,5, 4, 1))


#pour ajouter le modèle 2

G2<-function(k1,k2){
  return ( (1/(2*(x/2+1)^3))*(k1*(1-1/((x/2+1)^2))^(k1-1)+(k2)*(1-1/(x/2+1)^2)^(k2-1)))}


plot(x, gd, type = 'l', col = "black", lwd = 1, lty=2,
     main = "Densité modèle 2 comparé à la Gamma avec p=1/2", 
     xlab = "x (mm)", ylab = "Densité de probabilité")

lines(x,G2(2,2),col="red",lwd=2,lty=4)
lines(x,G2(2,5),col="purple",lwd=2,lty=1)
lines(x,G2(2,10),col="pink",lwd=2,lty=3)

legend("topright", legend = c("Gamma (1.4,1.4)", "k1=2,k2=2", "k1=1,k2=5", "k1=2,k2=10"),
      col = c("black","red", "purple","pink"), lwd = 2, lty = c(2, 4,1,3))

#pour le modele 3


G3t<- function(d){
  return ( ((1+d)/d)*(1+xi*x)^(-1/xi -1)*(1-(1+xi*x)^(-d/xi)))
}
plot(x, gp, type = 'l', col = "black", lwd = 1, lty=2,
     main = "Densité modèle 3 comparé à la GP et la Gamma", 
     xlab = "x(mm)", ylab = "Densité de probabilité")
lines(x,G3t(1),col="red",lwd=2,lty=4)
lines(x,G3t(5),col="purple",lwd=2,lty=1)

lines(x, gd, col = "blue", lwd = 2, lty = 3)
legend("topright", legend = c("GP (d=infini)", "d=1", "d=5", "Gamma(1.4,1.4)"),
       col = c("black","red", "purple","blue"), lwd = 2, lty = c(2, 4,1,3))


#MODELE 4

G4<- function(d,k){
  c1<- (k/ (2*sigma))*((1+d)/d) *(1+xi*x/sigma)^(-1/xi -1)*(1-(1+xi*x/sigma)^(-d/xi))
  c2<- (1-((1+d)/d)*(1+xi*x/sigma)^(-1/xi) * (1-(1+xi*x/sigma)^(-d/xi)/(1+d)))^(k/2 -1)
  return (c1*c2)
}
plot(x, gp, type = 'l', col = "black", lwd = 1, lty=2,
     main = "Densité modèle 4 comparé à la GP et la Gamma", 
     xlab = "x(mm)", ylab = "Densité de probabilité")
lines(x,G4(5,2),col="red",lwd=2,lty=4)
lines(x,G4(2,2),col="pink",lwd=2,lty=5)
lines(x,G4(5,10),col="purple",lwd=2,lty=1)

lines(x, gd, col = "blue", lwd = 2, lty = 3)
legend("topright", legend = c("GP (d=infini)", "d=5,k=2","d=2,k=2", "d=5,k=10", "Gamma(1.4,1.4)"),
       col = c("black","red","pink", "purple","blue"), lwd = 2, lty = c(2, 4,5,1,3))

