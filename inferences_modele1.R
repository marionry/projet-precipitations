library(evd) 

#---------------MODELE 1-------------

# param à estimer
k <- 2
sigma <- 1
xi <- 0.2


# Fonction G1
G1 <- function(k, xi, x) {
  return((1 - (1 + xi * x)^(-1/xi))^k)
}

# Inverse de G1
G1_inv <- function(k, xi, y) {
  return(((1 - y^(1/k))^(-xi) - 1) / xi)
}

# Densité g1
g1 <- function(k, xi, x) {
  return((k * (xi * x + 1)^(-1/xi - 1)) * (1 - (xi * x + 1)^(-1/xi))^(k - 1))
}


# Log-vraisemblance
log_g1 <- function(params, X){
  xi <- params[1]
  sigma <- params[2]
  k <- params[3]
  
  # Contraintes sur les params 
  if (xi < 0|| sigma <= 0|| k<=0){return(1e7)}  # pour eviter xi<0
  
  n <- length(X)
  y <- n*log(k)+(k-1)*sum(log(1-(1+xi*X/sigma)^(-1/xi)))-n*log(sigma)-(1+1/xi)*sum(log(1+xi*X/sigma))
  
  return(-y)  # On maximise la log-vraisemblance (i.e minimise -log_v)
}



#BOXPLOT

est_model1 <- function(m=20,n=100,xi=0.2,sigma=1,k=2,X){
  xi_hat <- c()
  sigma_hat <- c()
  k_hat <- c()
  
  for (i in 1:m){
    # Génération d'un échantillon
    U <- runif(n = 100, min = 0, max = 1)
    X <- G1_inv(k,xi,U)
    
    #estimations
    vecteur_initial <- c(0.1, 0.5, 1.5)  # Valeurs initiales cohérentes
    estimations_G1 <- optim(vecteur_initial,log_g1,X=X)
    
    xi_hat_G1 <- estimations_G1$par[1]
    sigma_hat_G1 <- estimations_G1$par[2]
    k_hat_G1 <- estimations_G1$par[3]
    
    #Implementation de la liste des estimations
    xi_hat <- c(xi_hat,xi_hat_G1)
    sigma_hat <- c(sigma_hat,sigma_hat_G1)
    k_hat <- c(k_hat,k_hat_G1)
    
    
    
    
  }
  return(list(xi_hat=xi_hat,sigma_hat=sigma_hat,k_hat=k_hat))
}

estimations_param_G1 <- est_model1(m=10^4,n=100,xi=0.2,sigma=1,k=2,X)

#On trace les boxplots

# Définir l'espace graphique pour 3 sous-plots (1 ligne, 3 colonnes)
par(mfrow=c(1, 3))

# Boxplot pour les valeurs de xi_hat
boxplot(estimations_param_G1$xi_hat, main=bquote("(Modèle 1)" ~ xi == .(xi)),ylab=expression(xi ~ "_hat"), col="lightblue")
abline(h=xi, col="red", lwd=2) # Ajouter une ligne rouge à xi=0.2

# Boxplot pour les valeurs de k_hat
boxplot(estimations_param_G1$k_hat, main=bquote("(Modèle 1)" ~ k == .(k)), ylab=expression(k ~ "_hat"), col="lightgreen")
abline(h=k, col="red", lwd=2)

# Boxplot pour les valeurs de sigma_hat
boxplot(estimations_param_G1$sigma_hat, main=bquote("(Modèle 1)" ~ sigma == .(sigma)),ylab=expression(sigma ~ "_hat"), col="lightcoral")
abline(h=sigma, col="red", lwd=2)
