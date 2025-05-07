## ALGORITHME ML MODELE 2 

# Fonction G
G_model2 <- function(v, p, k1, k2) {
  return(p * v^k1 + (1 - p) * v^k2)
}

# Inverse de G:
G_model2_inv <- function(y, p, k1, k2) {
  f <- function(v) { G_model2(v, p, k1, k2) - y }
  return(uniroot(f, c(0, 1), tol = 1e-10)$root) # Plus grande précision
}

# Inverse de H(x)
H_inv <- function(v, xi, sigma) {
  return(sigma / xi * (v^(-xi) - 1))
}

# Inverse de G2
G2_inv <- function(y, p, k1, k2, sigma, xi) {
  return(H_inv(G_model2_inv(y, p, k1, k2), xi, sigma))
}


# On vérifie la G2
x = seq(0, 7, by = 0.0001)
params = c(0.4, 2, 5, 1, 0.2)

p <- params[1]
k1 <- params[2]
k2 <- params[3]
sigma <- params[4]
xi <- params[5]


# Fonction de log-vraisemblance
log_g2 <- function(params, X) {
  p <- params[1]
  k1 <- params[2]
  k2 <- params[3]
  sigma <- params[4]
  xi <- params[5]
  
  if (xi <= 0 || p <= 0 || p >= 1 || sigma <= 0 || k1 <= 0 || k2 <= 0) {
    return(10^7)
  }
  
  n <- length(X)
  
  H_X <- pgpd(X, loc=0, scale=sigma, shape=xi)
  h_X <- dgpd(X, loc=0, scale=sigma, shape=xi)
  
  if (any(H_X <= 0 | H_X >= 1 | h_X <= 0)) {
    return(10^7)
  }
  
  g2_X <- (p * k1 * H_X^(k1-1) + (1 - p) * k2 * H_X^(k2-1)) * h_X
  
  if (any(g2_X <= 0)) {
    return(10^7)
  }
  
  log_vraisemblance <- sum(log(g2_X))
  
  return(-log_vraisemblance)
}


# BOXPLOT

est_model2 <- function(m=100, n=100, p=0.4, k1=2, k2=5, sigma=1, xi=0.2, n_start=5) {
  
  p_hat <- c()
  k1_hat <- c()
  k2_hat <- c()
  sigma_hat <- c()
  xi_hat <- c()
  
  for (i in 1:m){
    # Génération d'un échantillon aléatoire
    U <- runif(n = n, min = 0, max = 1)
    X2 <- sapply(U, function(u){G2_inv(u, p, k1, k2, sigma, xi)})
    
    best_value <- Inf
    best_par <- NULL
    
    for (start in 1:n_start) {
      # Conditions initiales aléatoires autour des vrais paramètres
      vecteur_initial <- c(
        runif(1, 0.3, 0.7),    
        runif(1, 1, 3),        
        runif(1, 4, 6),        
        runif(1, 0.8, 1.2),   
        runif(1, 0.1, 0.3)     
      )
      
      bornes_inf <- c(0.1, 0.5, 0.5, 0.1, 0.1)
      bornes_sup <- c(0.9, 10, 10, 10, 2)
      
      estimations_G2 <- optim(vecteur_initial, log_g2, method = "L-BFGS-B",
                              lower = bornes_inf, upper = bornes_sup, X = X2,
                              control = list(maxit = 1000))
      
      if (estimations_G2$convergence == 0 && estimations_G2$value < best_value) {
        best_value <- estimations_G2$value
        best_par <- estimations_G2$par
      }
    
      p_hat <- c(p_hat, best_par[1])
      k1_hat <- c(k1_hat, best_par[2])
      k2_hat <- c(k2_hat, best_par[3])
      sigma_hat <- c(sigma_hat, best_par[4])
      xi_hat <- c(xi_hat, best_par[5])

    }
  }
  
  return(list(p_hat=p_hat, k1_hat=k1_hat, k2_hat=k2_hat, sigma_hat=sigma_hat, xi_hat=xi_hat))
}



estimations_param_G2 <- est_model2(m=10^4, n=300, p=0.4, k1=2, k2=5, sigma=1, xi=0.2)

par(mfrow=c(1,5))

boxplot(estimations_param_G2$p_hat, main=bquote("(Modèle 2)" ~ p == .(p)), ylab=expression(hat(p)), col="lightblue")
abline(h=0.4, col="red", lwd=2)

boxplot(estimations_param_G2$k1_hat, main=bquote("(Modèle 2)" ~ k[1] == .(k1)), ylab=expression(hat(k)[1]), col="lightgreen")
abline(h=2, col="red", lwd=2)

boxplot(estimations_param_G2$k2_hat, main=bquote("(Modèle 2)" ~ k[2] == .(k2)), ylab=expression(hat(k)[2]), col="lightyellow")
abline(h=5, col="red", lwd=2)

boxplot(estimations_param_G2$sigma_hat, main=bquote("(Modèle 2)" ~ sigma == .(sigma)), ylab=expression(hat(sigma)), col="lightpink")
abline(h=1, col="red", lwd=2)

boxplot(estimations_param_G2$xi_hat, main=bquote("(Modèle 2)" ~ xi == .(xi)), ylab=expression(hat(xi)), col="lightcoral")
abline(h=0.2, col="red", lwd=2)
