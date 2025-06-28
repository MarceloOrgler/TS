library(pacman)
p_load(dplyr,tidyr,splines,glarma)

setwd('C:\\Users\\Marcelo\\OneDrive\\Área de Trabalho\\ts\\trabalho 2 - TS\\Trabalho7')
getwd()

Like <- function(Parametros, y, nb, np, nq, lamb, x, Ind) {
  
  ny <- length(y)
  
  # Separating parameters: beta, fi, theta
  if (nb > 0) {
    beta <- Parametros[1:nb]
  }
  
  ar_ma_start <- nb + 1
  
  if (np > 0 && nq == 0) {
    fi <- Parametros[ar_ma_start:(ar_ma_start + np - 1)]
  }
  
  if (nq > 0 && np == 0) {
    theta <- Parametros[ar_ma_start:(ar_ma_start + nq - 1)]
  }
  
  if (np > 0 && nq > 0) {
    fi <- Parametros[ar_ma_start:(ar_ma_start + np - 1)]
    theta <- Parametros[(ar_ma_start + np):(ar_ma_start + np + nq - 1)]
  }
  
  # Initialize variables
  z    <- rep(0, ny)
  w    <- rep(0, ny)
  mi   <- rep(0, ny)
  erro <- rep(0, ny)
  
  # First observation: compute w[1]
  if (nb == 0) {
    w[1] <- 0
  } else if (nb == 1) {
    w[1] <- beta[1]
  } else {
    w[1] <- beta[1] + as.numeric(x[1, 1:(nb - 1)] %*% beta[2:nb])
  }
  
  mi[1]   <- exp(w[1])
  erro[1] <- (y[1] - mi[1]) / (mi[1]^lamb)
  L       <- y[1] * w[1] - exp(w[1])
  
  # Loop over t = 2 to ny
  for (t in 2:ny) {
    
    # AR component
    if (np > 0) {
      for (k in 1:min(np, t - 1)) {
        z[t] <- z[t] + fi[k] * (z[t - k] + erro[t - k])
      }
    }
    
    # MA component
    if (nq > 0) {
      for (k in 1:min(nq, t - 1)) {
        z[t] <- z[t] + theta[k] * erro[t - k]
      }
    }
    
    # Compute w[t] based on z[t] and x[t,]
    if (nb == 0) {
      w[t] <- z[t]
    } else if (nb == 1) {
      w[t] <- beta[1] + z[t]
    } else {
      w[t] <- beta[1] + as.numeric(x[t, 1:(nb - 1)] %*% beta[2:nb]) + z[t]
    }
    
    mi[t]   <- exp(w[t])
    erro[t] <- (y[t] - mi[t]) / (mi[t]^lamb)
    L       <- L + y[t] * w[t] - exp(w[t])
  }
  
  L <- -L  # negative log-likelihood
  
  valores <- cbind(erro, z, mi)
  
  if (Ind == 0) return(L)
  if (Ind == 1) return(valores)
}