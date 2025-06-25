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


select_best_k <- function(y, z, other_covariates = NULL,
                          k_range = 2:8,
                          np = 0, nq = 0, lamb = 0.5,
                          Ind = 0) {
  
  best_AIC <- Inf
  best_fit <- list()
  
  for (k in k_range) {
    # Build spline basis
    spline_basis <- ns(z, df = k)
    n_spline <- ncol(spline_basis)
    
    # Combine with other covariates
    if (!is.null(other_covariates)) {
      x <- cbind(spline_basis, other_covariates)
      n_other <- ncol(other_covariates)
    } else {
      x <- spline_basis
      n_other <- 0
    }
    
    nb <- 1 + n_spline + n_other  # intercept + spline + other covariates
    n_params <- nb + np + nq
    start_vals <- rep(0, n_params)
    
    # Fit model
    opt <- optim(par = start_vals,
                 fn = Like,
                 y = y,
                 nb = nb,
                 np = np,
                 nq = nq,
                 lamb = lamb,
                 x = x,
                 Ind = 0,
                 method = "BFGS",
                 hessian = FALSE)
    
    logLik <- -opt$value
    AIC <- 2 * n_params - 2 * logLik
    
    if (AIC < best_AIC) {
      best_AIC <- AIC
      best_k <- k
      best_coef <- opt$par
      best_nb <- nb
      
      # Extract and split coefficients
      beta0 <- best_coef[1]
      beta_spline <- best_coef[2:(1 + n_spline)]
      beta_other <- if (n_other > 0) best_coef[(2 + n_spline):(1 + n_spline + n_other)] else numeric(0)
      fi   <- if (np > 0) best_coef[(nb + 1):(nb + np)] else numeric(0)
      theta <- if (nq > 0) best_coef[(nb + np + 1):(nb + np + nq)] else numeric(0)
      
      best_fit <- list(
        best_k = best_k,
        best_AIC = best_AIC,
        best_logLik = logLik,
        beta0 = beta0,
        beta_spline = beta_spline,
        beta_other = beta_other,
        fi = fi,
        theta = theta
      )
    }
  }
  
  return(best_fit)
}


df <- read.csv("Falencia.csv",header = TRUE,sep=';')
df<- df %>% select(Falencia,INDPRO,PERMIT,BAA,GS10,FEDFUNDS)

df[] <- lapply(df, function(x) {
  if (is.character(x)) as.numeric(gsub(",", ".", x)) else x
})

df <- df %>%
  mutate(across(-1, as.numeric))

df <- df %>% select(-c(9,10,11))
df<- df %>% drop_na()

summary(df)

X <- as.matrix(df[ , setdiff(names(df), "Falencia")])

fit =glarma(df$Falencia, X=X[,c('BAA','GS10')] , phiLags = c(1,2,3) ,thetaLags = NULL, type = "Poi")
summary(fit)

# Example: plot residuals vs. covariate
plot(X[, "BAA"], df$Falencia, main = "Falencia x BAA")

other_covs <- as.matrix(df$GS10)

resposta<-select_best_k(y=df$Falencia,z=df$BAA,other_covariates = other_covs,
                        k_range=2:8,np=3,nq=0,lamb =0.5,Ind=0)
resposta
