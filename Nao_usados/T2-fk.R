library(pacman)
p_load(ggplot2,forecast,dlm,numDeriv,plotly)
setwd('C:\\Users\\Marcelo\\OneDrive\\Área de Trabalho\\ts\\')
source('funcoes_auxiliares.R')
getwd()


##### StructTS

#####1)

# Fit the local level model

y_mnl <- simul_y_mnl(T=100,10,0.5,n_seed=1)

fit <- StructTS(y_mnl, "level")

# Extract the estimated level
df <- data.frame(
  time = 1:length(y_mnl),
  observed = y_mnl,
  estimated_level = fitted(fit)[, "level"]
)

# Extract estimated variances
sigma2_epsilon <- fit$coef["epsilon"]
sigma2_eta <- fit$coef["level"]


ggplot(df, aes(x = time)) +
  geom_line(aes(y = observed, color = "Observed y_mnl"), size = 1) +
  geom_line(aes(y = estimated_level, color = "Estimated Level"), size = 1) +
  scale_color_manual(values = c("Observed y_mnl" = "black", "Estimated Level" = "red")) +
  labs(
    title = "Observed Data vs. Estimated Level",
    x = "Time",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal()

# Print estimated variances
cat("Estimated Observation Noise Variance (sigma^2_epsilon):", sigma2_epsilon, "\n")
cat("Estimated State Noise Variance (sigma^2_eta):", sigma2_eta, "\n")

######### simulação com T = 10000

y_mnl <- simul_y_mnl(T=1000,10,0.5,1)

fit <- StructTS(y_mnl, "level")

# Extract the estimated level
df <- data.frame(
  time = 1:length(y_mnl),
  observed = y_mnl,
  estimated_level = fitted(fit)[, "level"]
)

# Extract estimated variances
sigma2_epsilon <- fit$coef["epsilon"]
sigma2_eta <- fit$coef["level"]

ggplot(df, aes(x = time)) +
  geom_line(aes(y = observed, color = "Observed y_mnl"), size = 1) +
  geom_line(aes(y = estimated_level, color = "Estimated Level"), size = 1,alpha=0.3) +
  scale_color_manual(values = c("Observed y_mnl" = "black", "Estimated Level" = "red")) +
  labs(
    title = "Observed Data vs. Estimated Level",
    x = "Time",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal()

# Print estimated variances
cat("Estimated Observation Noise Variance (sigma^2_epsilon):", sigma2_epsilon, "\n")
cat("Estimated State Noise Variance (sigma^2_eta):", sigma2_eta, "\n")

# O estimador de MLE é consistente, 
#então para T alto, o estimador MLe das variância é próximo ao verdadeiro

# recuperando manualmente a séria do nível ajsutada:

fitted_mu_manual<-mnl_fk(T=1000,sigma2_epsilon,sigma2_eta,a0=0,p0=1)

df_fitted_mnl <- data.frame(manual = fitted_mu_manual,structts=df$estimated_level)
tail(df_fitted_mnl)
head(df_fitted_mnl)


############ MTL structts

T <- 100

sigma2_epsilon <- 1   # Observation noise variance
sigma2_eta <- 1   # State noise variance for the level
sigma2_qsi <- 0.1    # State noise variance for the slope
n_seed <- 4          # Seed for reproducibility

y_mtl <- simul_y_mtl(T, sigma2_eta=sigma2_eta, sigma2_qsi=sigma2_qsi, sigma2_epsilon=sigma2_epsilon, n_seed=42)

fit_mtl <- StructTS(y_mtl, "trend")

df_mtl <- data.frame(
  time = 1:T,
  observed = y_mtl,
  adjusted_level = fit_mtl$fitted[, "level"],
  adjusted_trend = fit_mtl$fitted[, "slope"]
)

# Plotando os valores observados e ajustados

ggplot(df_mtl, aes(x = time)) +
  geom_line(aes(y = observed, color = "Observed"), size = 1) +
  geom_line(aes(y = adjusted_level, color = "Adjusted Level"), size = 1,alpha=0.3) +
  scale_color_manual(values = c("Observed" = "black", "Adjusted Level" = "red")) +
  labs(
    title = "Observed vs Adjusted Level StructTS (MTL)",
    x = "Time",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal()

ggplot(df_mtl, aes(x = time)) +
  geom_line(aes(y = adjusted_trend, color = "Adjusted Trend"), size = 1) +
  scale_color_manual(values = c("Adjusted Trend" = "black")) +
  labs(
    title = "Trend using StructTS (MTL)",
    x = "Time",
    y = "Value",
    color = "Legend"
  ) +
  theme_minimal()

cat("Estimated Observation Noise Variance (sigma^2_epsilon):", fit_mtl$coef["epsilon"], "\n")
cat("Estimated State Noise Variance (sigma^2_eta):", fit_mtl$coef["level"], "\n")
cat("Estimated State Noise Variance (sigma^2_qsi):", fit_mtl$coef["slope"], "\n")


# Previsão e suavização no StructTs:

y_mtl <- simul_y_mtl(T=100,10,0.1,1,n_seed=1)

# Ajustando o modelo com StructTS
fit <- StructTS(y_mtl, type = "trend")

# Obtendo previsões

h <- 20  # Horizontes futuros

forecast_obj <- forecast(fit, h = h, level = 90)

# Extraindo estados suavizados
smoothed <- tsSmooth(fit)

# Criando dataframe para o ggplot

df <- data.frame(
  Tempo = 1:T,
  Observado = y_mtl,
  Suavizado = smoothed[, 1]
)

df_pred <- data.frame(
  Tempo = (T + 1):(T + h),
  Previsao = forecast_obj$mean,
  Lower = forecast_obj$lower[, 1],
  Upper = forecast_obj$upper[, 1]
)

# Plotando
ggplot() +
  geom_line(data = df, aes(x = Tempo, y = Observado), color = "black", size = 1) +
  geom_line(data = df, aes(x = Tempo, y = Suavizado), color = "blue", size = 1) +
  geom_line(data = df_pred, aes(x = Tempo, y = Previsao), color = "red", size = 1) +
  geom_ribbon(data = df_pred, aes(x = Tempo, ymin = Lower, ymax = Upper), fill = "red", alpha = 0.2) +
  labs(title = "Modelo de Nível Local: Observado, Suavizado e Previsão",
       x = "Tempo", y = "Valor") +
  theme_minimal()
   
# exercícios: 
#1) Modelo MNL - gerar 50 com T=100 e 50 com T = 1000 com sigma2_epsilon=1 e sigma2_eta=0.5 
#2) estimar o modelo para cada amostra
#3) comparar o boxplot da estimativa sigma2_epsilon para os 2 tamanhos de amostra.

################# DLM
# parm Initial values for the parameters to be estimated.
# DLM (a local level model is a polynomial DLM of order 1, a local linear trend
# is a polynomial DLM of order 2)
# 
#  The arguments dV and dW are used to specify the diagonal of the observation and evolution
# covariance matrices respectively

T=100
y_mnl <- simul_y_mnl(T=T,1,0.5,1)

#  The function dlmModPoly in the dlm package in R creates
# a state-space model for a polynomial trend process. It helps define local level (MNL) and local trend (MTL) models 
# by setting up the system matricesrequired for the Kalman filter.

build_mnl <- function(theta) {
    dlmModPoly(order = 1, dV = theta[1], dW = theta[2])
}


# parm Initial values for the parameters to be estimated.
# 
# lower <- constraint in param values


fit <- dlmMLE(y_mnl, parm = c(100, 2), build_mnl, lower = rep(1e-4, 2))

mod_mnl <- build_mnl(fit$par)
drop(V(mod_mnl))
drop(W(mod_mnl))

# The inverse of the Hessian matrix of the negative loglikelihood function evaluated at the
# MLEs is, by standard maximum likelihood theory, an estimate of the asymptotic variance
# matrix of the maximum likelihood estimators

#obtendo a hessiana
hs <- hessian(function(x) dlmLL(y_mnl, build_mnl(x)), fit$par)
#Checando se é positiva definida
all(eigen(hs, only.values = TRUE)$values > 0)


aVar <- solve(hs)

cat("Standard Error of dV:", sqrt(diag(aVar)[1]), "\n")
cat("Standard Error of dW:", sqrt(diag(aVar)[2]), "\n")

# Apply Kalman & smoothing
smooth_mnl <- dlmSmooth(y_mnl, mod_mnl)
filter_mnl <- dlmFilter(y_mnl,mod_mnl)
# Extract smoothed state estimates
mu_hat <- drop(smooth_mnl$s)

# Extracting standard deviations (square roots of variances) from the covariance matrix
std_dev <- sqrt(diag(aVar))

# Forecast using dlmForecast

# Gerando as previsões e o erro padrão com dlmForecast
fore_mnl <- dlmForecast(filter_mnl, nAhead = 10)

# Extraímos os resultados das previsões
#matriz de valores esperados para futuras obs
f <- fore_mnl$f
#lista de variâncias de futuras observações
Q <- fore_mnl$Q



# Calculando o erro padrão das previsões (banda de 50%)
hwidth <- qnorm(0.25, lower = FALSE) * sqrt(unlist(Q))  # 50% intervalo de confiança

# Criando o intervalo de previsão
lower <- f - hwidth
upper <- f + hwidth

# Criar uma data frame para facilitar o uso no ggplot2
df_forecast <- data.frame(
  Time = (length(y_mnl) + 1):(length(y_mnl) + 10),
  Forecasted = f,
  Lower = lower,
  Upper = upper
)

# Criando o gráfico
# Criando o gráfico
ggplot(df_forecast, aes(x = Time)) +
  # Adiciona as previsões como uma linha vermelha
  geom_line(aes(y = Forecasted, color = "Forecasted"), size = 1) +
  
  # Adiciona as bandas de previsão (intervalo de 50%)
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "red", alpha = 0.3) +
  
  # Adiciona as observações reais da série (y_mnl)
  geom_line(data = data.frame(Time = 1:length(y_mnl), Observed = y_mnl),
            aes(y = Observed, color = "Observed"), size = 1, linetype = "solid") +
  
  # Adiciona a linha de suavização (suavizado)
  geom_line(data = data.frame(Time = time(smooth_mnl$s), Smoothed = smooth_mnl$s),
            aes(y = Smoothed, color = "Smoothed"), size = 1, linetype = "dashed") +
  
  # Personalizando os eixos e título
  labs(
    title = "Forecasted Values with Prediction Intervals",
    x = "Time",
    y = "Level",
    color = "Legend"
  ) +
  
  # Customizando a paleta de cores
  scale_color_manual(values = c("Forecasted" = "red", "Observed" = "black", "Smoothed" = "blue")) +
  
  # Tema minimalista
  theme_minimal() +
  theme(legend.position = "top")

df_residual <- data.frame(resid = residuals(filter_mnl, sd = FALSE),zero=rep(0,T),tempo<-seq(1:T))

ggplot() +
  geom_line(data = df_residual, aes(x = tempo, y = resid), color = "black", size = 1) +
  geom_line(data=df_residual, aes(x=tempo,y=zero))+
  geom_point(data=df_residual,aes(x = tempo , y = resid), color = "black", size = 2 ,alpha=0.3) +
  
  labs(title = "Modelo de Nível Local: Erro de previsão",
       x = "Tempo", y = "resid") +
  theme_minimal()                        

shapiro.test(df_residual$resid)
cat("media_residuos:",mean(df_residual$resid))
cat("var_residuos:",var(df_residual$resid))


T=1000
y_mnl <- simul_y_mnl(T=T,1,0.5,1)
fit <- dlmMLE(y_mnl, parm = c(100, 2), build_mnl, lower = rep(1e-4, 2))

# Define search space
sigma_epsilon2_vals <- seq(0.1,4, length.out = 100)
sigma_eta2_vals <- seq(0.1, 4, length.out = 100)

# Compute log-likelihood over the grid
likelihood_matrix <- outer(sigma_epsilon2_vals, sigma_eta2_vals, 
                           Vectorize(function(se, sw) logLik_mnl(c(se, sw))))

# Convert to dataframe for ggplot

likelihood_df <- expand.grid(sigma_epsilon2 = sigma_epsilon2_vals, 
                             sigma_eta2 = sigma_eta2_vals)

likelihood_df$logLik <- as.vector(likelihood_matrix)

# Plot the likelihood surface
ggplot(likelihood_df, aes(x = sigma_epsilon2, y = sigma_eta2, fill = logLik)) +
  geom_tile() +
  scale_fill_viridis_c()+
  labs(title = "Log-Likelihood Surface for MNL Model",
       x = expression(sigma[epsilon]^2),
       y = expression(sigma[eta]^2),
       fill = "Log-Likelihood") +
  theme_minimal()

# Compare with dlmMLE results
cat("Estimated sigma^2_epsilon:", fit$par[1], "\n")
cat("Estimated sigma^2_eta:", fit$par[2], "\n")

mle_estimates <- likelihood_df %>%
  filter(logLik == max(logLik)) %>%
  select(sigma_epsilon2, sigma_eta2, logLik)

names(mle_estimates)<-c("sigma2_epsilon",'sigma2_eta','Max LogLikelihood')

print(mle_estimates)

fig <- plot_ly(
  x = unique(likelihood_df$sigma_epsilon2),  # X-axis (sigma_epsilon2 values)
  y = unique(likelihood_df$sigma_eta2),      # Y-axis (sigma_eta2 values)
  z = matrix(likelihood_df$logLik, 
             nrow = length(unique(likelihood_df$sigma_epsilon2)), 
             ncol = length(unique(likelihood_df$sigma_eta2))),  # Reshape log-likelihood data
  type = "surface",
  colorscale = "Inferno"  # Makes high values vivid red
)

fig <- fig %>%
  layout(
    title = "Log-Likelihood Surface for MNL Model",
    scene = list(
      xaxis = list(title = expression(sigma[epsilon]^2)),
      yaxis = list(title = expression(sigma[eta]^2)),
      zaxis = list(title = "Log-Likelihood")
    )
  )


fig  # Display the plo

  
  