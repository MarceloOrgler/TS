library(pacman)
p_load(ggplot2,dplyr,lmtest,forecast,dlm)

rm(list=ls())
#ajustando a área de trabalho:

setwd("C:\\Users\\Marcelo\\OneDrive\\Área de Trabalho\\ts\\trabalho5\\")
getwd()

#1)a) série temporalde produção de energia,	jan/1985  até jan/2018. Vinda de kaggle. 
data<- read.csv("Electric_Production.csv")
Y<- data$IPG2211A2N

#1)b)
H=12
n<-length(Y)
Y_train<-ts(Y[1:(n-H)],start=c(1985,1),frequency=12)  
Y_test <- Y[(n-H+1):n]

#1)c)
x11(width = 10, height = 12)  # Ou use `windows()` no Windows
plot(Y_train)

#Série não estacionária com tendência e sazonalidade.

#2)a)
par(mfrow=c(2,1),mar = c(4, 4, 2, 1))
acf(Y_train, lag = 40)
pacf(Y_train, lag = 40)

#série não estacionária
log_Y<-log(Y_train)
YDif<-diff(log_Y)
x11(width = 10, height = 12)
plot(YDif)

x11(width = 10, height = 12)  
par(mfrow = c(4, 2))
acf(YDif,lag = 40)
pacf(YDif, lag = 40)

YDif<-diff(log_Y,lag=12)
x11(width = 10, height = 12)
plot(YDif)

x11(width = 10, height = 12)  # Ou use `windows()` no Windows
par(mfrow = c(4, 2))
acf(YDif,lag = 40)
pacf(YDif,, lag = 40)

#2)B)

M1<-arima(log_Y, order = c(1, 1, 1),   seasonal = list(order = c(0, 1, 1)))   
coeftest(M1)
AIC(M1)

M3<-arima(log_Y, order = c(1, 1 , 1),   seasonal = list(order = c(1, 1, 1)))   
coeftest(M3)
AIC(M3)

M4<-arima(log_Y, order = c(1, 1, 1),   seasonal = list(order = c(0, 1, 2)))   
M4
coeftest(M4)
AIC(M4)

M5<-arima(log_Y, order = c(2, 1, 1),   seasonal = list(order = c(1, 1, 1)))   
coeftest(M5)
AIC(M5)

M6<-arima(log_Y, order = c(2, 1, 1),   seasonal = list(order = c(0, 1, 1)))   
coeftest(M6)
AIC(M6)

M7<-arima(log_Y, order = c(2, 2, 1),   seasonal = list(order = c(0, 0, 1)))   
coeftest(M7)
AIC(M7)

M8<-arima(log_Y, order = c(2, 2, 1),   seasonal = list(order = c(0, 1, 1)))   
coeftest(M8)
AIC(M8)

M9<-arima(log_Y, order = c(1, 2, 2),   seasonal = list(order = c(0, 1, 1)))   
coeftest(M9)
AIC(M9)

M10<-arima(log_Y, order = c(1, 2, 2),   seasonal = list(order = c(0, 1, 1)))   
coeftest(M10)
AIC(M10)

M_ARIMA<-M10


#2)C)
# Graficos de 
x11(width = 10, height = 12)
plot(as.numeric(M_ARIMA$res),type='l')

x11(width = 10, height = 12)
hist(M_ARIMA$res)

std_resid <- (M_ARIMA$res - mean(M_ARIMA$res))/sd(M_ARIMA$res)

x11(width = 10, height = 12)
qqnorm(std_resid, 
       ylab="Residuos", 
       xlab="Escores Normais", 
       main="QQ-Plot para os Resíduos do modelo M_ARIMA") 
qqline(std_resid)

par(mfrow=c(2,1),mar = c(4, 4, 2, 1))
acf(M_ARIMA$res)
pacf(M_ARIMA$res)

# Teste de normalidade para os residuos
shapiro.test(M_ARIMA$res)

# Teste Box-pierce até a ordem 12
Box.test(M_ARIMA$resid, lag = 12, type = c("Box-Pierce", "Ljung-Box"))

#2)d)

buildmod<-function(par){
    trend <- dlmModPoly(order=2,dV=exp(par[1]),dW=exp(par[2:3]))
    seasonal <- dlmModSeas(frequency=12,dV=0,dW=c(exp(par[4]),rep(0,10)))
    model<- trend+seasonal
    return(model)
}

fit<-dlmMLE(Y_train,parm=rep(0,4),build=buildmod)
fit

var_estimates <- exp(fit$par)
var_estimates

smoothed <- dlmSmooth(dlmFilter(Y_train, buildmod(exp(fit$par))))

plot_data <- data.frame(
  Data = time(Y_train),
  Observado = as.numeric(Y_train),
  Suavizado = smoothed$s[-1,1] + smoothed$s[-1,3],  # Nível + sazonalidade
  Tendencia = smoothed$s[-1,1],
  Sazonalidade = smoothed$s[-1,3]
)

ggplot(plot_data, aes(x = Data)) +
  geom_line(aes(y = Observado, color = "Observado"), linewidth = 0.7) +
  geom_line(aes(y = Suavizado, color = "Suavizado"), linewidth = 0.5) +
  labs(title = "Série Temporal: Observado vs Suavizado",
       y = "Valor",
       x = "Tempo") +
  scale_color_manual(name = "Série",
                     values = c("Observado" = "Black", 
                                "Suavizado" = "blue")) +
  theme_minimal()

#3)

Prev=forecast(M_ARIMA, H, level=c(95))

time_index <- time(Y)

df <- data.frame(
  Time = time_index,
  Observed = as.numeric(Y),
  Test = c(rep(NA, n-H), Y_test),
  Fitted = c(as.numeric(exp(fitted(M_ARIMA))),rep(NA,H)),
  Forecast = rep(NA, n),
  Lower = rep(NA, n),
  Upper = rep(NA, n)
)

# Fill in forecast values
for(i in 1:H) {
  df$Forecast[n-H+i] <- exp(Prev$mean[i])
  df$Lower[n-H+i] <- exp(Prev$lower[i])
  df$Upper[n-H+i] <- exp(Prev$upper[i])
}

# Create the plot
ggplot(df, aes(x = Time)) +
  geom_line(aes(y = Observed, color = "Observed")) +
  geom_line(aes(y = Test, color = "Test"), na.rm = TRUE) +
  geom_line(aes(y = Fitted, color = "Fitted"), na.rm = TRUE) +
  geom_line(aes(y = Forecast, color = "Forecast"), na.rm = TRUE) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), 
              fill = "red", alpha = 0.2, na.rm = TRUE) +
  scale_color_manual(values = c("Observed" = "black", 
                                "Test" = "green",
                                "Fitted" = "blue",
                                "Forecast" = "red")) +
  labs(x = "Year", y = "Y", color = "Series") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Calculo do EQMP
EQMP <- sum((as.numeric(Y_test)-exp(Prev$mean))^2)/H
EQMP

fit<-dlmFilter(Y_train,buildmod(exp(fit$par)))
forecast<-dlmForecast(fit,nAhead=12)

q_value <- qnorm(0.975)
lower <- forecast$f - q_value * sqrt(unlist(forecast$Q))
upper <- forecast$f + q_value * sqrt(unlist(forecast$Q))


df2<- data.frame(
  Time = time(Y),
  Observed = as.numeric(Y),
  Test = c(rep(NA, n-H), Y_test),
  Forecast = c(rep(NA, n-H),as.numeric(forecast$f)),
  Lower = c(rep(NA, n-H),lower),
  Upper = c(rep(NA, n-H),upper)
)

ggplot(df2, aes(x = Time)) +
  geom_line(aes(y = Observed, color = "Observed")) +
  geom_line(aes(y = Test, color = "Test"), na.rm = TRUE) +
  geom_line(aes(y = Forecast, color = "Forecast"), na.rm = TRUE) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), 
              fill = "red", alpha = 0.2, na.rm = TRUE) +
  scale_color_manual(values = c("Observed" = "black", 
                                "Test" = "green",
                                "Fitted" = "blue",
                                "Forecast" = "red")) +
  labs(x = "Year", y = "Y", color = "Series") +
  theme_minimal() +
  theme(legend.position = "bottom")

EQMP <- sum((as.numeric(Y_test)-forecast$f)^2)/H
EQMP

