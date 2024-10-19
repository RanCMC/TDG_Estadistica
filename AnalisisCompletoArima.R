# Cargar las librerías necesarias
library(forecast)
library(tseries)

set.seed(123)
serieTiempo <- ts(rnorm(100), frequency = 12)  

# 1. Analizar la serie de tiempo: visualización y prueba de estacionariedad
plot(serieTiempo, main="Serie de Tiempo", ylab="Valores", xlab="Tiempo")

# Graficar la ACF y PACF de la serie de tiempo original
par(mfrow = c(1, 2)) 
acf(serieTiempo, main="ACF de la Serie de Tiempo")
pacf(serieTiempo, main="PACF de la Serie de Tiempo")

# Realizar la prueba de Dickey-Fuller aumentada (ADF) para verificar la estacionariedad
adf_test <- adf.test(serieTiempo, alternative = "stationary")
print(adf_test)

# 2. Identificar el mejor modelo ARIMA usando auto.arima
mejor_arima <- auto.arima(serieTiempo)
summary(mejor_arima)

# 3. Graficar el ajuste del modelo
par(mfrow = c(1, 1))  # Volver a un gráfico por ventana
plot(forecast(mejor_arima, h=12), main="Forecast con ARIMA")

# 4. Análisis de los residuales
# Graficar los residuales del modelo
residuals_arima <- residuals(mejor_arima)
par(mfrow = c(2, 2))
plot(residuals_arima, main="Residuales del Modelo ARIMA", ylab="Residuales")
acf(residuals_arima, main="ACF de los Residuales")
pacf(residuals_arima, main="PACF de los Residuales")
qqnorm(residuals_arima); qqline(residuals_arima, col="red")  # QQ-plot

# Prueba de Ljung-Box para evaluar la autocorrelación en los residuales
box_test <- Box.test(residuals_arima, type="Ljung-Box")
print(box_test)

# Verificar la normalidad de los residuales
shapiro_test <- shapiro.test(residuals_arima)
print(shapiro_test)

# 5. Interpretación
# Si el ADF test no rechaza la hipótesis nula, es posible que necesites diferenciar la serie.
# Si los residuales no son ruido blanco (Ljung-Box p-valor bajo), considera ajustar el modelo.
