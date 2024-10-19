############ LIBRERÍAS Y RECURSOS ############
library("forecast")
library("zoo")
library("boot")

dfDemanda <- read.csv() # Dataframe.csv
dfTension <- read.csv() # Dataframe.csv

############ GRAFICAR SERIE ############

graficarSerie <- function(listaVariable){
  plot(listaVariable,type='l',ylab='Variable',main='Comportamiento de varaible sin evento')
  grid()
}

############ RESUMEN ESTADISTICO ############

resumen <- function(lista){
  summary(lista)
  media <- mean(lista)
  std <- sd(lista)
  resultadosResumen_estadistico <- list(media=media,std = std)
  return(resultadosResumen_estadistico)
}

############ MODELO ARIMA ############
testArima <- function(tsSerie_tiempo) {
  modArima <- auto.arima(tsSerie_tiempo)
  iAR <- modArima$arma[1]
  iMA <- modArima$arma[2]
  iI <- ndiffs(tsSerie_tiempo)
  listaResiduos_arima <- modArima$residuals
  listaCoef_arima <- modArima$coef
  resultadosArima <- list(iAR=iAR,iI=iI,iMA=iMA,
                     listaResiduos_arima=listaResiduos_arima,
                     listaCoef_arima=listaCoef_arima)
  return(resultadosArima)
}

############ GRADOS DE POLINOMIOS ############
MejorModeloTendencia <- function(tsSerie_tiempo) {
  dfPoly <- data.frame(tsVariable = tsSerie_tiempo,
                       tiempo=c(1:length(tsSerie_tiempo)))
  listaCV <- vector()
  for (iGrado in 1:13){
    modeloPoly <- glm(tsVariable ~ poly(tiempo,iGrado), data = dfPoly)
    cv.err.k1 <- cv.glm(dfPoly,modeloPoly, K=10)
    listaCV[iGrado] <- cv.err.k1$delta[1]
  }
  plot(listaCV,type='l',main="Validación cruzada",xlab="Grado",
       
       ylab="Error")
  grid()
  iGrado <- which.min(diff(listaCV))+1
  resultadosGrados = list(iGrado=iGrado,listaCV=listaCV, dfPoly=dfPoly)
  
  return(resultadosGrados)
}

EliminarTendencia <- function(iGrado, dfPoly) {
  
  modeloPoly = lm(tsVariable ~ poly(tiempo,iGrado),data=dfPoly)
  listaResiduos_modelo_poly <- modeloPoly$residual
  return(listaResiduos_modelo_poly)
}
