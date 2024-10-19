############ LIBRERÍAS Y RECURSOS ############
library("strucchange")
library("changepoint")

############ GENERAR SIMULACIONES ############

############ MUESTREO DE LOS DATOS ############
# La siguiente función es general para todos los 
# datos simulados. Esta función lo que hará es obtener
# datos cada minuto, 5 minutos y 10 minutos de la
# serie original

muestreoSerie <- function(listaSimulada){
  
  listaSimulada_1 <- listaSimulada[seq(1,length(listaSimulada),by=60)]
  listaSimulada_5 <- listaSimulada[seq(1,length(listaSimulada),by=300)]
  listaSimulada_10 <- listaSimulada[seq(1,length(listaSimulada),by=600)]
  
  resultadosMuestreos <- list(listaSimulada_1 = listaSimulada_1,
                         listaSimulada_5 = listaSimulada_5,
                         listaSimulada_10 = listaSimulada_10)
  
  return(resultadosMuestreos)
}

############ GENERAR EVENTO ############
# Con estas funciones lo que se pretende es simular
# el comportamiento de las variables eléctricas
# ante un evento en el sistema.

generarEvento <- function(listaSimulada,iCambio,iAmplitud,iDuracion, mu, sigma2){
  
    # Definir los coeficientes del polinomio
    b <- 0  # coeficiente de x
    c <- listaSimulada[iCambio]-iAmplitud  # término constante
    a <- (listaSimulada[iCambio]-c)/(iDuracion)^2  # coeficiente de x^2
    
    # Definir la función del polinomio
    polinomio_cuadratico <- function(x) {
      a* (x-iCambio)^2 + b * x + c
    }
    
    # Generar una secuencia de puntos x
    x <- seq(iCambio, iCambio+iDuracion, length.out = iDuracion)
    
    # Calcular los valores del polinomio en esos puntos
    y_polinomio <- polinomio_cuadratico(x)
    
    # Definir los parámetros del ruido blanco
    n <- length(x)  # Número de puntos de ruido blanco a generar
    mean <- 0       # Media del ruido blanco
    sd <- sigma2         # Desviación estándar del ruido blanco
    
    # Generar ruido blanco
    ruido_blanco <- rnorm(n, mean, sd)
    listaEvento_simulado_con_ruido <- y_polinomio + ruido_blanco

  
  return(listaEvento_simulado_con_ruido)
}

serieCon_evento <- function(listaSimulada,iMultiplicador,iCambio,iDuracion,mu,sigma_2){

  n <- length(listaSimulada)
  iAmplitud <- iMultiplicador*sigma_2
  y <- generarEvento(listaSimulada,iCambio,iAmplitud,
                    iDuracion, mu, sigma_2)
  
  m <- c(rep(0, iCambio),y,rep(0, (n)-(iCambio+length(y))))
  
  listaSimulada_ <- c(listaSimulada[0:iCambio],rep(0, iDuracion),
                    listaSimulada[(iCambio+iDuracion):(n-1)])
  ax <- listaSimulada_+ m
  return(ax)
}

############ SERIE ~N(mu,sigma^2) ############
# Aquí se tendrá en cuenta el valor de mu y sigma^2
# de los residuales del modelo ARIMA de los datos

# Los residuales de demanda siguen una distribución ~ N(mu_d,sigma2_d)
# Los residuales de tensión siguen una distribución ~ N(mu_t,sigma2_t)

# La longitud de los datos completos es de 60*60*24
# lo que implica que es un dato cada segundo, osea
# 86400 valores en una serie de tiempo. Esta misma cantidad
# se debe simular

simularDatos_normal<- function(listaVariable,mu,sigma2){
  
  iLongitud_datos <- length(listaVariable)
  
  listaNormal_simulada_completa <- rnorm(iLongitud_datos,
                                         mean = mu,
                                         sd = sigma2)
  
  return(listaNormal_simulada_completa)
}

############ SERIE REGRESION LINEAL ############
# En esta función se tendrá en cuenta los residuales
# del polinomio del conjunto de datos.

# Se debe realizar una regresión lineal y a esta
# añadirle ruido blanco con el sd de los residuales
# del polinomio de la varaible

regresionLineal <- function(listaResiduos_modelo_poly){
  n <- length(listaResiduos_modelo_poly)
  t <- 1:n
  
  modeloLineal <- lm(listaResiduos_modelo_poly~t)
  
  listaPrediccion_modelo_lineal <- predict(modeloLineal)
  return(listaPrediccion_modelo_lineal)
}

simularDatos_lineal <- function(listaResiduos_modelo_poly, listaPrediccion_modelo_lineal){
  n <- length(listaResiduos_modelo_poly)
  t <- 1:n
  
  sigma_2_residuales_poly <- sd(listaResiduos_modelo_poly)
  
  listaRuido_blanco_variable <- rnorm(n,0,sigma_2_residuales_poly)
  listaModelo_con_ruido <- listaPrediccion_modelo_lineal + listaRuido_blanco_variable
  return(listaModelo_con_ruido)
}

############ SERIE ARIMA PDQ ############
simularDatos_ARIMA <- function(listaVariable,iAR,iI,iMA,ar_coef,ma_coef){
  
  n <- length(listaVariable)
  listaARIMA_simulada <- arima.sim(n = n, 
               list(order = c(iAR,iI,iMA), 
                    ar = ar_coef, 
                    ma = ma_coef))
  return(listaARIMA_simulada)
}
