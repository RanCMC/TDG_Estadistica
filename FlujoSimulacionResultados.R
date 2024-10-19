############ LIBRERÍAS Y RECURSOS ############
library("data.table")
library("forecast")
library("strucchange")
library("changepoint")
library("fastcpd")
library("waveslim")
library("zoo")
library("boot")
source('ObtenerCaracteristicasSeries.R')
source('GenerarSimulacionesSeries.R')

##############################################
# IMPORTANTE:                                #
#     Llamar primero el SOURCE               #
#     'ObtenerCaracteristicasSeries.R'       #
#     'GenerarSimulacionesSeries.R'          #
##############################################

############ DATOS VARIABLES BASE ############
dfDemanda <- read.csv() # Dataframe.csv
dfTension <- read.csv() # Dataframe.csv

######################################################
# PRIMERA PARTE: OBTENER CARACTERÍTICAS DE LOS DATOS #
######################################################

# Esta información proviene de 'ObtenerCaracteristicasSeries.R'

# Pasos a seguir para obtener las características de las variables
# 1. Obtener los datos de las variables
listaDemanda <- dfDemanda$TotalDemanda
listaTension <- dfTension$TotalTensiones

# 2. Graficar las variables
graficarSerie(listaDemanda)
graficarSerie(listaTension)

# 3. Realizar un resumen estadístico de las variables (datos crudos)
resultadosResumen_estadistico_demanda <- resumen(listaDemanda)
resultadosResumen_estadistico_tension <- resumen(listaTension)

# 4. Modelo ARIMA
resultadosARIMA_demanda <- testArima(listaDemanda)
listaResiduos_demanda <- resultadosARIMA_demanda$listaResiduos_arima
resultadosARIMA_tension <- testArima(listaTension)
listaResiduos_tension <- resultadosARIMA_tension$listaResiduos_arima

# El modelo ARIMA de demanda tiene los p_d,d_d,q_d:
p_d <- resultadosARIMA_demanda$iAR
d_d <- resultadosARIMA_demanda$iI
q_d <- resultadosARIMA_demanda$iMA
listaCoeficientes_ARIMA_demanda <- resultadosARIMA_demanda$listaCoef_arima
ar_coef_d <- listaCoeficientes_ARIMA_demanda[grep("^ar", names(listaCoeficientes_ARIMA_demanda))]
ma_coef_d <- listaCoeficientes_ARIMA_demanda[grep("^ma", names(listaCoeficientes_ARIMA_demanda))]

# El modelo ARIMA de tensión tiene los p_t,d_t,q_t:
p_t <- resultadosARIMA_tension$iAR
d_t <- resultadosARIMA_tension$iI
q_t <- resultadosARIMA_tension$iMA
listaCoeficientes_ARIMA_tension <- resultadosARIMA_tension$listaCoef_arima
ar_coef_t <- listaCoeficientes_ARIMA_tension[grep("^ar", names(listaCoeficientes_ARIMA_tension))]
ma_coef_t <- listaCoeficientes_ARIMA_tension[grep("^ma", names(listaCoeficientes_ARIMA_tension))]

# 5. Realizar un resumen estadístico de las variables (Modelo ARIMA)

# Los residuales de demanda siguen una distribución ~ N(mu_d,sigma2_d)
resultadosResumen_estadistico_residuos_demanda <- resumen(listaResiduos_demanda)
mu_d <- resultadosResumen_estadistico_residuos_demanda$media # Lista proveniente de SOURCE
sigma2_d <- resultadosResumen_estadistico_residuos_demanda$std # Lista proveniente de SOURCE

# Los residuales de tensión siguen una distribución ~ N(mu_t,sigma2_t)
resultadosResumen_estadistico_residuos_tension <- resumen(listaResiduos_tension)
mu_t <- resultadosResumen_estadistico_residuos_tension$media # Lista proveniente de SOURCE
sigma2_t <- resultadosResumen_estadistico_residuos_tension$std # Lista proveniente de SOURCE

# 6.Obtener el mejor grado de polinómio que se ajuste a los datos
# Nota: Para los datos anteriores se halló ya los valores
# Se deja la metodología propuesta

resultadosPoly_demanda <- MejorModeloTendencia(listaDemanda)
listaResiduales_poly_demanda <- EliminarTendencia(resultadosPoly_demanda$iGrado, 
                                                  resultadosPoly_demanda$dfPoly)
resultadosPoly_tension <- MejorModeloTendencia(listaTension)
listaResiduales_poly_tension <- EliminarTendencia(resultadosPoly_tension$iGrado, 
                                                  resultadosPoly_tension$dfPoly)

# Los residuales de polydemanda siguen una distribución ~ N(mu_d,sigma2_d)
summary(listaResiduales_poly_demanda)
mu_poly_d <- mean(listaResiduales_poly_demanda)
sigma2_poly_d <- sd(listaResiduales_poly_demanda)
listaPrediccion_modelo_lineal_demanda <- regresionLineal(listaResiduales_poly_demanda)

# Los residuales de polytensión siguen una distribución ~ N(mu_t,sigma2_t)
summary(listaResiduales_poly_tension)
mu_poly_t <- mean(listaResiduales_poly_tension)
sigma2_poly_t <- sd(listaResiduales_poly_tension)
listaPrediccion_modelo_lineal_tension <- regresionLineal(listaResiduales_poly_tension)

######################################################
# SEGUNDA PARTE: FUNCIONES SIMULADORAS DE DATOS      #
######################################################

# Esta información proviene de 'GenerarSimulacionesSeries.R'

############ SERIE ~N(mu,sigma^2) ############
requerimientosNormal <- function(listaVariable,mu,sigma_2){
  datosRequerimientos_normal <- list(listaVariable=listaVariable,
                                     mu=mu,
                                     sigma_2=sigma_2,
                                     sTipo_datos='Norm')
  return(datosRequerimientos_normal)
}

############ SERIE REGRESION LINEAL ############
requerimientosLm <- function(listaResiduos_modelo_poly,listaPrediccion_modelo_lineal){
  datosRequerimientos_lm <- list(listaResiduos_modelo_poly=listaResiduos_modelo_poly,
                                 listaPrediccion_modelo_lineal=listaPrediccion_modelo_lineal,
                                 sTipo_datos='lm')
return(datosRequerimientos_lm)
}

############ SERIE ARIMA PDQ ############
requerimientosARIMA <- function(listaVariable,iAR,iI,iMA,ar_coef,ma_coef){
  datosRequerimientos_ARIMA <- list(listaVariable=listaVariable,
                                    iAR=iAR,
                                    iI=iI,
                                    iMA=iMA,
                                    ar_coef=ar_coef,
                                    ma_coef=ma_coef,
                                    sTipo_datos='ARIMA')
  return(datosRequerimientos_ARIMA)
}
############ FUNCIÓN SIMULAR ############

generarSimulaciones <- function(datosRequerimientos,iDuracion,iMultiplicador){
  
  sTipo_datos <- datosRequerimientos$sTipo_datos
  if (sTipo_datos == 'Norm'){
    listaVariable <- datosRequerimientos$listaVariable
    mu <- datosRequerimientos$mu
    sigma_2 <- datosRequerimientos$sigma_2
    listaSimulada <- simularDatos_normal(listaVariable,mu,sigma_2)
    iAmplitud <- sigma_2*iMultiplicador
  
    } else if (sTipo_datos == 'lm'){
      listaResiduos_modelo_poly <- datosRequerimientos$listaResiduos_modelo_poly
      listaPrediccion_modelo_lineal <- datosRequerimientos$listaPrediccion_modelo_lineal
      mu <- mean(listaResiduos_modelo_poly)
      sigma_2 <- sd(listaResiduos_modelo_poly)
      listaSimulada <- simularDatos_lineal(listaResiduos_modelo_poly,listaPrediccion_modelo_lineal)
      iAmplitud <- sigma_2*iMultiplicador
      
    } else if (sTipo_datos == 'ARIMA'){
      listaVariable <- datosRequerimientos$listaVariable
      iAR <- datosRequerimientos$iAR
      iI <- datosRequerimientos$iI
      iMA <- datosRequerimientos$iMA
      ar_coef <- datosRequerimientos$ar_coef
      ma_coef <- datosRequerimientos$ma_coef
      listaSimulada <- simularDatos_ARIMA(listaVariable,iAR,iI,iMA,ar_coef,ma_coef)
      sigma_2 <- sd(diff(listaSimulada)) 
      mu <- mean(listaSimulada)
      iAmplitud <- (abs(max(listaVariable))-abs(min(listaVariable)))*iMultiplicador/50 # si 5 -> 10%, si 10 -> 20%

  }
  iCambio <- length(listaSimulada)*2/3
  
  listaSimulada_con_eventos <- serieCon_evento(listaSimulada,iMultiplicador,iCambio,iDuracion,mu,sigma_2)
  
  resultadosMuestreos <- muestreoSerie(listaSimulada_con_eventos)
  return(resultadosMuestreos)
}


######################################################
# TERCERA PARTE: FUNCIONES DE METODOLOGÍAS           #
######################################################

############ CHANGEPOINT ############

funChangepoint <- function(listaSimulada){
  start_time <- Sys.time()
  
  CP <- cpt.meanvar(listaSimulada)
  listaHallazgos_CP <- CP@cpts
  
  end_time <- Sys.time()
  time_diff <- end_time - start_time
  
  resultadoCP <- list(listaHallazgos_CP = listaHallazgos_CP,
                      time_diff = time_diff)
  return(resultadoCP)
}
############ STRUCCHANGE ############

funStrucchange <- function(listaSimulada){
  
  n <- length(listaSimulada)
  t <-  1:n
  
  start_time <- Sys.time()
  
  SC <- breakpoints(listaSimulada~t)
  listaHallazgos_SC <- SC$breakpoints
  
  end_time <- Sys.time()
  time_diff <- end_time - start_time
  
  resultadoSC <- list(listaHallazgos_SC = listaHallazgos_SC,
                      time_diff = time_diff)
  return(resultadoSC)
}
############ FASTCPD ############

funFastcpd_ARIMA <- function(listaSimulada){
  
  start_time <- Sys.time()
  
  modArima <- auto.arima(listaSimulada)
  iAR <- modArima$arma[1]
  iMA <- modArima$arma[2]
  iI <- ndiffs(listaSimulada)
  if (iAR == 0 && iMA == iMA && iI == 0){
    # iAR <- 1
    # iI <- 1
    iMA <- 1
  }
  
  fcpARIMA <- fastcpd.arima(
    data = listaSimulada,
    order = c(iAR, iI, iMA),
    segment_count = 3,
    lower = c(-1, -1, -1, -1, -1, 1e-10),
    upper = c(1, 1, 1, 1, 1, Inf),
    line_search = c(1, 0.1, 1e-2)
  )

  end_time <- Sys.time()
  time_diff <- end_time - start_time
  
  listaHallazgos_Fcpd <- fcpARIMA@cp_set
  
  resultadoFcpd <- list(listaHallazgos_Fcpd = listaHallazgos_Fcpd,
                      time_diff = time_diff)
  return(resultadoFcpd)
}


######################################################
# CUARTA PARTE: FUNCIÓN PRINCIPAL                    #
######################################################
# Con esta función se espera realizar un número X de
# iteraciones para realizar las simulaciones de las 
# series de tiempo y verificar el comportamiento y la
# potencia de las metodologías.


##########################################################
##########################################################
##########################################################
#
# ------------- EJECUCIÓN DEL CÓDIGO COMPLETO ---------- #
#
##########################################################
##########################################################
##########################################################



######################################################
# QUINTA PARTE: FUNCIÓN SIMULACIÓN PARALELIZACIÓN    #
######################################################
# En esta parte se pretende realizar la generación de
# gráficas para presentar los resultados obtenidos con
# las simulaciones

library(foreach)
library(doParallel)


funcionPrincipal_iteraciones <- function(argumentos) {  
  iIteraciones <- argumentos$iIteraciones
  iDuracion <- argumentos$iDuracion
  iMultiplicador <- argumentos$iMultiplicador
  datosRequerimientos <- argumentos$datosRequerimientos
  listaVariable <- argumentos$datosRequerimientos$listaVariable
  
  set.seed(5686798)
  sTipo_datos <- datosRequerimientos$sTipo_datos
  listaColumnas_resultados_iteraciones <- c('Iteración', 'Tipo Modelo', 'Muestreo [min]', 
                                            'Amplitud', 'Duración CP [s]', 'Duración SC [s]', 
                                            'Duración FCPD [s]', 'Detección CP', 'Detección SC', 
                                            'Detección FCPD')
  
  funcionDeteccion_metodologia <- function(listaHallazgos, listaRango) {
    interseccion <- intersect(unlist(listaHallazgos), unlist(listaRango))
    return(ifelse(length(interseccion) > 0, 1, 0))
  }
  
  # Configurar el cluster de procesamiento paralelo
  num_cores <- detectCores() - 1  # Usar un núcleo menos que el total disponible
  cl <- makeCluster(num_cores)
  # cl <- makeCluster(1)
  registerDoParallel(cl)
  
  # Cargar librerías en el entorno paralelo
  clusterEvalQ(cl, {
    library("data.table")
    library("forecast")
    library("strucchange")
    library("changepoint")
    library("fastcpd")
    library("waveslim")
    library("zoo")
    library("boot")
    source('ObtenerCaracteristicasSeries.R')
    source('GenerarSimulacionesSeries.R')
  })
  
  # Exportar funciones necesarias al entorno paralelo
  clusterExport(cl, c("generarSimulaciones", "funChangepoint", "funStrucchange", 
                      "funFastcpd_ARIMA","simularDatos_normal","regresionLineal",
                      "simularDatos_lineal","simularDatos_ARIMA",
                      "serieCon_evento","generarEvento","muestreoSerie"))
  
  
  # Procesar iteraciones en paralelo usando foreach
  
  resultado_paralelo <- foreach(i = 1:iIteraciones, .packages = c("data.table",
                                                                  "forecast",
                                                                  "strucchange",
                                                                  "changepoint",
                                                                  "fastcpd",
                                                                  "waveslim",
                                                                  "zoo",
                                                                  "boot")) %dopar% {
                                                                    resultadosMuestreos_simulaciones <- generarSimulaciones(datosRequerimientos, iDuracion, iMultiplicador)
                                                                    listaMinutos_1 <- resultadosMuestreos_simulaciones$listaSimulada_1
                                                                    listaMinutos_5 <- resultadosMuestreos_simulaciones$listaSimulada_5
                                                                    listaMinutos_10 <- resultadosMuestreos_simulaciones$listaSimulada_10
                                                                    obtener_resultados <- function(iMuestreo, listaMinutos, iRango) {
                                                                      iLongitud_lista <- length(listaMinutos)
                                                                      iPosicion_evento <- iLongitud_lista * 2 / 3
                                                                      listaRango_deteccion <- seq(iPosicion_evento - iRango, iPosicion_evento + iRango)
                                                                      resultadoCP <- funChangepoint(listaMinutos)
                                                                      boolHallazgo_cp <- funcionDeteccion_metodologia(resultadoCP$listaHallazgos_CP, listaRango_deteccion)
                                                                      resultadoSC <- funStrucchange(listaMinutos)
                                                                      boolHallazgo_sc <- funcionDeteccion_metodologia(resultadoSC$listaHallazgos_SC, listaRango_deteccion)
                                                                      resultadoFcpd <- funFastcpd_ARIMA(listaMinutos)
                                                                      boolHallazgo_fcpd <- funcionDeteccion_metodologia(resultadoFcpd$listaHallazgos_Fcpd, listaRango_deteccion)
                                                                      
                                                                      c(i, sTipo_datos, iMuestreo, iMultiplicador, resultadoCP$time_diff, resultadoSC$time_diff,
                                                                        resultadoFcpd$time_diff, boolHallazgo_cp, boolHallazgo_sc, boolHallazgo_fcpd)
                                                                    }
                                                                    
                                                                    resultados_1 <- obtener_resultados(1, listaMinutos_1, 10)
                                                                    resultados_5 <- obtener_resultados(5, listaMinutos_5, 7)
                                                                    resultados_10 <- obtener_resultados(10, listaMinutos_10, 3)
                                                                    
                                                                    return(list(dfResultados_1 = resultados_1, dfResultados_5 = resultados_5, dfResultados_10 = resultados_10))
                                                                  }
  
  # Detener el cluster
  stopCluster(cl)
  
  # Convertir resultados en data frames
  dfResultados_1 <- do.call(rbind, lapply(resultado_paralelo, `[[`, "dfResultados_1"))
  dfResultados_5 <- do.call(rbind, lapply(resultado_paralelo, `[[`, "dfResultados_5"))
  dfResultados_10 <- do.call(rbind, lapply(resultado_paralelo, `[[`, "dfResultados_10"))
  
  colnames(dfResultados_1) <- listaColumnas_resultados_iteraciones
  colnames(dfResultados_5) <- listaColumnas_resultados_iteraciones
  colnames(dfResultados_10) <- listaColumnas_resultados_iteraciones
  
  resultadosMetodologias_iteraciones <- list(dfResultados_1 = as.data.frame(dfResultados_1),
                                             dfResultados_5 = as.data.frame(dfResultados_5),
                                             dfResultados_10 = as.data.frame(dfResultados_10))
  
  # Si quieres graficar en la primera iteración, puedes usar código similar al anterior
  # Debes agregar un paso para manejar gráficos o visualizaciones según tu necesidad
  
  return(resultadosMetodologias_iteraciones)
}

