# Archivo: R/funciones.R

library(igraph)
library(MASS)
library(ppcor)
library(qgraph)
library(dplyr)
library(tidyverse)

#' Divide y elimina la variable de respuesta
#'
#' @param df Data frame con los datos.
#' @param respuesta_var Nombre de la variable de respuesta.
#' @return Una lista con dos data frames: uno para cada valor de la variable de respuesta.
#' @export
divide_y_elimina <- function(df, respuesta_var) {
  data_1 <- df[df[[respuesta_var]] == 1, -which(names(df) == respuesta_var)]
  data_0 <- df[df[[respuesta_var]] == 0, -which(names(df) == respuesta_var)]
  return(list(data_1 = data_1, data_0 = data_0))
}

#' Calcula la matriz de correlación parcial
#'
#' @param df Data frame con los datos.
#' @param nivel_significacion Nivel de significación para los valores p (por defecto 0.05).
#' @return Matriz de correlación parcial ajustada.
#' @export
matriz_correlacion <- function(df, nivel_significacion = 0.05) {
  correlaciones <- pcor(df)
  mat_cor <- correlaciones$estimate
  p_values <- correlaciones$p.value
  mat_cor[p_values > nivel_significacion] <- 0
  for (i in 1:ncol(mat_cor)) mat_cor[i,i] <- 0
  return(mat_cor)
}

#' Realiza el test de permutaciones para la diferencia máxima de matrices de correlación
#'
#' @param df Data frame con los datos.
#' @param respuesta_var Nombre de la variable de respuesta.
#' @param n_permutaciones Número de permutaciones (por defecto 1000).
#' @param nivel_significacion Nivel de significación para los valores p (por defecto 0.05).
#' @return Vector de diferencias máximas de las matrices de correlación para cada permutación.
#' @export
test_permutaciones <- function(df, respuesta_var, n_permutaciones = 1000, nivel_significacion = 0.05) {
  resultados_maximos <- numeric(n_permutaciones)
  for (i in 1:n_permutaciones) {
    df_permutado <- df
    df_permutado[[respuesta_var]] <- sample(df[[respuesta_var]])
    division <- divide_y_elimina(df_permutado, respuesta_var)
    matriz_cor_1 <- matriz_correlacion(division$data_1, nivel_significacion)
    matriz_cor_0 <- matriz_correlacion(division$data_0, nivel_significacion)
    diferencia_maxima <- max(matriz_cor_1 - matriz_cor_0)
    resultados_maximos[i] <- diferencia_maxima
  }
  return(resultados_maximos)
}

#' Realiza el test de permutaciones para la diferencia total de matrices de correlación
#'
#' @param df Data frame con los datos.
#' @param respuesta_var Nombre de la variable de respuesta.
#' @param n_permutaciones Número de permutaciones (por defecto 1000).
#' @param nivel_significacion Nivel de significación para los valores p (por defecto 0.05).
#' @return Vector de diferencias totales de las matrices de correlación para cada permutación.
#' @export
test <- function(df, respuesta_var, n_permutaciones = 1000, nivel_significacion = 0.05) {
  resultados_IGS <- numeric(n_permutaciones)
  for (i in 1:n_permutaciones) {
    df_permutado <- df
    df_permutado[[respuesta_var]] <- sample(df[[respuesta_var]])
    division <- divide_y_elimina(df_permutado, respuesta_var)
    matriz_cor_1 <- matriz_correlacion(division$data_1, nivel_significacion)
    matriz_cor_0 <- matriz_correlacion(division$data_0, nivel_significacion)
    diferencia_cor <- sum(abs(matriz_cor_1 - matriz_cor_0))
    resultados_IGS[i] <- diferencia_cor
  }
  return(resultados_IGS)
}

#' Calcula diferencias estadísticas por nodo entre dos matrices de correlación
#'
#' @param matriz1 Primera matriz de correlación.
#' @param matriz2 Segunda matriz de correlación.
#' @param nodo Nombre del nodo para el cual se calcularán las diferencias.
#' @return Data frame con la diferencia media y el valor p del test t.
#' @export
calculo <- function(matriz1, matriz2, nodo) {
  if (!(nodo %in% colnames(matriz1))) 
    stop("El nodo no existe en la matriz1")
  if (!(nodo %in% colnames(matriz2))) 
    stop("El nodo no existe en la matriz2") 
  valores1 <- matriz1[, nodo]
  valores2 <- matriz2[, nodo]
  test <- t.test(valores1, valores2)
  p_valor <- test$p.value
  diferencia_media <- mean(valores1) - mean(valores2)
  results <- data.frame(Nodo = nodo, Diferencia_Media = round(diferencia_media, 3), P_valor = round(p_valor, 3))
  return(results)
}
# Archivo: R/funciones.R

library(igraph)
library(MASS)
library(ppcor)
library(qgraph)
library(dplyr)
library(tidyverse)

#' Divide y elimina la variable de respuesta
#'
#' @param df Data frame con los datos.
#' @param respuesta_var Nombre de la variable de respuesta.
#' @return Una lista con dos data frames: uno para cada valor de la variable de respuesta.
#' @export
divide_y_elimina <- function(df, respuesta_var) {
  data_1 <- df[df[[respuesta_var]] == 1, -which(names(df) == respuesta_var)]
  data_0 <- df[df[[respuesta_var]] == 0, -which(names(df) == respuesta_var)]
  return(list(data_1 = data_1, data_0 = data_0))
}

#' Calcula la matriz de correlación parcial
#'
#' @param df Data frame con los datos.
#' @param nivel_significacion Nivel de significación para los valores p (por defecto 0.05).
#' @return Matriz de correlación parcial ajustada.
#' @export
matriz_correlacion <- function(df, nivel_significacion = 0.05) {
  correlaciones <- pcor(df)
  mat_cor <- correlaciones$estimate
  p_values <- correlaciones$p.value
  mat_cor[p_values > nivel_significacion] <- 0
  for (i in 1:ncol(mat_cor)) mat_cor[i,i] <- 0
  return(mat_cor)
}

#' Realiza el test de permutaciones para la diferencia máxima de matrices de correlación
#'
#' @param df Data frame con los datos.
#' @param respuesta_var Nombre de la variable de respuesta.
#' @param n_permutaciones Número de permutaciones (por defecto 1000).
#' @param nivel_significacion Nivel de significación para los valores p (por defecto 0.05).
#' @return Vector de diferencias máximas de las matrices de correlación para cada permutación.
#' @export
test_permutaciones <- function(df, respuesta_var, n_permutaciones = 1000, nivel_significacion = 0.05) {
  resultados_maximos <- numeric(n_permutaciones)
  for (i in 1:n_permutaciones) {
    df_permutado <- df
    df_permutado[[respuesta_var]] <- sample(df[[respuesta_var]])
    division <- divide_y_elimina(df_permutado, respuesta_var)
    matriz_cor_1 <- matriz_correlacion(division$data_1, nivel_significacion)
    matriz_cor_0 <- matriz_correlacion(division$data_0, nivel_significacion)
    diferencia_maxima <- max(matriz_cor_1 - matriz_cor_0)
    resultados_maximos[i] <- diferencia_maxima
  }
  return(resultados_maximos)
}

#' Realiza el test de permutaciones para la diferencia total de matrices de correlación
#'
#' @param df Data frame con los datos.
#' @param respuesta_var Nombre de la variable de respuesta.
#' @param n_permutaciones Número de permutaciones (por defecto 1000).
#' @param nivel_significacion Nivel de significación para los valores p (por defecto 0.05).
#' @return Vector de diferencias totales de las matrices de correlación para cada permutación.
#' @export
test <- function(df, respuesta_var, n_permutaciones = 1000, nivel_significacion = 0.05) {
  resultados_IGS <- numeric(n_permutaciones)
  for (i in 1:n_permutaciones) {
    df_permutado <- df
    df_permutado[[respuesta_var]] <- sample(df[[respuesta_var]])
    division <- divide_y_elimina(df_permutado, respuesta_var)
    matriz_cor_1 <- matriz_correlacion(division$data_1, nivel_significacion)
    matriz_cor_0 <- matriz_correlacion(division$data_0, nivel_significacion)
    diferencia_cor <- sum(abs(matriz_cor_1 - matriz_cor_0))
    resultados_IGS[i] <- diferencia_cor
  }
  return(resultados_IGS)
}

#' Calcula diferencias estadísticas por nodo entre dos matrices de correlación
#'
#' @param matriz1 Primera matriz de correlación.
#' @param matriz2 Segunda matriz de correlación.
#' @param nodo Nombre del nodo para el cual se calcularán las diferencias.
#' @return Data frame con la diferencia media y el valor p del test t.
#' @export
calculo <- function(matriz1, matriz2, nodo) {
  if (!(nodo %in% colnames(matriz1))) 
    stop("El nodo no existe en la matriz1")
  if (!(nodo %in% colnames(matriz2))) 
    stop("El nodo no existe en la matriz2") 
  valores1 <- matriz1[, nodo]
  valores2 <- matriz2[, nodo]
  test <- t.test(valores1, valores2)
  p_valor <- test$p.value
  diferencia_media <- mean(valores1) - mean(valores2)
  results <- data.frame(Nodo = nodo, Diferencia_Media = round(diferencia_media, 3), P_valor = round(p_valor, 3))
  return(results)
}


