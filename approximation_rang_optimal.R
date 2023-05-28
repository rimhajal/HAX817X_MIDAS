# importation des packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(magrittr) 
library(Metrics)

# importer les données
données <- read.table("Regions.txt", header=TRUE, row.names=1)
données <- data.matrix(données[,1:34])

méthodeSvd <- function(données, M = 1) {
  résultat <<- svd(données)
  with(résultat,
       u[, 1:M, drop = FALSE] %*%
         (d[1:M] * t(v[, 1:M, drop = FALSE])))
}

approxRangOpt <- function(données, standardiser, K, proportionV, proportionF,méthode){
  if (standardiser==1) {standardiser(données)}
  else {données <- données}
  dupliquer(données)
  moyenne_colonne <<- as.vector(colMeans(données_dupliquées_bis, na.rm=TRUE))
  if (méthode==1) {suppression_aléatoire_1(données_dupliquées,données_dupliquées_bis, proportionV)}
  if (méthode==2) {suppression_aléatoire_2(données_dupliquées, données_dupliquées_bis, proportionV)}
  if (méthode==3) {suppression_aléatoire_3(données_dupliquées_bis, proportionV)}
  EQMP <<- c()
  for (j in 1:15){
      for (k in 1:K) {
       if (méthode==1) {suppression_aléatoire_1(données_dupliquées, données_dupliquées_bis, proportionF)}
       if (méthode==2) {suppression_aléatoire_2(données_dupliquées, données_dupliquées_bis, proportionF)}
       if (méthode==3) {suppression_aléatoire_3(données_dupliquées_bis, proportionF)}
       données_manquantes <<- is.na(données_dupliquées)
       imputer_colonne_par_moyenne(données_dupliquées_bis)
       données_imputées <<- méthodeSvd(données_dupliquées_bis, M=1)
       données_dupliquées[données_manquantes] <<- données_imputées[données_manquantes]
       EQMP <<- append(EQMP,rmse(données_dupliquées[,j],données_imputées[,j]))
      }
    cat("Rang : ", j, "Erreur : ", colMeans(matrix(EQMP,K))[j] , "\n")
  }
  visualiser_donnéees_imputées_svd(données,données_imputées)
}
