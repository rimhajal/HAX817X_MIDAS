# importation des packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(magrittr) 
library(Metrics)

# importer les données
données <- read.table("/Users/leeloo/Documents/SSD/semestre\ 2/projet/Regions.txt", header=TRUE, row.names=1)
données <- data.matrix(données)
X <- data.matrix(données[,1:34])
Y <- data.matrix(données[,35:65])

# standardiser les données
standardiser <- function (X, Y) {
  nombrelignes <- nrow(données)
  X <<- scale(X, center=TRUE) * sqrt((nombrelignes-1)/nombrelignes)
  Y <<- scale(Y, center=TRUE) * sqrt((nombrelignes-1)/nombrelignes)
  return(X, Y)
}

# supprimer aléatoirement les données
suppression_aléatoire_3 <- function (Y, proportion) {
  données <- as.vector(Y)
  indicatrices <- matrix(ncol=5, nrow=nrow(Y)*ncol(Y))
  indicatrices[,1] <- données
  for (i in 1:length(données)) {
    if (is.na(données[i])==TRUE) {indicatrices[i,2]<-0}
    else {indicatrices[i,2]<-runif(1, min = 1, max = 1.999999)} # données uniformes entre 1 et 2 runif
  }
  indicatrices[,3] <- rep(1:nrow(Y),ncol(Y))
  indicatrices[,4] <- rep(1:ncol(Y), each=nrow(Y))
  indicatrices <- indicatrices[order(indicatrices[,2], decreasing=FALSE),]
  nombre_données <- nrow(Y)*ncol(Y) - sum(is.na(Y))
  nombre_pseudomanquantes <- floor(proportion*nombre_données)
  indicatrices[(sum(is.na(Y))+1):(sum(is.na(Y))+nombre_pseudomanquantes+1),2] <- 0
  indicatrices[,5] <- indicatrices[,1]*floor(indicatrices[,2]) # partie entière de colonne 2
  indicatrices[,5]<-replace(indicatrices[,5], indicatrices[,5]==0, NA)
  for (k in 1:length(Y)) {
    Y[indicatrices[k,3],indicatrices[k,4]] <<- indicatrices[k,5]
  }
}

# visualiser les données manquantes
visualisation_données_manquantes <- function (Y) {
  visualisation <- (Y %>%
                      as.data.frame() %>%
                      mutate(row_id=row_number()) %>%
                      pivot_longer(-row_id, names_to="feature", values_to="valeurs") %>%
                      ggplot(aes(x=feature,y=row_id,fill=valeurs))+
                      geom_tile()+
                      theme(axis.text.x = element_text(angle = 90, size=5))+
                      labs(x = "variable", y = "numéro de ligne")+
                      ggtitle("Représentation Graphique des Données Manquantes")+
                      scale_fill_continuous(na.value = 'red'))
  ggsave("Visualisation des Données Manquantes.png", width=5, height=10)
  print(visualisation)
}

# remplacer les données manquantes par la moyenne de la colonne
imputer_colonne_par_moyenne <- function (données_dupliquées_bis) {
  moyenne_colonne <<- as.vector(colMeans(données_dupliquées_bis, na.rm=TRUE))
  for (i in 1:nrow(données_dupliquées_bis)) {
    for (j in 1:ncol(données)) {
      if (is.na(données_dupliquées_bis[i,j])==TRUE) {
        données_dupliquées_bis[i,j] <<- moyenne_colonne[j]
      }
    }
  }
}

# méthode de la SVD
méthode_svd <- function(données) {
  A <- as.matrix(données)
  TA <- t(A)
  ATA <- A %*% TA
  ATA.e <- eigen(ATA)
  u <- ATA.e$vectors
  u
  
  TAA <- TA %*% A
  TAA.e <- eigen(TAA)
  v <- TAA.e$vectors[1:ncol(données),1:nrow(données)]
  v
  
  r <- sqrt(ATA.e$values)
  r <- r * diag(length(r))[1:nrow(données),]
  r <- t(r)
  r
  
  svd_matrice <<- u %*% r %*% t(v)
  svd_matrice
}

# remplacer les donnés manquantes par la méthode de la SVD
imputation_svd <- function (données_dupliquées) {
  epsilon <<- 1e-7
  erreur <<- 1
  iteration <<- 0
  données_manquantes <<- is.na(données_dupliquées)
  mssold <<- mean((scale(données_dupliquées, moyenne_colonne, FALSE)[!données_manquantes])^2)
  mss0 <<- mean(données_dupliquées[!données_manquantes]^2)
  while(erreur > epsilon){           
    iteration <<- iteration + 1
    données_imputées <<- méthode_svd(données_dupliquées_bis)
    données_dupliquées[données_manquantes] <<- données_imputées[données_manquantes]
    mss <<- mean(((données_dupliquées-données_imputées)[!données_manquantes])^2)
    erreur <<- (mssold - mss)/mss0
    mssold <<- mss
    eqmp <- rmse(données, données_imputées)
    cat("Rang: ", iteration, "Erreur Quadratique Moyenne de Prédiction:", eqmp)
  }
}

# visualisation des données imputées
visualiser_donnéees_imputées_svd <- function (données_imputées, données) {
  imputeddata <<- tibble(imputed=données_imputées[données_manquantes],
                         truth = données[données_manquantes])
  visualisation <- (imputeddata %>%
                      ggplot(aes(truth, imputed))+
                      geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "blue", linewidth = 1.25) +
                      geom_point()+
                      theme_bw(16)+
                      labs(color = "y = x")+
                      ggtitle("Comparaison des Valeurs Réelles et Imputées")+
                      theme(plot.title = element_text(hjust = 0.5, size=12,face="bold"))+
                      theme (axis.title= element_text(size=rel(0.75)))+
                      labs(x="Valeurs Réelles", y="Valeurs Imputées"))
  ggsave("Régression Linéaire des Données Imputées.png", width=5, height=10)
  print(visualisation)
}

# compilation globale
imputation <- function (données, X, Y, standardiser, proportion = 0.2) {
  if (standardiser==1) {standardiser(X, Y)}
  suppression_aléatoire_3(Y, 0.2)
  visualisation_données_manquantes(Y)
  données_dupliquées <<- cbind(X,Y)
  données_dupliquées_bis <<- cbind(X,Y)
  imputer_colonne_par_moyenne(données_dupliquées_bis)
  imputation_svd(données_dupliquées)
  visualiser_donnéees_imputées_svd(données_imputées, données)
}
imputation(données, X, Y, 0, 0.2)
rmse(données[données_manquantes], données_imputées[données_manquantes])
