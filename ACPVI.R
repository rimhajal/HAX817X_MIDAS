# importation des packages
library(ggplot2)
library(tidyr)
library(dplyr)
library(magrittr) 
library(Metrics)
library(ade4)

# importer les données
données <- read.table("Regions.txt", header=TRUE, row.names=1)
X <- data.matrix(données[,1:34])
Y <- data.matrix(données[,35:65])
Y1 <- Y

#1er cas : X complet
suppression_aléatoire_3 <- function (tableau, proportion) {
  données <- as.vector(tableau)
  indicatrices <- matrix(ncol=5, nrow=nrow(tableau)*ncol(tableau))
  indicatrices[,1] <- données
  for (i in 1:length(données)) {
    if (is.na(données[i])==TRUE) {indicatrices[i,2]<-0}
    else {indicatrices[i,2]<-runif(1, min = 1, max = 1.999999)} # données uniformes entre 1 et 2 runif
  }
  indicatrices[,3] <- rep(1:nrow(tableau),ncol(tableau))
  indicatrices[,4] <- rep(1:ncol(tableau), each=nrow(tableau))
  indicatrices <- indicatrices[order(indicatrices[,2], decreasing=FALSE),]
  nombre_données <- nrow(tableau)*ncol(tableau) - sum(is.na(tableau))
  nombre_pseudomanquantes <- floor(proportion*nombre_données)
  indicatrices[(sum(is.na(tableau))+1):(sum(is.na(tableau))+nombre_pseudomanquantes+1),2] <- 0
  indicatrices[,5] <- indicatrices[,1]*floor(indicatrices[,2]) # partie entière de colonne 2
  indicatrices[,5]<-replace(indicatrices[,5], indicatrices[,5]==0, NA)
  for (k in 1:length(tableau)) {
    tableau[indicatrices[k,3],indicatrices[k,4]] <- indicatrices[k,5]
  }
  return(tableau)
}

imputer_colonne_par_moyenne <- function (données_dupliquées_bis) {
  moyenne_colonne <- as.vector(colMeans(données_dupliquées_bis, na.rm=TRUE))
  for (i in 1:nrow(données_dupliquées_bis)) {
    for (j in 1:ncol(données_dupliquées_bis)) {
      if (is.na(données_dupliquées_bis[i,j])==TRUE) {
        données_dupliquées_bis[i,j] <- moyenne_colonne[j]
      }
    }
  }
  return(données_dupliquées_bis)
}

# méthode de la SVD
méthode_svd <- function(données) {
  A <- as.matrix(données)
  TA <- t(A)
  ATA <- A %*% TA
  ATA.e <- eigen(ATA)
  u.mat <- ATA.e$vectors
  u.mat
  
  TAA <- TA %*% A
  TAA.e <- eigen(TAA)
  v.mat <- TAA.e$vectors[1:ncol(données),1:nrow(données)]
  v.mat
  
  r <- sqrt(ATA.e$values)
  r <- r * diag(length(r))[1:nrow(données),]
  r <- t(r)
  r
  
  svd.matrix <<- u.mat %*% r %*% t(v.mat)
  svd.matrix
}

Y <- suppression_aléatoire_3(Y, proportion=0.2)
données_manquantesY <- is.na(Y)
Y <- imputer_colonne_par_moyenne(Y)

inv <- solve(t(X)%*%X,tol = 1e-25)
A <- t(Y)%*%X%*%inv%*%t(X)%*%Y

eigenvalues <- eigen(A, symmetric = TRUE)$values
eigenvectors <- eigen(A, symmetric = TRUE)$vectors

D <- diag(eigenvalues)
P <- eigenvectors
v <- P %*% D %*% solve(P)

u <- matrix(nrow = ncol(X),ncol=ncol(Y))
for(h in 1:ncol(Y)){
  u[,h] <- (1/sqrt(abs(eigenvalues[h]))) * solve(t(X)%*%X,tol = 1e-30) %*% t(X) %*% Y %*% v[,h]
}
matT <- X%*%u%*%t(v)

c<-c()
for (h in 1:ncol(Y)){
  c <- append(c,sum(diag(t(Y)%*%matT[,h])))
}

Y_chapeau <- matrix(0,nrow = nrow(Y), ncol = ncol(Y))

for (h in 1:ncol(Y)){
  matrices <- c[h]*X%*%u[,h]%*%t(v[,h])
  Y_chapeau <- Y_chapeau + matrices
}

Y[données_manquantesY]<-Y_chapeau[données_manquantesY]
rmse(Y1,Y)

imputeddata <- tibble(imputed=Y[données_manquantesY],
                      truth = Y1[données_manquantesY])
visualisation <- (imputeddata %>%
                    ggplot(aes(truth, imputed))+
                    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "blue", linewidth = 1.25) +
                    geom_point()+
                    theme_bw(16)+
                    theme(plot.title = element_text(hjust = 0.5, size=12,face="bold"))+
                    theme (axis.title= element_text(size=rel(0.75)))+
                    labs(x="Valeurs Réelles", y="Valeurs Imputées"))

print(visualisation)

#2e cas : X incomplet
X <- data.matrix(données[,1:34])
X <- suppression_aléatoire_3(X, proportion=0.2)
données_manquantesX <- is.na(X)
X <- imputer_colonne_par_moyenne(X)
données_imputéesX <- méthode_svd(X)
X[données_manquantesX] <- données_imputéesX[données_manquantesX]

Y <- data.matrix(données[,35:65])
Y <- suppression_aléatoire_3(Y, proportion=0.2)
données_manquantesY <- is.na(Y)
Y <- imputer_colonne_par_moyenne(Y)

inv <- solve(t(X)%*%X,tol = 1e-25)
A <- t(Y)%*%X%*%inv%*%t(X)%*%Y

eigenvalues <- eigen(A, symmetric = TRUE)$values
eigenvectors <- eigen(A, symmetric = TRUE)$vectors

D <- diag(eigenvalues)
P <- eigenvectors
v <- P %*% D %*% solve(P)

u <- matrix(nrow = ncol(X),ncol=ncol(Y))
for(h in 1:ncol(Y)){
  u[,h] <- (1/sqrt(abs(eigenvalues[h]))) * solve(t(X)%*%X,tol = 1e-30) %*% t(X) %*% Y %*% v[,h]
}
matT <- X%*%u%*%t(v)

c<-c()
for (h in 1:ncol(Y)){
  c <- append(c,sum(diag(t(Y)%*%matT[,h])))
}

Y_chapeau <- matrix(0,nrow = nrow(Y), ncol = ncol(Y))

for (h in 1:ncol(Y)){
  matrices <- c[h]*X%*%u[,h]%*%t(v[,h])
  Y_chapeau <- Y_chapeau + matrices
}

Y[données_manquantesY]<-Y_chapeau[données_manquantesY]
rmse(Y1,Y)

imputeddata <- tibble(imputed=Y[données_manquantesY],
                      truth = Y1[données_manquantesY])
visualisation <- (imputeddata %>%
                    ggplot(aes(truth, imputed))+
                    geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "blue", linewidth = 1.25) +
                    geom_point()+
                    theme_bw(16)+
                    theme(plot.title = element_text(hjust = 0.5, size=12,face="bold"))+
                    theme (axis.title= element_text(size=rel(0.75)))+
                    labs(x="Valeurs Réelles", y="Valeurs Imputées"))

print(visualisation)