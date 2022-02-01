KB <- function(dmuLt,dmut,rho,C,nw,Phi,X) {
  
  nf<-dim(C)[1] #number of measurements
  nc<-dim(C)[2] #number of states
  dmuKt <- NA
  
  
  for (i in nf:nw){
    Y <- c(dmuLt[i], dmut[i-1])
    Q <- diag(nc) * rho
    R <- diag(nf) 
    
    #Kalman filter###############################################
    K <- (Phi+Q) %*% t(C) %*% solve(C %*% (Phi+Q) %*% t(C) + R)
    X <- X + K %*% (Y - C %*% X)
    Phi <- (diag(nc) - K %*% C) %*% (Phi+Q)
    #############################################################
    
    #Lo guardamos como condicion inicial para la proxima
    if (i==nf){
      Pp=Phi
      Xx=X
    }
    dmuKt <- c(dmuKt, dmuLt[i] + X[nc]) # almacena delta de kalman, que es lande + bias
    
  }
  
  output <- list(dmuKt[-1], Pp, Xx)
  return(output)
}