# This is the algorithm of the Kalman filter. The equations are given in Appendixes A and B.

KB <- function(dmuLt,dmut,rho,C,L.dyn,Phi,X) {
  
  dmuKt <- NA # Declare an array to store the predictions inside the window
  
  for (j in n.rows:L.dyn){
    Y <- c(dmuLt[j], dmut[j-1]) # Measurements
    I <- diag(2)  # Diagonal matrix of size 2x2       
    
    #Kalman filter###############################################
    K <- (Phi+I*rho) %*% t(C) %*% solve(C %*% (Phi+I*rho) %*% t(C) + I)
    X <- X + K %*% (Y - C %*% X)
    Phi <- (I - K %*% C) %*% (Phi+I*rho)
    #############################################################
    
    if (j == n.rows){ Pp=Phi ;  Xx=X  } # This is stored for the start of the next window
    dmuKt <- c(dmuKt, dmuLt[j] + X[2]) # Append the prediction of the new method for j,
                                       # which is given by the breeder's prediction plus the bias
    
  }
  
  output <- list(dmuKt[-1], Pp, Xx)
  return(output)
}