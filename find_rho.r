test.rho<-function(){
  for (t in 1:n.traits) {
    
    # dmuLt stores the prediction of the breeder's equation for trait t
    dmuLt <- Fs[t, ]
    
    # dmut stores the measured changes in trait means for trait t.
    # there is an extra zero because we have not measured the change
    # dmut[i]=mu[i+1]-mu[i] 
    dmut <- c(dmu[ , t], 0) 
    
    # Xb and Pb are taken from what was stored for each trait from the previous window
    Xb=XX[ ((t-1)*n.cols+1) : (t*n.cols) ]
    Pb=PP[ ((t-1)*n.cols+1) : (t*n.cols) , 1:n.cols ]
    
    # The Kalman filter is run, with the given value of rho 
    output<-KB(dmuLt,dmut,rho,C,L.dyn,Pb,Xb)
    dmuKt <- output [[1]]
    
    # Now we calculate the difference between the prediction using the method with a given rho 
    # inside the window, against the observed change inside the window.
    # The first element of dmu and the last element of dmuKt are removed so that the
    # difference is compatible
    dmurt <- c(dmu[ c(-1), t])
    dmuKt <- dmuKt[-length(dmuKt)]
    
    eK <- sum((dmuKt-dmurt)**2)
    error.rho[t] <- eK
    
  }
  return(error.rho)  
}
