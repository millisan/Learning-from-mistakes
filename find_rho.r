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
    dmuKt <- output [[1]] # dmuKt stores the predictions using the method inside the window
    
    # Now we calculate the difference between the predictions using the method with a given rho 
    # (i.e. dmuKt), against the observed changes inside the window (i.e. dmurt). 
    # For the teeth experiments, dmurt is directly the measured change in
    # trait means. For the fly wing experiments, dmurt is obtained by making a linear 
    # regression of the mean of the trait inside the window.
    # The first element of dmurt is removed from the comparison because there is no
    # prediction for the first generation in the window. The last element of dmuKt is removed
    # because it corresponds to generation i, for which there is no observed change to compare.
    
    dmurt <- c(dmu[ c(-1), t])
    dmuKt <- dmuKt[-length(dmuKt)]
    
    eK <- sum((dmuKt-dmurt)**2)
    error.rho[t] <- eK
    
  }
  return(error.rho)  
}
