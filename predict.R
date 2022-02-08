# This function runs the Kalman filter with the best rho

predict<-function(){
  
  for (t in 1:n.traits) {
    # dmuLt stores the prediction of the breeder's equation for trait t
    dmuLt <- Fs[t, ]
    
    # dmut stores the measured changes in trait means for trait t.
    # there is an extra zero because we have not measured the change
    # dmut[i]=mu[i+1]-mu[i] 
    dmut <- c(dmu[ , t], 0) 

    if (i==L.0){
      # The initial conditions for the state variables are the breeder's prediction
      # and bias=0
      Xb=c(dmuLt[L.0],0)
      # The initial conditions for the covariance of the error is the zero matrix
      Pb=diag(n.cols)*best.rho.t[t]
    } else {
      # For the rest of the generations, Xb and Pb are taken from what was stored for each trait
      Pb=PP[ ((t-1)*n.cols+1) : (t*n.cols) , 1:n.cols ]
      Xb=XX[ ((t-1)*n.cols+1) : (t*n.cols) ]
    }
    
    # Run the Kalman filter with the best rho for each trait
    output<-KB(dmuLt,dmut,best.rho.t[t],C,L.dyn,Pb,Xb)
    dmuKt <- output [[1]]
    Pb <- output [[2]]
    Xb <- output [[3]]
    
    # The last element inside the window has the prediction of interest,
    # which is the change in trait mean from i to i+1
    dmuKK[i,t] <<- dmuKt[length(dmuKt)]
    dmuLL[i,t] <<- dmuLt[length(dmuLt)]
    
    # The covariance of the error and the value of the states is 
    # stored for the next generation
    PP[ ((t-1)*n.cols+1) : (t*n.cols) , 1:n.cols ] <<- Pb 
    XX[ ((t-1)*n.cols+1) : (t*n.cols) , 1 ] <<- Xb
    
  }
}