predict<-function(){
  #Ahora se hace el Kalman con el mejor q con objectivo de predicción
  for (t in 1:nT) {
    
    dmut <- c(dmu[(j-nw+1):(j-1), t], 0) ##cambio raw medido; dmu[j]=mu[j+1]-mu[j] por eso no se incluye 
    dmuLt <- Fs[ t, (j-nw+1) : j ]
    
    if (j==nw.0){
      Xb=c(dmuLt[nf],rep(0,nc-1))
      Pb=diag(nc)*0
    } else {
      Pb=PP[ ((t-1)*nc+1) : (t*nc) , 1:nc ]
      Xb=XX[ ((t-1)*nc+1) : (t*nc) ]
    }
    
    output<-KB(dmuLt,dmut,best.rho.t[t],C,nw,Pb,Xb)
    dmuKt <- output [[1]]
    Pb <- output [[2]]
    Xb <- output [[3]]
    
    #El último elemento tiene la prediccion i -> i+1
    dmuKK[j,t] <<- dmuKt[length(dmuKt)]
    dmuLL[j,t] <<- dmuLt[length(dmuLt)]
    
    PP[ ((t-1)*nc+1) : (t*nc) , 1:nc ] <<- Pb 
    XX[ ((t-1)*nc+1) : (t*nc) , 1 ] <<- Xb
    
  }
}