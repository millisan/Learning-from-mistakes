test.rho<-function(){
  for (t in 1:nT) {
    #los deltas, le agregamos un zero al final porque no sabemos el deltamu real en j
    dmut <- c(dmu[ (j-nw+1) : (j-1), t], 0) ##cambio raw medido; dmu[j]=mu[j+1]-mu[j] por eso no se incluye 
    #breeder's hasta j, ese si tiene valor en j
    dmuLt <- Fs[t, (j-nw+1) : j]
    
    #inicializar para la primera genetacion
    if (j==nw.0){
      #Xb es el vector de estados. para inciiarlizarlo tenemos Lande y dsp xero
      Xb <- c(dmuLt[nf],rep(0,nc-1))
      #en la primera generacion obligamos que q sea cero porque no hay data 
      rho=0
      #inicializamos la covariance del error de estimacion de los estados (Xb-X). confiamos en lande
      Pb <- diag(nc)*0
      return() #we don't have enough data to calculate anything yet
    }else{
      #en el resto de generaciones usamos la informacion de la ventana que se uso en la generacion pasada
      Xb=XX[ ((t-1)*nc+1) : (t*nc) ]
      Pb=PP[ ((t-1)*nc+1) : (t*nc) , 1:nc ]
    }
    
    output<-KB(dmuLt,dmut,rho,C,nw,Pb,Xb)
    dmuKt <- output [[1]]
    Pb <- output [[2]]
    Xb <- output [[3]]
    
    #en el ala suavizamos el dmurt
    #j-nw:inicio de ventana, pero en ventana se arranca en nf-> arrancamos en j-nw+nf
    dmurt <- c(dmu[ (j-nw+nf) : (j-1), t])
    dmuKt <- dmuKt[-length(dmuKt)] #el ultimo no se puede comparar porque no tenemos el delta real
    
    dmuKt<-as.matrix(dmuKt)
    dmurt<-as.matrix(dmurt)
    
    #cuando j=nw.0 o nw.0+1
    eK <- t(dmuKt-dmurt) %*% (dmuKt-dmurt)
    error.rho[t] <- eK
    
  }
  return(error.rho)  
}
