source("predict.R")
source("find_rho.r")
source("Kalman_Filter.r")




# Global variables are defined

nmax=30  # Total number of generations in the experiment
n.traits=5   # Number of traits
L <- 5 # Maximum number of generations in the sliding time window

C <- rbind(c(1, -1), c(1, 0)) # This is the matrix that relates the
                              # measurements to the states
n.rows<-dim(C)[1] 
n.cols<-dim(C)[2]

# Initialize the matrices of covariances
XX<-matrix(,nrow=n.traits*n.cols,ncol=1)        # XX stores the estimates at the start of
                                                # sliding time window
PP<-matrix(,nrow=n.traits*n.cols, ncol=n.cols)  # PP stores the covariance matrix of the error
                                                # at the start of the window

#inicializar matrices que guardan los cambios  
dmuKK <- matrix(,ncol=n.traits, nrow= nmax)
dmuLL <- matrix(,ncol=n.traits, nrow= nmax)
dmuRR <- matrix(,ncol=n.traits, nrow= nmax)
#Definir primera generacion en la que empieza a correr, 
#no se puede la 1 porque se necesita registro  . 
#no se puede arrancar antes de tener nf datos que son las salidas. 
#en este caso necesitas la pasada, nf = 2
#cambiar nombre a nw.0

L.0 = n.rows # Initial number of generations in the sliding time window.
L.dyn = L.0 # Dynamic number of generations in the sliding time window. The 
            # maximum number is L. When there are less than L data points available
            # the window is as large as possible. This allows to correct predictions
            # even at the beginning of the experiment.

############################################################







####Cargar datos#############################################3


s <- as.matrix(read.table(paste("s_",id,".txt",sep="")))
mu <- as.matrix(read.table(paste("mu_",id,".txt",sep="")))
#PONER UNA G CORTA!!
G <- as.matrix(read.table(paste("G_",id,".txt",sep="")))
P <- as.matrix(read.table(paste("P_",id,".txt",sep="")))
Gi= G[ 1:5 , 1:5 ]    # Extraccion de G(i)
Pi= P[ 1:5 , 1:5 ]    # Extraccion de P(i)


#PONER A MU COMO /100 EN EL ARCHIVO
mu=mu/100
dmu.store=mu[2:(nmax+1),]-mu[1:nmax,]
##################################################3



######Breeder's
Fi=Gi %*% solve(Pi)
Fs.store=Fi%*%t(s)
################################




#Matrix para iniciar la ventana movil  
R<-NA

for (j in L.0:nmax){     # se establece la ventana deslizante
  dmu <- matrix(mu[2:j,]-mu[1:(j-1),],ncol=n.traits)
  Fs <- Fi%*%t(s[1:j,])
  
  store.error.rho <- rep(NA, n.traits+1)
  error.rho <- rep(0,n.traits)
  
  for (rho in c(0,c(1) %o% 10^(-5:3))){
    error.rho <- test.rho() 
    store.error.rho<-rbind(store.error.rho, c(rho,error.rho))
  }  
  
  
  #Quedarse con el rho que minimiza el error
  
  store.error.rho<-store.error.rho[-1,]
  best.rho.t <- rep(NA,n.traits)
  
  for (t in 1:n.traits){
    best.rho.t[t] <- store.error.rho[which.min(store.error.rho[,t+1]),1]      
  }
  
  predict()
  
  
  if(L.dyn<L){L.dyn=L.dyn+1}
  
  
}







#el delta real
for(t in 1:n.traits){
  L=nmax+1
  dmurtt<-mu[2:L,]-mu[1:(L-1),]
  dmuRR [, t] <- dmurtt[,t]
}

#MSE por generacion

euc.dist<-function(v){sqrt(sum(v**2))}

eB <- apply(dmuLL-dmuRR, MAR=1, FUN= euc.dist) / apply(dmuRR, MAR=1, FUN= euc.dist)
eK <- apply(dmuKK-dmuRR, MAR=1, FUN= euc.dist) / apply(dmuRR, MAR=1, FUN= euc.dist)

#RMSE total
RMSE.K <- norm((dmuKK-dmuRR)[-c(1),],type = "F") / norm(dmuRR[-c(1),],type = "F")
RMSE.L  <- norm((dmuLL-dmuRR)[-c(1),],type = "F") / norm(dmuRR[-c(1),],type = "F")



#Plot



par(mfrow=c(3,2),mar=c(2,2,0,0)+2, mgp=c(2,1,0))
for(trait in 1:5){
  
  y.max <- max(c(dmuKK[,trait],dmuLL[,trait],dmuRR[,trait]),na.rm=T)
  y.min <- min(c(dmuKK[,trait],dmuLL[,trait],dmuRR[,trait]),na.rm=T)
  
  plot(dmuRR[,trait],type="p",pch=16,xlim=c(1,30),ylim=c(y.min,y.max),
       ylab=bquote(Delta [.(trait)]),
       lwd=1,xlab="Generations")
  
  
  lines(dmuRR[,trait],type="l",ylim=c(-.1,.1))
  lines(dmuLL[,trait],col="blue",pch=16)
  points(dmuLL[,trait],col="blue",pch=16)
  
  lines(dmuKK[,trait],col="red",pch=16)
  points(dmuKK[,trait],col="red",pch=16)
  
  abline(h=0)
  
  
}
