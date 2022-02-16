source("predict.R")
source("find_rho.r")
source("Kalman_Filter.r")

# This is the main script file. 

# Global variables are defined#############################################################

nmax <- 30       # Total number of generations in the experiment
n.traits <- 5    # Number of traits
L <- 5          # Maximum number of generations in the sliding time window

C <- rbind(c(1, -1), c(1, 0)) # This is the matrix that relates the
                              # measurements to the states. See Appendix A.

n.rows <- dim(C)[1]       # Dimensions of matrix C
n.cols <- dim(C)[2]

L.0 = 2  # Initial generation. We need at least 1 recorded change, so we
         # start in generation 2.

L.dyn = L.0   # Dynamic number of generations in the sliding time window. The 
              # maximum number is L. When there are less than L data points available
              # the window is as large as possible. This allows to correct predictions
              # even at the beginning of the experiment.

NN=1  # For the first NN generations, there is not enough data in the time
      # window to find rho. Then, we use rho = 0 and the method gives the same
      # predictions as the breeder's equation. 

# Declare the matrices of covariances and states.

XX<-matrix(,nrow=n.traits*n.cols,ncol=1)        # XX stores the estimates at the start of
                                                # sliding time window
PP<-matrix(,nrow=n.traits*n.cols, ncol=n.cols)  # PP stores the covariance matrix of the error
                                                # at the start of the window

# Declare the matrices that store the predicted and observed changes

dmuKK <- matrix(,ncol=n.traits, nrow= nmax) # Store the predicted change using the new method
dmuLL <- matrix(,ncol=n.traits, nrow= nmax) # Store the predicted change using the breeder's equation
dmuRR <- matrix(,ncol=n.traits, nrow= nmax) # Store the observed change


###########################################################################################

###Load data###############################################################################

s <- as.matrix(read.table(paste("s.txt",sep="")))   # Selection differential
mu <- as.matrix(read.table(paste("mu.txt",sep=""))) # Mean of the traits
G <- as.matrix(read.table(paste("G.txt",sep="")))   # G-matrix at generation 1
P <- as.matrix(read.table(paste("P.txt",sep="")))   # P-matrix at generation 1

###########################################################################################

for (i in L.0:nmax){     # The method start in generation i
  
  # Measured changes in the trait inside time window
  dmu <- matrix(mu[ (i-L.dyn+2) : (i),] - mu[(i-L.dyn+1) : (i-1),], ncol=n.traits)
  
  # Breeder's predictions inside the time window
  Fs <- G %*% solve(P) %*%t (s[( i - L.dyn + 1 ) : i,])  
  
  # Find parameter rho in the time window (see Part III)
  
  if (i <= L.0+NN){
    best.rho.t <- c(0,0,0,0,0)
  }else{
    # Declare dummy arrays to store values during the identification of parameter rho (Part III)
    error.rho <- rep(0,n.traits)
    store.error.rho <- rep(NA, n.traits+1) 
    
    # Explore values of rho. The array store.error.rho stores the prediction error
    # for each trait using a given rho
  
    for (rho in c(0,c(1) %o% 10^(-5:3))){
      error.rho <- test.rho() # This runs the Kalman filter inside the window, and
                              # returns the prediction error for each trait, for a given rho
      store.error.rho<-rbind(store.error.rho, c(rho,error.rho))
    }  
  
    # For each trait we keep the value of rho that minimizes the prediction error in the
    # array best.rho.t[t]
    
    store.error.rho<-store.error.rho[-1,]
    best.rho.t <- rep(NA,n.traits)
    for (t in 1:n.traits){
      best.rho.t[t] <- store.error.rho[which.min(store.error.rho[,t+1]),1]      
    }
  }
  
  # Now obtain the prediction using the best rho for each trait
  predict()
  
  
  if(L.dyn<L){L.dyn=L.dyn+1}
  
  
}


# Obtain the observed change for all generations
for(t in 1:n.traits){
  L=nmax+1
  dmurtt<-mu[2:L,]-mu[1:(L-1),]
  dmuRR [, t] <- dmurtt[,t]
}


# Plot the observed and predicted changes. 

par(mfrow=c(3,2),mar=c(2,2,0,0)+2, mgp=c(2,1,0))
for(trait in 1:n.traits){
  y.max <- max(c(dmuKK[,trait],dmuLL[,trait],dmuRR[,trait]),na.rm=T)
  y.min <- min(c(dmuKK[,trait],dmuLL[,trait],dmuRR[,trait]),na.rm=T)
  
  plot(dmuRR[,trait],type="p",pch=16,xlim=c(1,30),ylim=c(y.min,y.max),
       ylab=bquote(Delta [.(trait)]),
       lwd=1,xlab="Generations")
  lines(dmuRR[,trait],type="l")
  
  lines(dmuLL[,trait],col="blue")
  points(dmuLL[,trait],col="blue",pch=16)
  
  lines(dmuKK[,trait],col="red")
  points(dmuKK[,trait],col="red",pch=16)
  
  abline(h=0)
}