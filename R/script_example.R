n <- 1200   # Sample size

# Generate some demand
d <- rnorm(n,mean=1000,sd=100)
d[d<0] <- 0

# Inventroy settings
R <- 5
L <- 3
csl <- 0.5

# Forecast model
mu <- mean(d)
sigma <- sd(d)

# The data are stored as
# Review period, demand, S, inventory position, net stock, stock on hand, orders

simDat <- array(NA,c(n,7))
colnames(simDat) <- c("R","Demand","S","InvPos","NetStock","StockOnHnd","Orders")
# Initial values
initInv <- (L+R-1)*quantile(d,csl) # (L+R)*mu+qnorm(csl)*sigma
simDat[1,] <- c(0,d[1],initInv,initInv,initInv,initInv,0)

# Add all review indicators
r <- rep(0,R)
r[1] <- 1
simDat[2:n,1] <- rep(r,ceiling((n-1)/R))[1:(n-1)]

for (i in 2:n){
  
  # Store demand for this period
  simDat[i,2] <- d[i]
  
  # Calculate order up to based on forecasting model
  S <- quantile(colSums(matrix(d,nrow=R+L-1)),csl) # (L+R-1)*mu + qnorm(csl)*sigma
  simDat[i,3] <- S
  
  if (simDat[i,1] == 1){ # Is it a review period?
    
    # Set inventory position to S - demand
    simDat[i,4] <- simDat[i,3] - simDat[i,2]
    # See if orders exist and retrieve them
    if (i-1 >= L){
      ordersTemp <- simDat[i-L,7]  
    } else {
      ordersTemp <- 0
    }
    # Calculate new net stock
    simDat[i,5] <- simDat[i-1,5] - simDat[i,2] + ordersTemp
    # Calculate stock on hand
    if (simDat[i,5] < 0){
      simDat[i,6] <- 0  
    } else {
      simDat[i,6] <- simDat[i,5]
    }
    # Calculate new orders (remenber to add demand to inv position)
    ordersTemp <- simDat[i,3] - simDat[i,6]
    if (ordersTemp < 0){ordersTemp <- 0}
    simDat[i,7] <- ordersTemp
  
  } else { # Not in a review period
    
    # Set inventory position to S
    simDat[i,4] <- simDat[i-1,4] - simDat[i,2]
    # See if orders exist and retrieve them
    if (i-1 >= L){
      ordersTemp <- simDat[i-L,7]  
    } else {
      ordersTemp <- 0
    }
    # Calculate new net stock
    simDat[i,5] <- simDat[i-1,5] - simDat[i,2] + ordersTemp
    # Calculate stock on hand
    if (simDat[i,5] < 0){
      simDat[i,6] <- 0  
    } else {
      simDat[i,6] <- simDat[i,5]
    }
    # Calculate new orders
    simDat[i,7] <- 0
      
  }
}

simDat
plot(simDat[,3],type="l",ylim=c(-1000,max(simDat[,3])*1.1))
lines(simDat[,6],col="red")
lines(simDat[,5],col="blue")
lines(simDat[,2],col="green")
# lines(d,col="magenta")
abline(h=0,col="grey")

k <- simDat[-1,6] == 0
kPad <- c(k,rep(0,(R*ceiling(length(k)/R) - length(k))))
kPad <- colSums(matrix(kPad,nrow=R,byrow=FALSE))
kPad[kPad > 0] <- 1


c(csl, 1 - mean(kPad), 1 - sum(simDat[-1,6]==0)/(n-1))
