library(smooth)
library(stringr)

invSim <- function(csl=0.9, R=3, L=2,
                   fcst.model=c("ETS","ETS-SHR","ARIMA")) {
  
  y <- sim.ssarima(orders = list(ar=1,i=0,ma=0), obs = 500, AR = 0.7)
  y <- 100 + window(y$data, start = 201)
  y <- ts(y, start = 1)
  # plot(y)
  
  csl <- csl
  R <- R
  L <- L
  fcst.model <- fcst.model
  
  yTrain <- window(y, end = 50)
  yTest <- window(y, start = 51)
  
  H <- R + L
  n <- length(yTest)
  
  if (fcst.model == "ETS") {
    fit <- adam(yTrain, loss = "MSE")
    fcst <- forecast(fit, h = H, interval = "empirical", level = csl)
  } else if (fcst.model == "ARIMA") {
    fit <- ssarima(yTrain)
    fcst <- forecast(fit, h = H, interval = "empirical", level = csl)
  } else if (fcst.model == "ETS-SHR") {
    fit0 <- adam(yTrain, loss = "MSE")
    if (str_length(fit0$model)==8) {
      modelStruc <- substr(fit0$model, 5, 7)
    } else if (str_length(fit0$model)==9) {
      modelStruc <- substr(fit0$model, 5, 8)
    }
    fit <- adam(yTrain, model = modelStruc, loss = "RIDGE", lambda = 0.2)
    fcst <- forecast(fit, h = H, interval = "empirical", level = csl)
  } else if (is.null(fcst.model)) {
    fit <- adam(yTrain, loss = "MSE")
    fcst <- forecast(fit, h = H, interval = "empirical", level = csl)
  }
  
  simDat <- array(NA,c(n,7))
  colnames(simDat) <- c("R","Demand","S","InvPos","NetStock","StockOnHnd","Orders")
  
  initInv <- sum(fcst$upper[1:(H-1)])
  simDat[1,] <- c(0, yTest[1], initInv, initInv, initInv, initInv, 0)
  
  r <- rep(0,R)
  r[1] <- 1
  simDat[2:n,1] <- rep(r,ceiling((n-1)/R))[1:(n-1)]
  
  S <- sum(fcst$upper)
  simDat[2:n,"S"] <- S
  
  for (i in 2:n){
    
    # Store demand for this period
    # simDat[i,"Demand"] <- d[i]
    simDat[i,"Demand"] <- yTest[i]
    
    # Calculate order up to based on forecasting model
    # Cumulative demand for each R+L-1 period
    # S <- sum(fcst$upper)
    # simDat[i,"S"] <- S
    
    if (simDat[i,"R"] == 1){ # Is it a review period?
      
      # Set inventory position to S - demand
      simDat[i,"InvPos"] <- simDat[i,"S"] - simDat[i,"Demand"]
      # See if orders exist and retrieve them
      if (i-1 >= L){
        ordersTemp <- simDat[i-L,"Orders"]  
      } else {
        ordersTemp <- 0
      }
      # Calculate new net stock
      simDat[i,"NetStock"] <- simDat[i-1,"NetStock"] - simDat[i,"Demand"] + ordersTemp
      # Calculate stock on hand, no backorders
      if (simDat[i,"NetStock"] < 0){
        simDat[i,"StockOnHnd"] <- 0  
      } else {
        simDat[i,"StockOnHnd"] <- simDat[i,"NetStock"]
      }
      # Calculate new orders (remenber to add demand to inv position)
      ordersTemp <- simDat[i,"S"] - simDat[i,"StockOnHnd"]
      if (ordersTemp < 0){ordersTemp <- 0}
      simDat[i,"Orders"] <- ordersTemp
      
    } else { # Not in a review period
      
      # Set inventory position to S
      simDat[i,"InvPos"] <- simDat[i-1,"InvPos"] - simDat[i,"Demand"]
      # See if orders exist and retrieve them
      if (i-1 >= L){
        ordersTemp <- simDat[i-L,"Orders"]  
      } else {
        ordersTemp <- 0
      }
      # Calculate new net stock
      simDat[i,"NetStock"] <- simDat[i-1,"NetStock"] - simDat[i,"Demand"] + ordersTemp
      # Calculate stock on hand
      if (simDat[i,"NetStock"] < 0){
        simDat[i,"StockOnHnd"] <- 0  
      } else {
        simDat[i,"StockOnHnd"] <- simDat[i,5]
      }
      # Calculate new orders
      simDat[i,"Orders"] <- 0
      
    }
  }
  
  k <- simDat[51:250,"StockOnHnd"] == 0
  kPad <- c(k,rep(0,(R*ceiling(length(k)/R) - length(k))))
  kPad <- colSums(matrix(kPad,nrow=R,byrow=FALSE))
  kPad[kPad > 0] <- 1
  
  result <- list(parameter = fit$B,
                 inventory = simDat,
                 metric = c(csl, 1 - mean(kPad), 1 - sum(simDat[51:250,"StockOnHnd"]==0)/(n-1)))
  
  return(result)
  
}

iter <- 500
set.seed(020192)
InvSim_ETS_CSL90 <- vector("list", iter)
for (i in 1:iter) {
  InvSim_ETS_CSL90[[i]] <- invSim(csl = 0.9, R = 1, L = 3, fcst.model = "ETS")
}
set.seed(020192)
InvSim_ETSshr_CSL90 <- vector("list", iter)
for (i in 1:iter) {
  InvSim_ETSshr_CSL90[[i]] <- invSim(csl = 0.9, R = 1, L = 3, fcst.model = "ETS-SHR")
}
set.seed(020192)
InvSim_ARIMA_CSL90 <- vector("list", iter)
for (i in 1:iter) {
  InvSim_ARIMA_CSL90[[i]] <- invSim(csl = 0.9, R = 1, L = 3, fcst.model = "ARIMA")
}

rowMeans(sapply(InvSim_ETS_CSL90, function(x) x$metric))
rowMeans(sapply(InvSim_ETSshr_CSL90, function(x) x$metric))
rowMeans(sapply(InvSim_ARIMA_CSL90, function(x) x$metric))

boxplot(cbind(colMeans(t(sapply(InvSim_ETS_CSL90, function(x) x$inventory[,"StockOnHnd"]))),
              colMeans(t(sapply(InvSim_ETSshr_CSL90, function(x) x$inventory[,"StockOnHnd"]))),
              colMeans(t(sapply(InvSim_ARIMA_CSL90, function(x) x$inventory[,"StockOnHnd"])))))

# plot(simDat[,"S"],type="l",ylim=c(min(simDat[,"NetStock"])-200,max(simDat[,"S"])*1.1))
# lines(simDat[,"StockOnHnd"],col="red")
# lines(simDat[,"NetStock"],col="blue")
# lines(simDat[,"Demand"],col="green")
# # lines(d,col="magenta")
# abline(h=0,col="grey")
# 
# k <- simDat[-1,"StockOnHnd"] == 0
# kPad <- c(k,rep(0,(R*ceiling(length(k)/R) - length(k))))
# kPad <- colSums(matrix(kPad,nrow=R,byrow=FALSE))
# kPad[kPad > 0] <- 1
# 
# c(csl, 1 - mean(kPad), 1 - sum(simDat[-1,"StockOnHnd"]==0)/(n-1))
