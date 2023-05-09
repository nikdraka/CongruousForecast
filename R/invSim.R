

# Congruent Forecasts -----------------------------------------------------

SumFcst <- vector("list", 500)
for (i in 1:500) {
  
  y <- sim.ssarima(orders = list(ar=1,i=0,ma=0), obs = 500, AR = 0.6)
  y <- 100 + window(y$data, start = 201)
  y <- ts(y, start = 1)
  
  train <- 50
  test <- 250
  
  horizon <- 5
  origins <- test/horizon
  
  temp <- matrix(NA, ncol = 4, nrow = test)
  colnames(temp) <- c("ETS-MSE", "ETS-MAE", "ETS-SHR", "ARIMA")
  for (j in 1:origins) {
    yTrain <- window(y, end=50+(j-1)*horizon)
    
    range <- (train+1+horizon*(j-1)):(train+horizon*(j))
    
    fitETSmse <- adam(yTrain, loss = "MSE")
    fcstETSmse <- forecast(fitETSmse, h = horizon, interval = "empirical", level = 0.95)
    
    if (nchar(fitETSmse$model)==8) {
      modelStruc <- substr(fitETSmse$model, 5, 7)
    } else if (nchar(fitETSmse$model) == 9) {
      modelStruc <- substr(fitETSmse$model, 5, 8)
    }
    
    fitETSshr <- adam(yTrain, model = modelStruc, loss = "RIDGE", lambda = 0.3)
    fcstETSshr <- forecast(fitETSshr, h = horizon, interval = "empirical", level = 0.95)
    
    fitETSmae <- adam(yTrain, loss = "MAE")
    fcstETSmae <- forecast(fitETSmae, h = horizon, interval = "empirical", level = 0.95)
    
    fitARIMA <- ssarima(yTrain)
    fcstARIMA <- forecast(fitARIMA, h = horizon, interval = "empirical", level = 0.95)
    
    temp[range-train,1] <- sum(fcstETSmse$upper)
    temp[range-train,2] <- sum(fcstETSmae$upper)
    temp[range-train,3] <- sum(fcstETSshr$upper)
    temp[range-train,4] <- sum(fcstARIMA$upper)
    
  }
  
  # plot(temp[,1], type = "l", ylim = range(temp))
  # lines(temp[,2], type = "l", col = "red")
  # lines(temp[,3], type = "l", col = "blue")
  
  # meanSumFcst[i,] <- apply(temp, 2, mean)
  # sdSumFcst[i,] <- apply(temp, 2, sd)
  
  SumFcst[[i]] <- temp
  
}

save(SumFcst, file = "SumFcst.RData")

plot(ecdf(rowMeans(sapply(SumFcst, function(x) x[,1]))), xlim = c(640, 700))
plot(ecdf(rowMeans(sapply(SumFcst, function(x) x[,2]))), xlim = c(640, 700))
plot(ecdf(rowMeans(sapply(SumFcst, function(x) x[,3]))), xlim = c(640, 700))
plot(ecdf(rowMeans(sapply(SumFcst, function(x) x[,4]))), xlim = c(640, 700))

hist(c(sapply(SumFcst, function(x) x[,1])), breaks = 100, xlim = c(300, 1200), ylim = c(0, 8000),
     main = "ETS-MSE", xlab = "95% Quantile Forecasts")
hist(c(sapply(SumFcst, function(x) x[,2])), breaks = 100, xlim = c(300, 1200), ylim = c(0, 8000),
     main = "ETS-MAE", xlab = "95% Quantile Forecasts")
hist(c(sapply(SumFcst, function(x) x[,3])), breaks = 100, xlim = c(300, 1200), ylim = c(0, 8000),
     main = "ETS-SHR", xlab = "95% Quantile Forecasts")
hist(c(sapply(SumFcst, function(x) x[,4])), breaks = 100, xlim = c(300, 1200), ylim = c(0, 8000),
     main = "AutoARIMA", xlab = "95% Quantile Forecasts")

print(matrix(c(mean(sapply(SumFcst, function(x) x[,1])), sd(sapply(SumFcst, function(x) x[,1])),
               mean(sapply(SumFcst, function(x) x[,2])), sd(sapply(SumFcst, function(x) x[,2])),
               mean(sapply(SumFcst, function(x) x[,3])), sd(sapply(SumFcst, function(x) x[,3])),
               mean(sapply(SumFcst, function(x) x[,4])), sd(sapply(SumFcst, function(x) x[,4]))),
             ncol = 2, byrow = TRUE,
             dimnames = list(c("ETS-MSE", "ETS-MAE", "ETS-SHR", "AutoARIMA"),c("Mean", "SD"))))


# Simple experiments ------------------------------------------------------

library(smooth)
library(stringr)

invSim <- function(csl=0.9, R=3, L=2,
                   fcst.model=c("ETS-MSE", "ETS-MAE", "ETS-SHR", "ARIMA")) {
  
  # Generate AR(1) time series with alpha=0.7, e~N(0,1)
  y <- sim.ssarima(orders = list(ar=1,i=0,ma=0), obs = 500, AR = 0.7)
  y <- 100 + window(y$data, start = 201)
  y <- ts(y, start = 1)
  
  # Make sure that we dont have negative numbers
  while (sum(y < 0) != 0) {
    y <- sim.ssarima(orders = list(ar=1,i=0,ma=0), obs = 500, AR = 0.7)
    y <- 100 + window(y$data, start = 201)
    y <- ts(y, start = 1)
  }
  
  # Set up CSL, Review period, and Lead time, and what forecasting model we use
  csl <- csl
  R <- R
  L <- L
  fcst.model <- fcst.model
  
  # Set up the training and the test set
  # I do not employ any origin because we iterate it for 500 repetitions
  # Out of 300 obs, we use 50 as training set
  # 250 obs as test set
  yTrain <- window(y, end = 50)
  yTest <- window(y, start = 51)
  
  # H = 5, because R = 3 and L = 2
  # Produce 1-5 step ahead forecasts
  H <- R + L
  n <- length(yTest)
  
  # Choose the forecasting model
  # Produce forecasts on the mean and the quantile
  if (fcst.model == "ETS-MSE") {
    fit <- adam(yTrain, loss = "MSE")
    fcst <- forecast(fit, h = H, interval = "empirical", level = csl)
  } else if (fcst.model == "ARIMA") {
    fit <- msarima(yTrain)
    fcst <- forecast(fit, h = H, interval = "empirical", level = csl)
  } else if (fcst.model == "ETS-SHR") {
    fit0 <- adam(yTrain, loss = "MSE")
    if (nchar(fit0$model)==8) {
      modelStruc <- substr(fit0$model, 5, 7)
    } else if (nchar(fit0$model)==9) {
      modelStruc <- substr(fit0$model, 5, 8)
    }
    fit <- adam(yTrain, model = modelStruc, loss = "RIDGE", lambda = 0.2)
    fcst <- forecast(fit, h = H, interval = "empirical", level = csl)
  } else if (fcst.model == "ETS-MAE") {
    fit <- adam(yTrain, loss = "MAE")
    fcst <- forecast(fit, h = H, interval = "empirical", level = csl)
  } else if (is.null(fcst.model)) {
    fit <- adam(yTrain, loss = "MSE")
    fcst <- forecast(fit, h = H, interval = "empirical", level = csl)
  }
  
  # Measure forecast accuracy
  scaled <- mean(abs(diff(yTrain)))
  accMatrix <- matrix(c(abs(ME(yTest[1:H], fcst$mean, na.rm = TRUE)),
                        sqrt(MSE(yTest[1:H], fcst$mean, na.rm = TRUE)),
                        pinball(yTest[1:H], fcst$upper, level = 0.95, na.rm = TRUE),
                        abs(sCE(yTest[1:5], fcst$mean, scaled))),
                      ncol = 4, nrow = 1,
                      dimnames = list(fcst.model, c("AME", "RMSE", "PinballUp", "|sCE|")))
  
  # Inventory simulation!
  simDat <- array(NA,c(n,7))
  colnames(simDat) <- c("R","Demand","S","InvPos","NetStock","StockOnHnd","Orders")
  
  initInv <- sum(fcst$upper[1:(H-1)]) # followed your line
  simDat[1,] <- c(0, yTest[1], initInv, initInv, initInv, initInv, 0)
  
  # Set up review periods
  r <- rep(0,R)
  r[1] <- 1
  simDat[2:n,1] <- rep(r,ceiling((n-1)/R))[1:(n-1)]
  
  ### CAUTION!!!!! NEED ATTENTION!!!!
  # Set up the safety stock level the same as sum(fcst$upper) across the test set
  # Is it ok? I am not sure it's a good approach. How should we incorporate rolling origins?
  # How often does a company change their safety stock? Every review period?
  S <- sum(fcst$upper)
  simDat[2:n,"S"] <- S
  
  for (i in 2:n){
    
    # Store demand for this period
    simDat[i,"Demand"] <- yTest[i]
    
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
                 metric = c(csl, 1 - mean(kPad), 1 - sum(simDat[51:250,"StockOnHnd"]==0)/(n-1)),
                 accuracy = accMatrix)
  
  return(result)
  
}

set.seed(020192)
abc <- invSim(csl = 0.9, R = 3, L = 2, fcst.model = "ARIMA")

iter <- 500
set.seed(020192)
InvSim_ETSmse_CSL90 <- vector("list", iter)
for (i in 1:iter) {
  InvSim_ETSmse_CSL90[[i]] <- invSim(csl = 0.9, R = 3, L = 2, fcst.model = "ETS-MSE")
}
set.seed(020192)
InvSim_ETSmae_CSL90 <- vector("list", iter)
for (i in 1:iter) {
  InvSim_ETSmae_CSL90[[i]] <- invSim(csl = 0.9, R = 3, L = 2, fcst.model = "ETS-MAE")
}
set.seed(020192)
InvSim_ETSshr_CSL90 <- vector("list", iter)
for (i in 1:iter) {
  InvSim_ETSshr_CSL90[[i]] <- invSim(csl = 0.9, R = 3, L = 2, fcst.model = "ETS-SHR")
}
set.seed(020192)
InvSim_ARIMA_CSL90 <- vector("list", iter)
for (i in 1:iter) {
  InvSim_ARIMA_CSL90[[i]] <- invSim(csl = 0.9, R = 3, L = 2, fcst.model = "ARIMA")
}

save(InvSim_ETSmse_CSL90, file = "InvSim_ETSmse_CSL90.RData")
save(InvSim_ETSmae_CSL90, file = "InvSim_ETSmae_CSL90.RData")
save(InvSim_ETSshr_CSL90, file = "InvSim_ETSshr_CSL90.RData")
save(InvSim_ARIMA_CSL90, file = "InvSim_ARIMA_CSL90.RData")

InvSim_metric_CSL90 <- matrix(c(rowMeans(sapply(InvSim_ETSmse_CSL90, function(x) x$metric)),
                                rowMeans(sapply(InvSim_ETSmae_CSL90, function(x) x$metric)),
                                rowMeans(sapply(InvSim_ETSshr_CSL90, function(x) x$metric)),
                                rowMeans(sapply(InvSim_ARIMA_CSL90, function(x) x$metric))),
                              ncol = 3, byrow = TRUE, 
                              dimnames = list(c("ETS-MSE", "ETS-MAE", "ETS-SHR", "ARIMA"), c("CSL", "EmpCSL", "FillRate")))
round(InvSim_metric_CSL90*100, 2)

InvSim_StockOnHnd_CSL90 <- cbind(colMeans(t(sapply(InvSim_ETSmse_CSL90, function(x) x$inventory[,"StockOnHnd"]))),
                                 colMeans(t(sapply(InvSim_ETSmae_CSL90, function(x) x$inventory[,"StockOnHnd"]))),
                                 colMeans(t(sapply(InvSim_ETSshr_CSL90, function(x) x$inventory[,"StockOnHnd"]))),
                                 colMeans(t(sapply(InvSim_ARIMA_CSL90, function(x) x$inventory[,"StockOnHnd"]))))

boxplot(InvSim_StockOnHnd_CSL90, col = "white", names = c("ETS-MSE", "ETS-MAE", "ETS-SHR", "ARIMA"), las = 2,
        main = "Stock on hand")
lines(colMeans(InvSim_StockOnHnd_CSL90), type = "o", col = "red", pch = 20, cex = 1.5)

removeZero <- function(x) {
  x <- x[x!=0]
  return(x)
}

plot(rowMeans(sapply(InvSim_ETSmse_CSL90, function(x) removeZero(x$inventory[,"Orders"])))[20:83],
     type = "l", ylim = c(290, 310))
plot(rowMeans(sapply(InvSim_ETSmae_CSL90, function(x) removeZero(x$inventory[,"Orders"])))[20:83],
     type = "l", ylim = c(290, 310))
plot(rowMeans(sapply(InvSim_ETSshr_CSL90, function(x) removeZero(x$inventory[,"Orders"])))[20:83],
     type = "l", ylim = c(290, 310))
plot(rowMeans(sapply(InvSim_ARIMA_CSL90, function(x) removeZero(x$inventory[,"Orders"])))[20:83],
     type = "l", ylim = c(290, 310))

matrix(c(rowMeans(sapply(InvSim_ETSmse_CSL90, function(x) x$accuracy)),
         rowMeans(sapply(InvSim_ETSmae_CSL90, function(x) x$accuracy)),
         rowMeans(sapply(InvSim_ETSshr_CSL90, function(x) x$accuracy)),
         rowMeans(sapply(InvSim_ARIMA_CSL90, function(x) x$accuracy))),
       ncol = 4, byrow = TRUE,
       dimnames = list(c("ETS-MSE", "ETS-MAE", "ETS-SHR", "ARIMA"), c("AME", "RMSE", "PinballUp", "|sCE|")))


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
