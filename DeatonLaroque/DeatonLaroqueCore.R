# Based on Angus Deaton and Guy Laroque. Estimating a nonlinear rational expectations commodity price model with unobservable state variables
# Journal of Applied Econometrics, Vol. 10, S9-S40 (1995)
# https://doi.org/10.1002/jae.3950100503

library(cubature)
library(ggplot2)

# Import and organize data

dat = read.csv2("TimeSeries.csv", sep=",")
com = list(cocoa=new.env(), coffee=new.env(), copper=new.env(), rice=new.env(), sugar=new.env(), tin=new.env())

com$cocoa$pDat = c(dat$P_Cocoa)/c(dat$US_CPI)
com$cocoa$years = c(dat$Year)

com$coffee$pDat = c(dat$P_Coffee)/c(dat$US_CPI)
com$coffee$years = c(dat$Year)

com$copper$pDat = c(dat$P_Copper)/c(dat$US_CPI)
com$copper$years = c(dat$Year)

com$rice$pDat = c(dat$P_Rice)/c(dat$US_CPI)
com$rice$years = c(dat$Year)

com$sugar$pDat = c(dat$P_Sugar)/c(dat$US_CPI)
com$sugar$years = c(dat$Year)

com$tin$pDat = c(dat$P_Tin)/c(dat$US_CPI)
com$tin$years = c(dat$Year)

calculateT = function(e){
  e$t = length(e$years)
}

for (e in com) {
  calculateT(e)
}


# 3.1 Discretization of the Harvest

generateZTransitionMatrix = function(e) {
  f = function(x) {
    exp(-(x[1]^2-2*e$rho*x[1]*x[2]+x[2]^2)/(2*(1-e$rho^2)))/(2*pi*sqrt(1-e$rho^2))
  }
  theta = qnorm(0:e$nY/e$nY)
  e$Tz = matrix(nrow = e$nY, ncol = e$nY)
  for(i in 1:e$nY) {
    for(j in 1:i){
      e$Tz[i,j]=e$Tz[j,i]=adaptIntegrate(f, lowerLimit=c(theta[i],theta[j]),upperLimit=c(theta[i+1],theta[j+1]))$integral*e$nY
    }
  }
}

calculateZLevels = function(e){
  theta = qnorm(0:e$nY/e$nY)
  e$Z = (-(dnorm(theta[-1])-dnorm(theta[-(e$nY+1)]))*e$nY)
}

setZDiscretization = function(n,rho,e){
  e$nY=n
  e$rho=rho
  generateZTransitionMatrix(e)
  calculateZLevels(e)
}

expectedHarvest = function(currentZLvl, e){
  return (sum(e$Tz[currentZLvl,]*e$Z))
}

varianceExpectedHarvest = function(currentZLvl, e){
  mean = expectedHarvest(currentZLvl, e)
  return (sum(e$Tz[currentZLvl,]*(e$Z)^2)-mean^2)
}

# Set (Inverse) Demand function and real interest

setLinearPriceFunction = function(e) {
  e$P = function(x){return (e$a+e$b*x)}
  e$D = function(p){return ((p-e$a)/e$b)}
}

setRealInterest = function(r, e) {
  e$r = 0.05
}


# Transform to and from generic parameter vector theta using formula (43) (not the (43) that's actually (42))

toTheta = function(a,b,delta,e){
  e$a = a
  e$b = b
  e$delta = delta
  e$theta = c(a, log(-b), log(delta+0.05))
}

fromTheta = function(theta, e){
  e$theta = theta
  e$a = theta[1]
  e$b = -exp(theta[2])
  e$delta = -0.05+exp(theta[3])
}

initTheta = function (e) {
  toTheta(3, -2, 0.12, e)
}


# 3.2 Discretization of the availability

setXDiscretization = function (m,e) {
  # minimum is 0 stored + minimum harvest, maximum is maximum harvest/delta with a safety factor allowing for delta to shrink
  safety_factor = 2
  e$X = (0:(m-1))/m*(e$Z[e$nY]/e$delta*safety_factor-e$Z[1])+e$Z[1]
  e$nX=m
}


# 3.2 Iterate to improve price function f
# f(x,z) is stored as a matrix with f[k,i] representing f(X[k],Z[i])

initF = function(e){
  e$f = matrix(nrow=e$nX, ncol=e$nY)
  for (k in 1:e$nX){
    e$f[k,]=max(e$P(e$X[k]), 0)
  }
}


# Using formula (28)

iterateFOnce = function (e){
  newF = matrix(nrow=e$nX, ncol=e$nY)
  beta = (1-e$delta)/(1+e$r)
  splines=list()
  for (i in 1:e$nY){
    splines=append(splines,splinefun(e$X, e$f[,i]))
  }
  futureValue = function(x,i){
    v = 0
    y = e$D(splines[[i]](x))
    for (j in 1:e$nY){
      v = v + e$Tz[i,j]*splines[[j]](e$Z[j]+(1-e$delta)*(x-y))
    }
    return (beta*v)
  }
  for (k in 1:e$nX){
    currentPrice=e$P(e$X[k])
    for (i in 1:e$nY){
      newF[k,i]=max(futureValue(e$X[k],i),currentPrice)
    }
  }
  e$f=newF
}

iterateF = function(e, dCutoff = 1e-5, nCutoff = 100) {
  initF(e)
  lastF = e$f+Inf
  for (i in 1:nCutoff) {
    deviation = max(abs(e$f - lastF))
    if(deviation<dCutoff) {
      # print(i)
      return (TRUE)
    }
    lastF = e$f
    iterateFOnce(e)
  }
  return (FALSE)
}


# 3.3 Generate the combined storage/harvest transition matrix

# calculates g(x,z) from formulas (36), (39)
expectedStorage = function(currentXLvl, currentZLvl, e) {
  return ((1-e$delta)*(e$X[currentXLvl]-e$D(e$f[currentXLvl,currentZLvl]))+e$rho*e$Z[currentZLvl])
}


# XZ Transition matrix will be ordered by Z on a large scale, X on a small scale
# X1Z1, X2Z1, ..., XmZ1, X1Z2, X2Z2, ..., xmZn

generateXZTransitionMatrix = function(e){
  theta = qnorm(0:e$nY/e$nY)
  critL = matrix(nrow = e$nX*e$nY, ncol = e$nX*e$nY)
  critR = matrix(nrow = e$nX*e$nY, ncol = e$nX*e$nY)
  xd = append(append(Inf,(e$X[2:e$nX]-e$X[1:(e$nX-1)])/2),Inf) # delta/2 between X, xd[i] is left of X[i], xd[i+1] right
  for (i in 1:e$nX) {
    for (j in 1:e$nY) {
      c_row = (i-1)*e$nY + j
      for (k in 1:e$nX) {
        for (l in 1:e$nY) {
          c_column = (k-1)*e$nY + l
          xdc = e$X[i] - expectedStorage(k,l,e)
          critXL = xdc - xd[i]
          critXR = xdc + xd[i+1]
          critZL = theta[j]-e$rho*e$Z[l]
          critZR = theta[j+1]-e$rho*e$Z[l]
          critL[c_row, c_column] = max(critXL,critZL)
          critR[c_row, c_column] = min(critXR,critZR)
        }
      }
    }
  }
  e$Txz = matrix(pmax(0, pnorm(critR)-pnorm(critL)),nrow = e$nX*e$nY, ncol = e$nX*e$nY)
}

# Calculate invariant distribution

calculateInvariantDistribution = function(e) {
  e$id = eigen(e$Txz)$vectors[,1]
  if(Re(sum(e$id)) < 0) {
    e$id = - e$id
  }
}

# 5.1 Estimate harvest level probabilities from price series

calculateGamma = function(e) {
  e$gamma = matrix(nrow=e$nY, ncol=e$t)
  helper = matrix(nrow=e$nY, ncol=e$t)    # contains non-normalized probabilities
  for (i in 1:e$nY) {
    inversePriceFunction = splinefun(e$f[,i], e$X)
    probabilitySpline = splinefun(e$X, e$id[((i-1)*e$nX+1):(i*e$nX)])
    for (t in 1:e$t) {
      assumedX = inversePriceFunction(e$pDat[t])
      helper[i,t] = max(0, probabilitySpline(assumedX))
      # print(paste(i,t,assumedX,probabilitySpline(assumedX)))
    }
  }
  for (i in 1:e$nY) {
    for (t in 1:e$t) {
      e$gamma[i,t] = helper[i,t]/sum(helper[,t])
    }
  }
  for (t in e$t){
    stopifnot(round(sum(e$gamma[,t]), digits=6) == 1)
  }
}


# 5.1 Calculate 1-period-ahead expectations and variances

# Using formula (45), calculate a matrix condM that has m(pt, Zi) in condM[t,i]
calculateCondM = function(e) {
  e$condM = matrix(nrow=e$t, ncol=e$nY)
  forwardSplines=list()
  inverseSplines=list()
  for (i in 1:e$nY){
    forwardSplines=append(forwardSplines,splinefun(e$X, e$f[,i]))
    inverseSplines=append(inverseSplines,splinefun(e$f[,i], e$X))
  }
  for (t in 1:e$t){
    for (i in 1:e$nY){
      value = 0
      for (j in 1:e$nY){
        value=value+ e$Tz[j,i]*forwardSplines[[j]]((1-e$delta)*(inverseSplines[[i]](e$pDat[t])-e$D(e$pDat[t]))+e$Z[j])
      }
      e$condM[t,i] = value
    }
  }
}


# Using formula (48), calculate a vector m that has m(pt) in m[t]
calculateM = function (e) {
  e$m = diag(e$condM %*% e$gamma)
}


# Using formula (46), calculate a matrix condS that has s(pt, Zi) in condS[t,i]
calculateCondS = function (e) {
  e$condS = matrix(nrow=e$t, ncol=e$nY)
  forwardSplines=list()
  inverseSplines=list()
  for (i in 1:e$nY){
    forwardSplines=append(forwardSplines,splinefun(e$X, e$f[,i]))
    inverseSplines=append(inverseSplines,splinefun(e$f[,i], e$X))
  }
  for (t in 1:e$t){
    for (i in 1:e$nY){
      value = 0
      for (j in 1:e$nY){
        value=value+ e$Tz[j,i]*(forwardSplines[[j]]((1-e$delta)*(inverseSplines[[i]](e$pDat[t])-e$D(e$pDat[t]))+e$Z[j]))^2
      }
      e$condS[t,i] = value - (e$condM[t,i])^2
    }
  }
}


# Using formula (49), calculate a vector s that has s(pt) in s[t]
calculateS = function (e) {
  e$s = diag(e$condS %*% e$gamma) + diag((e$condM-e$m)^2 %*% e$gamma)
}



# Calculate (Log) Pseudo Likelihood Function

calculatePLF = function (e) {
  e$PLF = -0.5*((e$t-1)*log(2*pi) + sum(log(e$s[-e$t])) + sum((e$pDat[-1]-e$m[-e$t])^2/e$s[-e$t]))
}


# Complete model initialization

initEverything = function (e, r = 0.05, a = 3, b = -2, delta = 0.12, rho = 0.7, nY = 10, nX = 20) {
  setLinearPriceFunction(e)
  setRealInterest(r, e)
  toTheta(a, b, delta, e)
  setZDiscretization(nY,rho,e)
  setXDiscretization(nX,e)
}

# Calculate the entire model from given nY, nX, theta, rho, pDat, t

calculateEverything = function (e, rho, theta) {
  setZDiscretization(e$nY,rho,e)
  setXDiscretization(e$nX,e)
  fromTheta(theta,e)
  iterateF(e)
  generateXZTransitionMatrix(e)
  calculateInvariantDistribution(e)
  calculateGamma(e)
  calculateCondM(e)
  calculateM(e)
  calculateCondS(e)
  calculateS(e)
  calculatePLF(e)
  print(c(round(e$a, digits=2), round(e$b, digits=2), round(e$delta, digits=2), round(e$rho, digits=2), round(e$PLF)))
  return (e$PLF)
}

