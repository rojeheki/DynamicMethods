# Setup sample data and benchmark time requirement

initEverything(com$coffee)

system.time(calculateEverything(com$coffee,3,-2,0.12,0))

# Using transformations so that the optimization can be unbounded

funIID = function (data, par) {
  a = par[1]
  b = -0.01-exp(par[2])
  delta = -0.05+exp(par[3])
  return (- calculateEverything(data[[1]], a, b, delta, 0))
}

funAR = function (data, par) {
  a = par[1]
  b = -0.01-exp(par[2])
  delta = -0.05+exp(par[3])
  rho = 0.99-exp(par[4])
  return (- calculateEverything(data[[1]], a, b, delta, rho))
}

convertFromPar = function(par){
  a = par[1]
  b = -0.01-exp(par[2])
  delta = -0.05+exp(par[3])
  rho = 0.99-exp(par[4])
  return (c(a,b,delta,rho))
}

convertToPar = function(a,b,delta,rho){
  par = c(a, log(-b+0.01), log(delta+0.05), log(0.99-rho))
  return (par)
}


# Nelder-Mead

optIID = optim(par = convertToPar(15, -8, 0.12, 0), fn = funIID, data = c(com$coffee), control=list(maxit=200, alpha=1.2, gamma=2))

optAR = optim(par = convertToPar(15, -8, 0.12, 0.7), fn = funAR, data = c(com$coffee), control=list(maxit=200))


# BFGS (tries to set parameters to Inf or -Inf)

optIID = optim(par = convertToPar(15, -8, 0.12, 0), fn = funIID, data = c(com$coffee), method = "BFGS", control=list(maxit=200))




# Using bounded optimization

initEverything(com$coffee)

funIID = function (data, par) {
  return (- calculateEverything(data[[1]], par[1], par[2], par[3], 0))
}

funAR = function (data, par) {
  return (- calculateEverything(data[[1]], par[1], par[2], par[3], par[4]))
}

optIID = optim(par = c(15, -1, 0.12), lower = c(0.1,-50,0), upper = c(100,0,1), fn = funIID, data = c(com$coffee), method="L-BFGS-B", control=list(maxit=200))

