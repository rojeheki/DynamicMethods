funIID = function (data, par) {
  return (- calculateEverything(data[[1]], data[[2]], par))
}

funAR = function (data, par) {
  return (- calculateEverything(data[[1]], 1-exp(par[4]), par[1:3]))
}


initEverything(com$coffee)

system.time(calculateEverything(com$coffee,0,c(3,0.693,-1.772)))

optIID = optim(par = c(2, 1, 0.6), fn = funIID, data = c(com$coffee, 0), control=list(maxit=200))

optAR = optim(par = c(3, 0.693, -1.772, -1.2), fn = funAR, data = c(com$coffee), control=list(maxit=200))

