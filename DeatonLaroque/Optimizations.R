funIID = function (data, par) {
  return (- calculateEverything(data[[1]], data[[2]], par))
}

funAR = function (data, par) {
  return (- calculateEverything(data[[1]], par[4], par[1:3]))
}

initEverything(com$cocoa)
initEverything(com$coffee)

system.time(calculateEverything(com$cocoa,0,c(3,0.693,-1.772)))

optIID = optim(par = c(3, 0.693, -1.772), fn = funIID, data = c(com$coffee, 0), control=list(maxit=100))

optAR = optim(par = c(3, 0.693, -1.772, 0.7), fn = funIID, data = c(com$cocoa))

