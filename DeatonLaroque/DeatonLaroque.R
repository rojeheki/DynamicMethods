library(cubature)

# 3.1

generateTransitionMatrix = function(n,rho) {
  f = function(x) {
    exp(-(x[1]^2-2*rho*x[1]*x[2]+x[2]^2)/(2*(1-rho^2)))/(2*pi*sqrt(1-rho^2))
  }
  theta = qnorm(0:n/n)
  T = matrix(nrow = n, ncol = n)
  for(i in 1:n) {
    for(j in 1:i){
      T[i,j]=T[j,i]=adaptIntegrate(f, lowerLimit=c(theta[i],theta[j]),upperLimit=c(theta[i+1],theta[j+1]))$integral*n
    }
  }
  return (T)
}

calculateZLevels = function(n){
  theta = qnorm(0:n/n)
  return (-(dnorm(theta[-1])-dnorm(theta[-(n+1)]))*n)
}

s = new.env()

setDiscretization = function(n,rho){
  s$T=generateTransitionMatrix(n,rho)
  s$Z=calculateZLevels(n)
  s$n=n
  s$rho=rho
}

setDiscretization(10,0.7)

expectedHarvest = function(currentLvl){
  return (sum(s$T[currentLvl,]*s$Z))
}

varianceExpectedHarvest = function(currentLvl){
  mean = expectedHarvest(currentLvl)
  return (sum(s$T[currentLvl,]*(s$Z)^2)-mean^2)
}

# Figure 1
setDiscretization(10,0.7)
e = c()
for (i in 1:10){
  e[i]=expectedHarvest(i)
}

plot(-2:2, -2:2*0.7, type="l")
points(s$Z[1:10],e, type="l")

# Figure 2 - Adjustment of the actual conditional variance by 1/(1-rho^2) necessary to replicate the chart
plot(c(-2,2),c(0,1.5))
for (rho in c(0.1,0.3,0.5,0.7,0.9)){
  setDiscretization(10,rho)
  v = c()
  for (i in 1:10){
    v[i]=varianceExpectedHarvest(i)/(1-rho^2)
  }
  points(s$Z[1:10],v, type="l")
}

# Figure 3
setDiscretization(10,0.7)
lvls = c(5)
r = runif(200)
for (i in 1:200) {
  j=1
  while (r[i] > s$T[lvls[i],j]) {
    r[i] = r[i]-s$T[lvls[i],j]
    j = j+1
  }
  lvls[i+1] = j
}
plot(1:201, s$Z[lvls], type="l")

# (Inverse) Demand function

setDemandFunction = function(type){
  s$P = switch(type, "linear" = function(a,b,x){return (a*x+b)})
  s$D = switch(type, "linear" = function(a,b,x){return ((p-b)/a)})
}

setDemandFunction("linear")
