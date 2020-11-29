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

setZDiscretization = function(n,rho){
  s$T=generateTransitionMatrix(n,rho)
  s$Z=calculateZLevels(n)
  s$n=n
  s$rho=rho
}

setZDiscretization(10,0.7)

expectedHarvest = function(currentLvl){
  return (sum(s$T[currentLvl,]*s$Z))
}

varianceExpectedHarvest = function(currentLvl){
  mean = expectedHarvest(currentLvl)
  return (sum(s$T[currentLvl,]*(s$Z)^2)-mean^2)
}

# Figure 1
setZDiscretization(10,0.7)
e = c()
for (i in 1:10){
  e[i]=expectedHarvest(i)
}

plot(-2:2, -2:2*0.7, type="l")
points(s$Z[1:10],e, type="l")

# Figure 2 - Adjustment of the actual conditional variance by 1/(1-rho^2) necessary to replicate the chart
plot(c(-2,2),c(0,1.5))
for (rho in c(0.1,0.3,0.5,0.7,0.9)){
  setZDiscretization(10,rho)
  v = c()
  for (i in 1:10){
    v[i]=varianceExpectedHarvest(i)/(1-rho^2)
  }
  points(s$Z[1:10],v, type="l")
}

# Figure 3
setZDiscretization(10,0.7)
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
  s$P = switch(type, "linear" = function(a,b,x){return (a+b*x)})
  s$D = switch(type, "linear" = function(a,b,p){return ((p-a)/b)})
}

setDemandFunction("linear")

# set static parameters of the QMLE

q = new.env()     # q is for parameters of the QMLE so that the QMLE can be implemented agnostically
q$constants["r"] = 0.05

getR = function(){
  return (q$constants["r"])
}

# Transform to and from generic parameter vector theta using formula (43) (not the (43) that's actually (42))

toTheta = function(a,b,delta){
  return (c(a, log(-b), log(delta+0.05)))
}

fromTheta = function(theta){
  return (data.frame(a=theta[1],b=-exp(theta[2]),delta=-0.05+exp(theta[3])))
}

q$theta = toTheta(0.2, -0.15, 0.12)

getA = function(){
  return (fromTheta(q$theta)[[1]])
}

getB = function(){
  return (fromTheta(q$theta)[[2]])
}

getDelta = function(){
  return (fromTheta(q$theta)[[3]])
}

# Discretize X

setXDiscretization = function (m) {
  # minimum is 0 stored + minimum harvest, maximum is maximum harvest/delta with a safety factor allowing for delta to shrink
  safety_factor = 2
  s$X = (0:(m-1))/m*(s$Z[s$n]/getDelta()*safety_factor-s$Z[1])+s$Z[1]
  s$m=m
}

setXDiscretization(20)

# Iterate to improve price function f
# f(x,z) is stored as a matrix with f[k,i] representing f(X[k],Z[i])

initF = function(){
  s$f = matrix(nrow=s$m, ncol=s$n)
  a = getA()
  b = getB()
  for (k in 1:s$m){
    s$f[k,]=max(s$P(a,b,s$X[k]), 0)
  }
}

initF()


# Using formula (28)

iterateF = function (){
  a = getA()
  b = getB()
  delta = getDelta()
  newF = matrix(nrow=s$m, ncol=s$n)
  beta = (1-delta)/(1+getR())
  splines=list()
  for (i in 1:s$n){
    splines=append(splines,splinefun(s$X, s$f[,i]))
  }
  futureValue = function(x,i){
    v = 0
    y = s$D(a,b,splines[[i]](x))
    for (j in 1:s$n){
      v = v + s$T[i,j]*splines[[j]](s$Z[j]+(1-delta)*(x-y))
    }
    return (beta*v)
  }
  for (k in 1:s$m){
    currentPrice=s$P(a,b,s$X[k])
    for (i in 1:s$n){
      newF[k,i]=max(futureValue(s$X[k],i),currentPrice)
    }
  }
  s$f=newF
}


for(i in 1:50){
  iterateF()
}


plot(s$X,s$f[,1],type="l")
for(i in 2:s$n){
  points(s$X,s$f[,i],type="l")
}
