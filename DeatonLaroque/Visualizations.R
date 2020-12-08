# Figure 1
setZDiscretization(10,0.7, com$cocoa)
eH = c()
for (i in 1:10){
  eH[i]=expectedHarvest(i, com$cocoa)
}

plot(-2:2, -2:2*0.7, type="l")
points(com$cocoa$Z[1:10],eH, type="l")

# Figure 2 - Adjustment of the actual conditional variance by 1/(1-rho^2) necessary to replicate the chart
plot(c(-2,2),c(0,1.5))
for (rho in c(0.1,0.3,0.5,0.7,0.9)){
  setZDiscretization(10,rho,com$cocoa)
  v = c()
  for (i in 1:10){
    v[i]=varianceExpectedHarvest(i,com$cocoa)/(1-rho^2)
  }
  points(com$cocoa$Z[1:10],v, type="l")
}

# Figure 3
setZDiscretization(10,0.7,com$cocoa)
lvls = c(5)
r = runif(200)
for (i in 1:200) {
  j=1
  while (r[i] > com$cocoa$Tz[lvls[i],j]) {
    r[i] = r[i]-com$cocoa$Tz[lvls[i],j]
    j = j+1
  }
  lvls[i+1] = j
}
plot(1:201, com$cocoa$Z[lvls], type="l")


# Figure 4
plot(com$cocoa$X,com$cocoa$f[,1],type="l")
for(i in 2:com$cocoa$nY){
  points(com$cocoa$X,com$cocoa$f[,i],type="l")
}


# Visualize invariant distribution

s = com$cocoa

initEverything(s)
calculateEverything(s,s$rho,s$theta)

idVis = data.frame()
for (i in 1:s$nX){
  for (j in 1:s$nY){
    nl = data.frame(x=s$X[i], z=s$Z[j], p=as.numeric(s$id[(i-1)*s$nY+j]))
    idVis = rbind(idVis,nl)
  }
}


idPlot = ggplot(idVis, aes(x=z,y=x,z=p)) + geom_contour_filled()
idPlot

# Generate Monte-Carlo data

s$mc = data.frame(z = 1, x = 1, p=s$f[5,10])
r = runif(s$t)
for (i in 1:(s$t-1)) {
  j=1
  while (r[i] > s$Txz[j,(s$mc[i,1]-1)*s$nY+s$mc[i,2]]) {
    r[i] = r[i]-s$Txz[j,(s$mc[i,1]-1)*s$nY+s$mc[i,2]]
    j = j+1
  }
  s$mc[i+1,1] = (j-1)%%s$nY+1
  s$mc[i+1,2] = floor((j-1)/s$nY)+1
  s$mc[i+1,3] = s$f[s$mc[i+1,1],s$mc[i+1,2]]
}

plot(1:s$t,s$mc$p, type="l")