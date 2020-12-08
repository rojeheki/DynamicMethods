# Figure 1
s = com$cocoa
initEverything(s)
calculateEverything(s)

eH = c()
for (i in 1:10){
  eH[i]=expectedHarvest(i, s)
}

plot(-2:2, -2:2*0.7, type="l")
points(s$Z[1:10],eH, type="l")

# Figure 2 - Adjustment of the actual conditional variance by 1/(1-rho^2) necessary to replicate the chart
plot(c(-2,2),c(0,1.5))
for (rho in c(0.1,0.3,0.5,0.7,0.9)){
  setZDiscretization(10,rho,s)
  v = c()
  for (i in 1:10){
    v[i]=varianceExpectedHarvest(i,s)/(1-rho^2)
  }
  points(s$Z[1:10],v, type="l")
}

# Figure 3
setZDiscretization(10,0.7,s)
lvls = c(5)
r = runif(200)
for (i in 1:200) {
  j=1
  while (r[i] > s$Tz[lvls[i],j]) {
    r[i] = r[i]-s$Tz[lvls[i],j]
    j = j+1
  }
  lvls[i+1] = j
}
plot(1:201, s$Z[lvls], type="l")


# Figure 4
s = com$cocoa
initEverything(s)
calculateEverything(s)

plot(s$X,s$f[,1],type="l")
for(i in 2:s$nY){
  points(s$X,s$f[,i],type="l")
}


# Figure 5
idVis = data.frame()
for (i in 1:s$nX){
  for (j in 1:s$nY){
    nl = data.frame(x=s$X[i], z=s$Z[j], p=as.numeric(s$id[(i-1)*s$nY+j]))
    idVis = rbind(idVis,nl)
  }
}

idPlot = ggplot(idVis, aes(x=z,y=x,z=p)) + geom_contour_filled()
idPlot

# Figure 6
s = com$rice
initEverything(s, rho=0)
calculateEverything(s)

fi = splinefun(s$f[,1], s$X)
curve(fi(x), from=0, to=max(s$pDat))

# Figure 7
fm = splinefun(s$pDat, s$m)
curve(fm(x), from=0, to=max(s$pDat))

# Figure 8
fs = splinefun(s$pDat, s$s)
curve(fs(x), from=0, to=max(s$pDat))
