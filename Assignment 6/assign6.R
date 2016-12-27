## Initialize parameters ##
c = 0.02
n = 10^3
y = runif(1,0,1)

## Initialize vectors ##
xval = rep(NA,n)
yval = rep(NA,n)

yval[1] = y

x = runif(1,0,1)
while(abs(x-y)>=c){
  x = runif(1,0,1)
}
xval[1] = x;

## Gibbs Sampler ##
for(i in 2:n){
  y = runif(1,0,1)
  while(abs(xval[i-1]-y)>=c){
    y = runif(1,0,1)
  }
  yval[i] = y
  
  x = runif(1,0,1)
  while(abs(yval[i]-x)>=c){
    x = runif(1,0,1)
  }
  xval[i] = x
}

## Plot ##
##png("scat_002.png")
plot(xval, yval, pch=20, main = "Scatterplot (c=0.02)", xlab = "x", ylab = "y")
##dev.off()

##png("trace_002.png")
plot(1:n,xval,pch=20, main = "Traceplot (c=0.02)", xlab = "n", ylab = "x")
##dev.off()