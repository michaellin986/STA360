## Initialize parameters ##
c = 0.25
n = 10^3
v = runif(1,-c/2,c/2)

## Initialize vectors ##
uval = rep(NA,n)
vval = rep(NA,n)

vval[1] = v
unif.min = max(0,abs(v))
unif.max = min(1,1-abs(v))

u = runif(1,unif.min,unif.max)

uval[1] = u;

## Gibbs Sampler ##
for(i in 2:n){
  v = runif(1,-c/2,c/2)
  vval[i] = v
  
  unif.min = max(0,abs(vval[i]))
  unif.max = min(1,1-abs(vval[i]))
  
  u = runif(1, unif.min, unif.max)
  
  uval[i] = u;
}

xval = uval + vval
yval = uval - vval

## Plot ##
png("scat_025x.png")
plot(xval, yval, pch=20, main = "Scatterplot (c=0.25)", xlab = "x", ylab = "y")
dev.off()

png("trace_025x.png")
plot(1:n,xval,pch=20, main = "Traceplot (c=0.25)", xlab = "n", ylab = "x")
dev.off()