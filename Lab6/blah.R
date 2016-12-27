n=1000
x=rep(0,n)
y=rep(0,n)
z=rep(0,n)
xdist<- function(y, z, i){
  rnorm(1, mean= 0.899*y[i]+0.1*z[i], sd = 0.1899^0.5)
}
ydist <- function(x, z, i){
  rnorm(1, mean= 0.899*x[i]+0.1*z[i], sd = 0.1899^0.5)
}
zdist<- function(x, y, i){
  rnorm(1, mean= 0.053*x[i]+0.053*y[i], sd=0.989^0.5)
  
}
for (i in 2:n) {
  x[i]=xdist(y, z, i-1)
  y[i]=ydist(x, z, i)
  z[i]=zdist(x, y, i)
}
plot(x, type="l", main= "traceplot for x")
plot(x, type="l", main= "traceplot for y")