library(MASS)
N = 1000

xy.samp = matrix(rep(0,2*N), nrow=N, ncol=2, byrow=T)
colnames(xy.samp) = c("x","y")

z = rep(0,N)
z[1] = rnorm(1, mean = 0, sd = 0.98^0.5)

for(i in 2:N){
  xy.mean = c(0.1*z[i-1],0.1*z[i-1])
  xy.var = matrix(c(0.99,0.89,0.89,0.99), nrow=2, ncol=2, byrow = T)
  
  xy.samp[i,] = mvrnorm(1, xy.mean, xy.var)
  z[i] = rnorm(1, mean = 0.052632*xy.samp[i,1]+0.052632*xy.samp[i,2], sd = 0.98947^0.5)
}
x = xy.samp[,1]
y = xy.samp[,2]

png("traceblock.png")
plot(x, pch=20, xlab="N", ylab="X", main="Traceplot for X (Block Update)")
dev.off()

png("acfblock.png")
acf(x)
dev.off()