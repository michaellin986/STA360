N = 1000
x = rep(0,N)
y = rep(0,N)
z = rep(0,N)

x[1] = rnorm(1, mean = 0, sd = 0.18990^0.5);
y[1] = rnorm(1, mean = 0, sd = 0.18990^0.5);
z[1] = rnorm(1, mean = 0.052632*x[1]+0.052632*y[1], sd = 0.98947^0.5)


for(i in 2:N){
  x[i] = rnorm(1, mean = 0.89899*y[i-1]+0.01010*z[i-1], sd = 0.18990^0.5)
  y[i] = rnorm(1, mean = 0.89899*x[i]+0.01010*z[i-1], sd =0.18990^0.5)
  z[i] = rnorm(1, mean = 0.052632*x[i]+0.052632*y[i], sd = 0.98947^0.5)
}

png("tracenorm.png")
plot(x, pch=20, xlab="N", ylab="X", main="Traceplot for X (Normal Gibbs)")
dev.off()

png("acfnorm.png")
acf(x)
dev.off()