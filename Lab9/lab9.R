# load data
x = read.table("data.txt", header = F)[,1]
library(MCMCpack)

## Define prior parameters ##
a = 0; b = 0;
mu.0 = 0; sig2.0 = 10;
n = length(x)

## Define prior ##
mu.temp = rnorm(1, mean = (mu.0*sig2.0 + sig2.0*sum(log(x)))/(sig2.0 + n*sig2.0),
                sd = sqrt((sig2.0*sig2.0)/(sig2.0 + n*sig2.0)))
sig2.temp = rinvgamma(1, shape = a + n/2, 
                      scale = b + 1/2*sum( (log(x) - mu.temp)^2 ))

N = 10000
MU = rep(0, N)
SIG2 = rep(0, N)

## Gibbs Sampler ##
for(i in 1:N){
  mu.temp = rnorm(1, mean = (mu.0*sig2.temp + sig2.0*sum(log(x)))/(sig2.temp + n*sig2.0),
             sd = sqrt((sig2.0*sig2.temp)/(sig2.temp + n*sig2.0)))
  sig2.temp = rinvgamma(1, shape = a + n/2, 
                        scale = b + 1/2*sum( (log(x) - mu.temp)^2 ))
  MU[i] = mu.temp
  SIG2[i] = sig2.temp  
}

## Running Avg ##
MU.sum = rep(0, N)
MU.avg = rep(0, N)
MU.sum[1] = MU[1]
for(i in 2:N){
  MU.sum[i] = MU[i] + MU.sum[i-1]
  MU.avg[i] = MU.sum[i]/i
}

SIG2.sum = rep(0, N)
SIG2.avg = rep(0, N)
SIG2.sum[1] = SIG2[1]
for(i in 2:N){
  SIG2.sum[i] = SIG2[i] + SIG2.sum[i-1]
  SIG2.avg[i] = SIG2.sum[i]/i
}

## Plot ##
png("MUplot.png")
plot(MU, pch = '.', xlab = "N", ylab = expression(mu), 
     main = "Traceplot and Running Average")
lines(MU.avg, lwd = 2)
dev.off()

png("SIGplot.png")
plot(SIG2, pch = '.', xlab = "N", ylab = expression(sigma^{2}), 
     main = "Traceplot and Running Average")
lines(SIG2.avg, lwd = 2)
dev.off()

## Confidence Interval ##
mean.vec = exp(MU + SIG2/2)
var.vec = (exp(SIG2)-1)*exp(2*MU + SIG2)

mean.ci = quantile(mean.vec, c(0.025, 0.975))
var.ci = quantile(var.vec, c(0.025, 0.975))