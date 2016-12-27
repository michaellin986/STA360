## load data ##
data = read.table("data.txt", header = T)
library(MCMCpack)
data$weekend = rep(0,nrow(data))
data$weekend[data$day=='Saturday' | data$day=='Sunday'] = 1
X = data[,-2]
Y1 = subset(X, weekend==0)
row.names(Y1) = NULL
Y1 = Y1[,1]
Y2 = subset(X, weekend==1)
row.names(Y2) = NULL
Y2 = Y2[,1]

## set priors ##
m1 = m2 = 0
s21 = s22 = 1
a1 = a2 = 1
b1 = b2 = 1
n1 = length(Y1)
n2 = length(Y2)
B = 2000
N = 12000

## set seed ##
sig1.temp = 1 #for delta variance
sig2.temp = 1

## sample priors ##
mu2.temp = rnorm(1, mean = (s22*sum(log(Y2))+sig2.temp*m2)/(s22*n2 + sig2.temp),
                 sd = sqrt((sig2.temp*s22)/(s22*n2+sig2.temp)))
del.temp = -1
while(del.temp < 0){
del.temp = rnorm(1, mean = (s21*sum(log(Y1))+m1*sig1.temp-n1*mu2.temp*s21)/(n1*s21+sig1.temp),
                 sd = sqrt(sig1.temp*s21/(n1*s21+sig1.temp)) )
}
mu1.temp = mu2.temp + del.temp


DEL = rep(0, N)
MU1 = rep(0, N)
MU2 = rep(0, N)
SIG1 = rep(0, N)
SIG2 = rep(0, N)


for(i in 1:N){
  sig1.temp = rinvgamma(1, shape = a1+n1/2, scale = b1 + 0.5*sum((log(Y1)-mu1.temp)^2))
  
  sig2.temp = rinvgamma(1, shape = a2+n2/2, scale = b2 + 0.5*sum((log(Y2)-mu2.temp)^2))
  
  mu2.temp = rnorm(1, mean = (s22*sum(log(Y2))+sig2.temp*m2)/(s22*n2 + sig2.temp),
                   sd = sqrt((sig2.temp*s22)/(s22*n2+sig2.temp)))
  
  del.temp = -1
  while(del.temp <= 0){
    del.temp = rnorm(1, mean = (s21*sum(log(Y1))+m1*sig1.temp-n1*mu2.temp*s21)/(n1*s21+sig1.temp),
                     sd = sqrt(sig1.temp*s21/(n1*s21+sig1.temp)) )
  }
  
  mu1.temp = mu2.temp + del.temp
  
  DEL[i] = del.temp
  MU1[i] = mu1.temp
  MU2[i] = mu2.temp
  SIG1[i] = sig1.temp
  SIG2[i] = sig2.temp
}


MU1 = MU1[(B+1):N]
MU2 = MU2[(B+1):N]
SIG1 = SIG1[(B+1):N]
SIG2 = SIG2[(B+1):N]

N = 10000
## Running Avg ##
MU1.sum = rep(0, N)
MU1.avg = rep(0, N)
MU1.avg[1] = MU1[1]
MU1.sum[1] = MU1[1]
for(i in 2:N){
  MU1.sum[i] = MU1[i] + MU1.sum[i-1]
  MU1.avg[i] = MU1.sum[i]/i
}

SIG1.sum = rep(0, N)
SIG1.avg = rep(0, N)
SIG1.avg[1] = SIG1[1]
SIG1.sum[1] = SIG1[1]
for(i in 2:N){
  SIG1.sum[i] = SIG1[i] + SIG1.sum[i-1]
  SIG1.avg[i] = SIG1.sum[i]/i
}


## Plot MU1 and SIG1 ##
png("MU1plot.png")
traceplot(mcmc(MU1), xlab = "N", ylab = expression(mu[1]), 
     main = "Traceplot")
dev.off()

png("SIG1plot.png")
traceplot(mcmc(SIG1), xlab = "N", ylab = expression(sigma[1]^{2}), 
     main = "Traceplot")
dev.off()

png("MU1MCMC.png")
plot(MU1, pch = '.', xlab = "N", ylab = expression(mu[1]), 
          main = "MCMC")
lines(MU1.avg, lwd = 2)
dev.off()

png("SIG1MCMC.png")
plot(SIG1, pch = '.', xlab = "N", ylab = expression(sigma[1]^{2}), 
          main = "MCMC")
lines(SIG1.avg, lwd = 2)
dev.off()

## Mean ##
MU1.mean = mean(MU1)
MU2.mean = mean(MU2)
SIG1.mean = mean(SIG1)
SIG2.mean = mean(SIG2)

## CI ##
MU1.ci = quantile(MU1, c(0.025, 0.975))
MU2.ci = quantile(MU2, c(0.025, 0.975))
SIG1.ci = quantile(SIG1, c(0.025, 0.975))
SIG2.ci = quantile(SIG2, c(0.025, 0.975))

## Posterior Probability ##
mean(MU1>MU2)
mean(SIG1>SIG2)

## Posterior Predictive ##
wk = we = rep(0, 10000)
for(i in 1:10000){
  wk[i] = rlnorm(1, meanlog = MU1[i], sdlog = sqrt(SIG1[i]))
  we[i] = rlnorm(1, meanlog = MU2[i], sdlog = sqrt(SIG2[i]))
}
mean(wk>we)