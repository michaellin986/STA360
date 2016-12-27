## load data ##
data = read.table("data1.txt", header = T)
library(MCMCpack)
library(truncnorm)
data$weekday = rbinom(100,1,0.5)

X = data
Y1 = subset(X, weekday==1)
row.names(Y1) = NULL
Y1 = Y1[,1]
Y2 = subset(X, weekday==0)
row.names(Y2) = NULL
Y2 = Y2[,1]
Y = X[,1]
zero.mat = matrix(0, nrow = length(Y), ncol = 4)
Y.dat = cbind(Y, zero.mat)
colnames(Y.dat) = c("Y","w","Y*w","Y*(1-w)", "pw=1")
#Y[,1] =Y | Y[,2] = w | Y[,3] = Yw | Y[,4] = Y(1-w) | Y[,5] = p(w=1)#

## set prior parameters ##
p.a = 5; p.b = 2
m1 = m2 = 0
s21 = s22 = 1
a1 = a2 = 1
b1 = b2 = 1
n = nrow(Y.dat)
B = 2000
N = 10000 + B

## set seed priors ##
mu1.temp = 0 #mean(Y1)
sig1.temp = 10 #var(Y1)
mu2.temp = 0 #mean(Y2)
sig2.temp = 10 #var(Y2)


## sample priors ##
# sample p #
wlist = list()
wvec = Y.dat[,2]
wvec.temp = rep(0, n) #new w's
pw1 = rep(0, n) #p(w=1)
pyw1 = rep(0, n) #p(y|w=1)
pyw0 = rep(0, n) #p(y|w=0)

p.temp = rbeta(1, p.a + sum(wvec), n + p.b - sum(wvec))


for(j in 1:n){
  pyw1[j] = dlnorm(Y.dat[j,1], meanlog = mu1.temp, sdlog = sqrt(sig1.temp))
  pyw0[j] = dlnorm(Y.dat[j,1], meanlog = mu2.temp, sdlog = sqrt(sig2.temp))
  pw1[j] = pyw1[j]*p.temp/(pyw1[j]*p.temp + pyw0[j]*(1-p.temp))
  wvec.temp[j] = rbinom(1, 1, pw1[j])
}



Y.dat[,5] = pw1
Y.dat[,2] = wvec.temp
Y.dat[,3] = Y.dat[,1]*Y.dat[,2]
Y.dat[,4] = Y.dat[,1]*(1-Y.dat[,2])

n1 = sum(Y.dat[,2])
n2 = n-n1

mu2.temp = rnorm(1, mean = (s22*sum(log(Y.dat[,4][Y.dat[,4]!=0]))+sig2.temp*m2)/(s22*n2 + sig2.temp),
                 sd = sqrt((sig2.temp*s22)/(s22*n2+sig2.temp)))

del.temp = rtruncnorm(1, a = 0, b = Inf, mean = (s21*sum(log(Y.dat[,3][Y.dat[,3]!=0]))+m1*sig1.temp-n1*mu2.temp*s21)/(n1*s21+sig1.temp),
                 sd = sqrt(sig1.temp*s21/(n1*s21+sig1.temp)) )

mu1.temp = mu2.temp + del.temp


DEL = rep(0, N)
MU1 = rep(0, N)
MU2 = rep(0, N)
SIG1 = rep(0, N)
SIG2 = rep(0, N)



## gibbs sampler ##

for(i in 1:N){
  wvec = Y.dat[,2]
  wvec.temp = rep(0, n) #new w's
  pw1 = rep(0, n) #p(w=1)
  pyw1 = rep(0, n) #p(y|w=1)
  pyw0 = rep(0, n) #p(y|w=0)
  
  
  p.temp = rbeta(1, p.a + sum(wvec), n + p.b - sum(wvec))
  
  
  for(j in 1:n){
    pyw1[j] = dlnorm(Y.dat[j,1], meanlog = mu1.temp, sdlog = sqrt(sig1.temp))
    pyw0[j] = dlnorm(Y.dat[j,1], meanlog = mu2.temp, sdlog = sqrt(sig2.temp))
    pw1[j] = pyw1[j]*p.temp/(pyw1[j]*p.temp + pyw0[j]*(1-p.temp))
    wvec.temp[j] = rbinom(1, 1, pw1[j])
  }
  
  wlist[[i]] = wvec.temp
  
  Y.dat[,5] = pw1
  Y.dat[,2] = wvec.temp
  Y.dat[,3] = Y.dat[,1]*Y.dat[,2]
  Y.dat[,4] = Y.dat[,1]*(1-Y.dat[,2])
  
  n1 = sum(Y.dat[,2]==1)
  n2 = n-n1  
  
  mu2.temp = rnorm(1, mean = (s22*sum(log(Y.dat[,4][Y.dat[,4]!=0]))+sig2.temp*m2)/(s22*n2 + sig2.temp),
                   sd = sqrt((sig2.temp*s22)/(s22*n2+sig2.temp)))
    
  sig2.temp = rinvgamma(1, shape = a2+n2/2, scale = b2 + 0.5*sum((log(Y.dat[,4][Y.dat[,4]!=0])-mu2.temp)^2))
  
  sig1.temp = rinvgamma(1, shape = a1+n1/2, scale = b1 + 0.5*sum((log(Y.dat[,3][Y.dat[,3]!=0])-mu1.temp)^2))
  
  
  
  
  del.temp = rtruncnorm(1, a = 0, b = Inf, mean = (s21*sum(log(Y.dat[,3][Y.dat[,3]!=0]))+m1*sig1.temp-n1*mu2.temp*s21)/(n1*s21+sig1.temp),
                   sd = sqrt(sig1.temp*s21/(n1*s21+sig1.temp)) )
  
  
  mu1.temp = mu2.temp + del.temp
  
  DEL[i] = del.temp
  MU1[i] = mu1.temp
  MU2[i] = mu2.temp
  SIG1[i] = sig1.temp
  SIG2[i] = sig2.temp
  
  if(i %% 1200 ==0){
    print(i/N)
  }
}


MU1 = MU1[(B+1):N]
MU2 = MU2[(B+1):N]
SIG1 = SIG1[(B+1):N]
SIG2 = SIG2[(B+1):N]
N = N-B

# traceplots
traceplot(mcmc(MU1))
traceplot(mcmc(MU2))
traceplot(mcmc(SIG1))
traceplot(mcmc(SIG2))

# running avg
MU1.avg = cumsum(MU1)/seq(1,N)
MU2.avg = cumsum(MU2)/seq(1,N)
SIG1.avg = cumsum(SIG1)/seq(1,N)
SIG2.avg = cumsum(SIG2)/seq(1,N)

# No. 2: MCMC plot
png('mu1.png')
plot(MU1, pch = '.', xlab = "N", ylab = expression(mu[1]), 
     main = "MCMC and Running Average")
lines(MU1.avg, lwd = 2)
dev.off()

png('mu2.png')
plot(MU2, pch = '.', xlab = "N", ylab = expression(mu[2]), 
     main = "MCMC and Running Average")
lines(MU2.avg, lwd = 2)
dev.off()

png('sig1.png')
plot(SIG1, pch = '.', xlab = "N", ylab = expression(sigma[1]^{2}), 
     main = "MCMC and Running Average")
lines(SIG1.avg, lwd = 2)
dev.off()

png('sig2.png')
plot(SIG2, pch = '.', xlab = "N", ylab = expression(sigma[2]^{2}), 
     main = "MCMC and Running Average")
lines(SIG2.avg, lwd = 2)
dev.off()

# No. 3
png('post.png')
plot(MU1, MU2, pch = 20, xlab = expression(mu[1]), ylab = expression(mu[2]),
     main = "Posterior mu_1 and mu_2")
dev.off()

# No. 4
wlist.avg = rep(0, N)
wlist.sum = rep(0, N)
for(i in 1:N){
  wlist.avg[i] = mean(wlist[[i]])
  wlist.sum[i] = sum(wlist[[i]])
}

mean(wlist.sum)
quantile(wlist.sum, c(0.025, 0.975))

mean(wlist.avg > 5/7)

# summary statistics
mean(MU1)
quantile(MU1, c(0.025, 0.975))

mean(MU2)
quantile(MU2, c(0.025, 0.975))

mean(SIG1)
quantile(SIG1, c(0.025, 0.975))

mean(SIG2)
quantile(SIG2, c(0.025, 0.975))
