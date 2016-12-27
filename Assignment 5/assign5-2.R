N = 10^6

samp = rcauchy(N,0,1)

samp.sum = rep(0,N)
samp.mean = rep(0,N)

samp.sum[1] = samp[1]
for(i in 2:N){
  samp.sum[i] = samp.sum[i-1] + samp[i]
  samp.mean[i] = samp.sum[i]/i
}
samp.mean[1] = samp.sum[1]

x = 1:N
plot(x,samp.mean,log = "x", type="l", lty=1, xlab="N", ylab="Expected Value",
     main="Monte Carlo Approximation")