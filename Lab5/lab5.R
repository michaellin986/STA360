## Distribution parameters ##
b = 0.05
n = 50
y = 20
lam = 25

## Gibbs Sampling ##
N = 11000
b.samp = rep(NA,N)
b.samp[1] = b
n.samp = rep(NA,N)
n.samp[1] = n

for(i in 2:N){
  b.samp[i] = rbeta(1, y+1, n.samp[i-1]-y+1)
  n.samp[i] = y + rpois(1, (1-b.samp[i])*lam)
}

## Plot for first 10 ##
png("gibb10.png")
plot(b.samp[1:10], n.samp[1:10], pch = 20, xlab = "beta", 
     ylab = "N", main = "Gibbs Sampler For 10 Points")

for(i in 1:(9)){
  lines(c(b.samp[i],b.samp[i+1]), c(n.samp[i],n.samp[i]))
  lines(c(b.samp[i+1],b.samp[i+1]), c(n.samp[i],n.samp[i+1]))
}
dev.off()

## Plot for all points ##
png("gibb.png")
plot(b.samp, n.samp, pch = 20, xlab = "beta", 
     ylab = "N", main = "Gibbs Sampler")
dev.off()

## Credible Interval ##
print(quantile(b.samp, prob = c(0.05,0.95)))

## Burn-in discard ##
n.burn = n.samp[seq(1001,11000,1)]

## Probability N=20 ##
N.20 = length(n.burn[n.burn==20])
prob = N.20/length(n.burn)
print(prob)