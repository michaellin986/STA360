## Monte Carlo Method ##

# Prior Definition
theta.0 = 36.07
s.0 = 0.02

# Other Definition
s = 0.0002
N = 10^6
x.samp = c(36.077916, 36.078032, 36.078129, 36.078048, 
           36.077942, 36.089612, 36.077789, 36.077563)

# Draw theta's
theta.samp = rcauchy(N, location = theta.0, scale = s.0)


prob = rep(0,length(x.samp))
marg = rep(0,N)
marg.sum = rep(0,N)

for (j in 1:N) {
  prob = dcauchy(x.samp, location = theta.samp[j], scale = s)
  marg[j] = prod(prob)
}

marg.sum[1] = marg[1]
for (i in 2:N) {
  marg.sum[i] = marg.sum[i-1] + marg[i]
}

for (k in 1:N){
  marg.sum[k] = marg.sum[k]/k
}

x = 1:N
plot(x, marg.sum, log = "x", type="l", lty=1, col="blue", xlab="N", 
     ylab="approx marginal distribution", main="Monte Carlo marginal distribution")




## Importance sampling
# Define "q" parameters
theta.q = median(x.samp)
s.q = 10^-4

# sample theta's from "q"
theta.qsamp = rcauchy(N, location = theta.q, scale = s.q)

# calculate likelihood of theta's for "p" and "q"
p.like = dcauchy(theta.qsamp, location = theta.0, scale = s.0)
q.like = dcauchy(theta.qsamp, location = theta.q, scale = s.q)

marg.imp = rep(0,N)
for (j in 1:N) {
  prob = dcauchy(x.samp, location = theta.qsamp[j], scale = s)
  marg.imp[j] = prod(prob)
}

marg.impsum = rep(0,N)
marg.impsum[1] = marg.imp[1]*p.like[1]/q.like[1]
for (i in 2:N) {
  marg.impsum[i] = marg.impsum[i-1] + marg.imp[i]*p.like[i]/q.like[i]
}

for (k in 1:N) {
  marg.impsum[k] = marg.impsum[k]/k
}

lines(x, marg.impsum, log = "x", type="l", lty=1, col="red")

legend(10^2.5,3.5*10^18, c("Monte Carlo","Importance Sampling"), 
       lty=c(1,1), 
       lwd=c(2.5,2.5),col=c("blue","red"))