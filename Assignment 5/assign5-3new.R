
## Harmonic Mean Function ##
# Courtesy of Radford Neal's Blog #
harmonic.mean.marg.lik <- function (x, s0, s1, n)
{ post.prec <- 1/s0^2 + 1/s1^2
  t <- rnorm(n,(x/s1^2)/post.prec,sqrt(1/post.prec))
  lik <- dnorm(x,t,s1)
  1/mean(1/lik)
}

## Define Variable ##
mu.0 = 0
lambda.0 = 1/100^2
lambda = 1

val = rep(0,5)
for(i in 1:5){
val[i]=harmonic.mean.marg.lik(2,sqrt(lambda.0^-1),sqrt(lambda),10^6)
}

print(val)