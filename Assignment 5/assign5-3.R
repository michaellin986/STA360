## Define paramters ##
mu.0 = 0
lambda.0 = 1/100^2
lambda = 1

N = 10^6
harm.mean = rep(0,5)
for(j in 1:5){
  theta.samp = rnorm(N, mean = mu.0, sd = sqrt(lambda.0^(-1)))
  
  harm = rep(0,N)
  for(i in 1:N) {
    harm[i] = 1/dnorm(2, mean = theta.samp[i], sd = sqrt(lambda^(-1)))
  }
  
  harm.mean[j] = sum(harm)/N
}

marg = 1/harm.mean

print(marg)

true.val = dnorm(2, mean = 0, sd = sqrt(lambda^-1+lambda.0^-1))
print(true.val)