## Load data ##
x = c(204, 215, 182, 225, 207, 188, 205, 227, 190, 211, 196, 203)
y = c(211, 233, 244, 241, 195, 252, 238, 249, 220, 213)

## Define parameters ##
a = 2.5     #(shape)
b = 0.01  #(rate)
n = length(x)
m = length(y)

## Define priors ##
H1.prior = 1/4
H0.prior = 3/4

## Marginal ##
w = rep(0, n) #w[i] = log(x[i]!)
for(i in 1:n){
  w[i] = sum(log(seq(1,x[i],1)))
}

v = rep(0, m)
for(i in 1:m){
  v[i] = sum(log(seq(1,y[i],1)))
}

f = sum(log(seq(1,a+sum(x)-1,1))) #f = log(gamma(a+sum(x)))
g = sum(log(seq(1,a+sum(y)-1,1))) #g = log(gamma(a+sum(y)))
h = sum(log(seq(1,a+sum(x)+sum(y)-1,1))) #h = log(gamma(a+sum(x)+sum(y)))

H1.loglik = 2*a*log(b) + f - 2*log(gamma(a)) - sum(w) -
  (a+sum(x))*log(b+n) + g - sum(v) - (a+sum(y))*log(b+m)
H0.loglik = a*log(b) + h - (a+sum(x)+sum(y))*log(n+m+b) - sum(w) - sum(v) - log(gamma(a))

## Posterior ##
H1.post = exp(H1.loglik)*H1.prior/(exp(H1.loglik)*H1.prior+exp(H0.loglik)*H0.prior)
H0.post = exp(H0.loglik)*H0.prior/(exp(H1.loglik)*H1.prior+exp(H0.loglik)*H0.prior)

## Bayes Factor ##
B10 = exp(H1.loglik - H0.loglik)

## Odds ##
odds.prior = H1.prior/H0.prior
odds.post = H1.post/H0.post
