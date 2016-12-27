y = c(2,1,9,4,3,3,7,7,5,7)
x = seq(0,10,0.01)

## Define priors ##
n = 10
prior1 = dgamma(x, shape = 50, rate = 10)
prior2 = rep(0,length(x))
prior2[x<=3 | x>7] = 0
prior2[x>3 & x<=4] = 0.07
prior2[x>4 & x<=5] = 0.45
prior2[x>5 & x<=6] = 0.39
prior2[x>6 & x<=7] = 0.09

## Graph priors ##
png("prior.png")
plot(x, prior1,type="l", col='red', lwd = 2, main = "prior1=red, prior2=black", 
     xlab = expression(theta), ylab="prior")
lines(x, prior2, lwd = 2)
dev.off()

## Posteriors ##
post1 = dgamma(x, shape = 50+sum(y), rate = 10+n)
c = gamma(sum(y)+1)/(n^(sum(y)+1))*(0.07*(pgamma(4, sum(y)+1,n) - pgamma(3, sum(y)+1,n))
                                    +0.45*(pgamma(5, sum(y)+1,n) - pgamma(4, sum(y)+1,n))
                                    +0.39*(pgamma(6, sum(y)+1,n) - pgamma(5, sum(y)+1,n))
                                    +0.09*(pgamma(7, sum(y)+1,n) - pgamma(6, sum(y)+1,n))) 
post2.func = function(t) prior2*t^sum(y)*exp(-n*t)/c
post2 = post2.func(x)

## Graph posteriors ##
png("post.png")
plot(x,post1,type="l", col='red', lwd = 2, main = "post1=red, post2=black", 
     xlab = expression(theta), ylab="posterior")
lines(x,post2, lwd = 2)
dev.off()

## 95% CI ##
qpost1 = qgamma(c(0.025, 0.975), shape = 50+sum(y), rate = 10+n)