##### Initialization ######
## Load Package ##
library(MASS)
library(MCMCpack)

## Load Data ##
bdata = read.table("bluecrab.dat")
odata = read.table("orangecrab.dat")

## Working Variables ##
n = 50
N = 10000


##### Blue Crab #####
## Blue Crab Posterior Sample ##
ybar = c(mean(bdata[,1]), mean(bdata[,2]))
mu0 = ybar
L0 = cov(bdata)

nun = 4+n
S0 = cov(bdata)

## Sampler ##
sigma.temp = S0
sigma.blue = list()
theta.blue = NULL
cor.blue = rep(NA,N)

for(i in 1:N){
  ## Update theta.blue ##
  Ln = solve(solve(L0) + n*solve(sigma.temp))
  mun = Ln%*%(solve(L0)%*%mu0 + n*solve(sigma.temp)%*%ybar)
  theta.temp = mvrnorm(1, mun, Ln)
  theta.blue = rbind(theta.blue, theta.temp)
  
  ## Update sigma.blue ##
  Sn = S0 + (t(bdata)-theta.temp)%*%t((t(bdata)-theta.temp))
  sigma.temp = solve(rwish(nun, solve(Sn)))
  sigma.blue[[i]] = sigma.temp
  
  ## Correlation Coefficient ##
  cor.blue[i] = sigma.temp[1,2]/sqrt(sigma.temp[1,1]*sigma.temp[2,2])
}


##### Orange Crab #####
## Orange Crab Posterior Sample ##
ybar = c(mean(odata[,1]), mean(odata[,2]))
mu0 = ybar
L0 = cov(odata)

nun = 4+n
S0 = cov(odata)

## Sampler ##
sigma.temp = S0
sigma.or = list()
theta.or = NULL
cor.or = NULL

for(i in 1:N){
  ## Update theta.or ##
  Ln = solve(solve(L0) + n*solve(sigma.temp))
  mun = Ln%*%(solve(L0)%*%mu0 + n*solve(sigma.temp)%*%ybar)
  theta.temp = mvrnorm(1, mun, Ln)
  theta.or = rbind(theta.or, theta.temp)
  
  ## Update sigma.or ##
  Sn = S0 + (t(odata)-theta.temp)%*%t((t(odata)-theta.temp))
  sigma.temp = solve(rwish(nun, solve(Sn)))
  sigma.or[[i]] = sigma.temp
  
  ## Correlation Coefficient ##
  cor.or[i] = sigma.temp[1,2]/sqrt(sigma.temp[1,1]*sigma.temp[2,2])  
}

##### Plot #####
b = 1
png("bluetheta.png")
plot(theta.blue[b:10000,], pch = '.', col='blue', xlab = expression(theta[1]), ylab = expression(theta[2]), main = 'Blue Crab')
#points(theta.or[b:10000,], pch = '.', col='orange')
dev.off()

png("ortheta.png")
plot(theta.or[b:10000,], pch = '.', col='orange', xlab = expression(theta[1]), ylab = expression(theta[2]), main = 'Orange Crab')
dev.off()

##### Correlation Densities #####
## Blue ##
png("bluedist.png")
hist(cor.blue[b:10000], freq=F, xlab = expression(rho), 
     main = 'Blue Crab: Correlation Between Y1 and Y2')
lines(density(cor.blue[b:10000]), col = 'blue', lwd = 3)
dev.off()

## Orange ##
png("ordist.png")
hist(cor.or[b:10000], freq=F,xlab = expression(rho), 
     main = 'Orange Crab: Correlation Between Y1 and Y2')
lines(density(cor.or[b:10000]), col = 'orange', lwd = 3)
dev.off()

## Correlation Difference ##
print(mean(cor.blue < cor.or))
