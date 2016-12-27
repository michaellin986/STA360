##### Initialization ######
## Load Package ##
library(MASS)
library(MCMCpack)

## Load Data ##
data = read.table('interexp.dat',header=T)
Y = data

##### Imputation #####
## Empirical estimates ##
theta.A = mean(na.omit(data$yA))
theta.B = mean(na.omit(data$yB))
var.A = var(na.omit(data$yA))
var.B = var(na.omit(data$yB))
cor.AB = cor(na.omit(data))
cov.AB = cov(na.omit(data))

## Imputation Calculation ##
n = dim(data)[1]

for(i in 1:n){
  # Impute A response
  if(is.na(data[i,1])){
    data[i,1] = theta.A + (data[i,2] - theta.B)*cor.AB[1,2]*sqrt(var.A/var.B)
  }
  
  # Impute B response
  if(is.na(data[i,2])){
    data[i,2] = theta.B + (data[i,1] - theta.A)*cor.AB[1,2]*sqrt(var.B/var.A)
  }
}

## T-Test ##
data.test = t.test(data[,1], data[,2], paired=T)

##### Gibbs Sampler #####
## Prior Parameters ##
N = 1000
n = dim(Y)[1]
p = dim(Y)[2]
mu0 = c(theta.A, theta.B)
sd0 = c(sqrt(var.A), sqrt(var.B))
L0 = cov(na.omit(Y))
S0 = L0
nu0 = p + 2

## Starting values ##
sigma.temp = S0
Y.full = Y
O = 1*(!is.na(Y))
for(i in 1:n){
  if(is.na(Y.full[i,1])){
    Y.full[i,1] = theta.A
  }
  
  if(is.na(Y.full[i,2])){
    Y.full[i,2] = theta.B
  }
}

theta = NULL
sigma = list()

## Sampler ##
for(i in 1:N){
  ## update theta ##
  ybar = apply(Y.full, 2, mean)
  Ln = solve(solve(L0) + n*solve(sigma.temp))
  mun = Ln %*% (solve(L0)%*%mu0 + n*solve(sigma.temp)%*%ybar)
  theta.temp = mvrnorm(1,mun,Ln)
  theta = rbind(theta, theta.temp)
  
  ## update sigma ##
  Sn = S0 + (t(Y.full)-c(theta.temp))%*%t(t(Y.full)-c(theta.temp))
  sigma.temp = solve(rwish(nu0 + n, solve(Sn)))
  sigma[[i]] = sigma.temp
  
  ## update missing data ##
  for(j in 1:n){
    Y.full.temp = O[j,]
    if(Y.full.temp[1]==0 & Y.full.temp[2]==1){
      mean.temp = theta.A + cov.AB[1,2]*(1/var.B)*(Y.full[j,2]-theta.B)
      sd.temp = sqrt(var.A - cov.AB[1,2]*(1/var.B)*cov.AB[1,2])
      Y.full[j,1]=rnorm(1, mean = mean.temp, sd = sd.temp)
    }
    
    if(Y.full.temp[1]==1 & Y.full.temp[2]==0){
      mean.temp = theta.B + cov.AB[1,2]*(1/var.A)*(Y.full[j,1]-theta.A)
      sd.temp = sqrt(var.B - cov.AB[1,2]*(1/var.A)*cov.AB[1,2])
      Y.full[j,2]=rnorm(1, mean = mean.temp, sd = sd.temp)
    }
  }  
}

## Posterior Mean and Conf. Int. ##
Y.full.test = t.test(Y.full[,1], Y.full[,2], paired=T)

plot(data[,1],data[,2], col = "red", pch = 20)
points(Y.full[,1], Y.full[,2], col = "blue", pch=20)
points(data[,1][data==Y.full],data[,2][data==Y.full], col = "purple", pch = 20)