library(MCMCpack)
set.seed(1)
## load data ##
data = read.table("http://www.stat.washington.edu/people/pdhoff/Book/Data/hwdata/azdiabetes.dat", sep='', header=T)
data = data[,c(1:7)]
data$intercept = rep(1,nrow(data))
y = data[,2]
X = as.matrix(data[,-2])

## initialize parameters ##
g = length(y); nu0 = 2; s20 = 1;
S = 1000;
n = dim(X)[1]; p = dim(X)[2];

## compute beta ##
Hg = (g/(g+1)) * X%*%solve(t(X)%*%X)%*%t(X)
SSRg = t(y)%*%( diag(1, nrow=n) - Hg) %*%y

s2 = 1/rgamma(S, (nu0+n)/2, (nu0*s20+SSRg)/2)
Vb = g*solve(t(X)%*%X)/(g+1)
Eb = Vb%*%t(X)%*%y

E = matrix(rnorm(S*p, 0, sqrt(s2)), S, p)
beta = t( t(E%*%chol(Vb)) + c(Eb))

## posterior (beta) confidence interval ##
beta.npreg.ci = quantile(beta[,1], c(0.025, 0.975))
beta.bp.ci = quantile(beta[,2], c(0.025, 0.975))
beta.skin.ci = quantile(beta[,3], c(0.025, 0.975))
beta.bmi.ci = quantile(beta[,4], c(0.025, 0.975))
beta.ped.ci = quantile(beta[,5], c(0.025, 0.975))
beta.age.ci = quantile(beta[,6], c(0.025, 0.975))
beta.int.ci = quantile(beta[,7], c(0.025, 0.975))

## function: compute marginal probability ##
lpy.X = function(y, X, g=length(y), nu0=1, s20=try(summary(lm(y~-1+X))$sigma^2, silent=T)){
  n = dim(X)[1]; p = dim(X)[2];
  if(p==0){
    Hg = 0
    s20 = mean(y^2)
  }
  if(p>0){
    Hg = (g/(g+1)) * X%*%solve(t(X)%*%X)%*%t(X)
  }
  SSRg = t(y)%*%( diag(1, nrow=n) - Hg) %*%y
  
  -0.5*(n*log(pi)+p*log(1+g)+(nu0+n)*log(nu0*s20+SSRg)-nu0*log(nu0*s20)) +
    lgamma((nu0+n)/2) - lgamma(nu0/2)
}

## MCMC setup ##
z = rep(1, dim(X)[2])
lpy.c = lpy.X(y, X[, z==1, drop=F])
S = 1000
Z = matrix(NA, S, dim(X)[2])

## Gibbs sampler ##
for(s in 1:S){
  for(j in sample(1:dim(X)[2])){
    zp = z
    zp[j] = 1-zp[j]
    lpy.p = lpy.X(y, X[, zp==1, drop=F])
    r = (lpy.p - lpy.c)*(-1)^(zp[j]==0)
    z[j] = rbinom(1, 1, 1/(1+exp(-r)))
    if(z[j]==zp[j]){
      lpy.c = lpy.p
    }
  }
  
  Z[s,] = z
}


## find new beta ##
BETA = matrix(0, S, p)
for(s in 1:S){
  z = Z[s,]
  X.z = NULL
  for(i in 1:p){
    if(z[i]!=0){
      X.z = cbind(X.z, X[,i])
    }
  }
  
  Hg.new = (g/(g+1)) * X.z%*%solve(t(X.z)%*%X.z)%*%t(X.z)
  SSRg.new = t(y)%*%( diag(1, nrow=n) - Hg.new) %*%y
  
  s2.new = 1/rgamma(1, (nu0+n)/2, (nu0*s20+SSRg.new)/2)
  Vb.new = g*solve(t(X.z)%*%X.z)/(g+1)
  Eb.new = Vb.new%*%t(X.z)%*%y
  
  E.new = matrix(rnorm(1*ncol(X.z), 0, sqrt(s2.new)), 1, ncol(X.z))
  beta.new = t( t(E.new%*%chol(Vb.new)) + c(Eb.new))
  
  for(j in 1:sum(z)){
    BETA[s, which(z!=0)[j] ] = beta.new[j]
  }
  
}

## Pr(beta!=0|y) ##
Z.mean = rep(0,p)

for(i in 1:p){
  Z.mean[i] = mean(Z[,i])
}

## new posterior(BETA) confidence interval ##
BETA.npreg.ci = quantile(BETA[,1], c(0.025, 0.975))
BETA.bp.ci = quantile(BETA[,2], c(0.025, 0.975))
BETA.skin.ci = quantile(BETA[,3], c(0.025, 0.975))
BETA.bmi.ci = quantile(BETA[,4], c(0.025, 0.975))
BETA.ped.ci = quantile(BETA[,5], c(0.025, 0.975))
BETA.age.ci = quantile(BETA[,6], c(0.025, 0.975))
BETA.int.ci = quantile(BETA[,7], c(0.025, 0.975))

## new BETA1 ##
BETA1 = beta*Z
quantile(BETA1[,1], c(0.025, 0.975))
quantile(BETA1[,2], c(0.025, 0.975))
quantile(BETA1[,3], c(0.025, 0.975))
quantile(BETA1[,4], c(0.025, 0.975))
quantile(BETA1[,5], c(0.025, 0.975))
quantile(BETA1[,6], c(0.025, 0.975))
quantile(BETA1[,7], c(0.025, 0.975))