data = read.table("divorce.dat")

x = data$V1
y = data$V2
l = length(x)

n = 1000
bval = rep(0,n)
cval = rep(0,n)
zval = list()

tb2=16
tc2=16

bval[1] = rnorm(1)
cval[1] = rnorm(1)
zval[[1]] = rnorm(l)

for(i in 2:n){
  # Sample for beta
  bval[i] = rnorm(1, mean = sum(x*zval[[i-1]])/(1/tb2+sum(x^2)), sd = sqrt((1/tb2 + sum(x^2))^-1))
  
  # Sample for c
  temp = sort(zval[[i-1]])
  maxList = temp[temp<cval[i-1]]
  maxList = sort(maxList, decreasing = T)
  if(is.na(maxList[1])){
    maxZ = -Inf  # maxZ is largest Z_i that is less than c
  }
  else{
  maxZ = maxList[1]
  }  
  
  minList = temp[temp>cval[i-1]]
  minList = sort(minList)
  if(is.na(minList[1])){
    minZ = Inf   # minZ is smallest Z_i that is greater than c
  }
  else{
    minZ = minList[1]
  } 
    
  umin = pnorm(maxZ, mean = 0, sd = sqrt(tc2))
  umax = pnorm(minZ, mean = 0, sd = sqrt(tc2))
  U = runif(1, umin, umax)
  cval[i] = qnorm(U, mean = 0, sd = sqrt(tc2))
  
  # Sample for z
  ztemp = rep(0,25)
  for(j in 1:25){
    if(y[j]==1){
      lbound = pnorm(cval[i], mean = bval[i]*x[j], sd = 1)
      U = runif(1, lbound, 1)
      ztemp[j] = qnorm(U, mean = bval[i]*x[j], sd = 1)
      
    }
    
    if(y[j]==0){
      ubound = pnorm(cval[i], mean = bval[i]*x[j], sd = 1)
      U = runif(1, 0, ubound)
      ztemp[j] = qnorm(U, mean = bval[i]*x[j], sd = 1)
    }
  }
  zval[[i]] = ztemp
  
}
nval = 1:n

## running averages
bavg = rep(0,n)
bavg[1] = bval[1]
for(k in 2:n){
  bavg[k] = (bavg[k-1]*(k-1) + bval[k])/k
}

cavg = rep(0,n)
cavg[1] = cval[1]
for(k in 2:n){
  cavg[k] = (cavg[k-1]*(k-1) + cval[k])/k
}

## plot
#png("bacf.png")
print(acf(bval))
#dev.off()

#png("acfc.png")
print(acf(cval))
#dev.off()

#png("acfz.png")
print(acf(zval[[1]]))
#dev.off()

#png("mixb3.png")
plot(1:n, bval, pch=20, main = "Traceplot with Running Avg. for Beta", xlab = "n", ylab = "beta")
for(i in 1:n-1){
  lines(c(nval[i],nval[i+1]),c(bavg[i],bavg[i+1]), lwd=3)
}
#dev.off()

#png("mixc3.png")
plot(1:n, cval, pch=20, main = "Traceplot with Running Avg. for c", xlab = "n", ylab = "c")
for(i in 1:n-1){
  lines(c(nval[i],nval[i+1]),c(cavg[i],cavg[i+1]),lwd=3)
}
#dev.off()

## Calculate other things
quantile(bval, probs=0.025)
quantile(bval, probs=0.975)

sum=0
for(i in 1:n){
  if(bval[i]>0)
    sum=sum+1
}
print(sum/length(bval))