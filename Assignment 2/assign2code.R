## Precision parameters
N.graph = 1000
c.val = seq(0,0.5,length=N.graph);
N = 1000; # number of partitions


## Distribution parameters
a = 1.05;
b = 30;


## Initialize vectors
pel.val = rep(0,N.graph);
theta.val = seq(0,1,length=N); # all partitions
loss.val = rep(0,N);


## Calculate posterior expected loss values
for(k in 1:N.graph) {
  
  for(i in 1:N) {
    if(c.val[k] < theta.val[i]) {
      loss.val[i] = 10*abs(theta.val[i]-c.val[k]);
    }
    else loss.val[i] = abs(theta.val[i]-c.val[k]);
  }
  
  fun.val = loss.val*dbeta(theta.val, a, b);
  
  int.val = rep(0,N-1);
  
  for(j in 1:N-1) {
    int.val[j] = fun.val[j]/(N-1);
  }
  
  pel.val[k] = sum(int.val);
}


## Graph
png("pel.png")
plot(c.val, pel.val, type="l", lty=1, xlab="c", 
     ylab=expression(rho(c,x)), main="Posterior Expected Loss");
dev.off()