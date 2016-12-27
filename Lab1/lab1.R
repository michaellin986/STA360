## Part 1
a.prior = 1/2;
b.prior = 1/2;
A.success = 11;
A.total = 14;
B.success = 5;
B.total = 6;

x.A = seq(0, 1, length=100);
hx.A = dbeta(x.A, a.prior+A.success, b.prior+A.total-A.success);

x.B = seq(0, 1, length=100);
hx.B = dbeta(x.B, a.prior+B.success, b.prior+B.total-B.success);

png("densityA.png");
plot(x.A, hx.A, type="l", lty=1, xlab=expression(theta), 
     ylab="Beta(11.5, 3.5)", main="Posterior Distribution A");
dev.off()

png("densityB.png");
plot(x.B, hx.B, type="l", lty=1, xlab=expression(theta),
     ylab="Beta(5.5, 1.5)", main="Posterior Distribution B");
dev.off()


## Part 2
prob.A = 1-pbeta(0.8, a.prior+A.success, b.prior+A.total-A.success);
prob.B = 1-pbeta(0.8, a.prior+B.success, b.prior+B.total-B.success);


## Part 3
numsample = 1000000
sample.A = rbeta(numsample, a.prior+A.success, b.prior+A.total-A.success);
sample.B = rbeta(numsample, a.prior+B.success, b.prior+B.total-B.success);
counter = 0
for(i in 1:numsample) {
  if (sample.A[i] < sample.B[i]) counter = counter + 1;
}

blah = counter/numsample;