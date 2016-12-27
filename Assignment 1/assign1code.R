x.data=c(20.9, 69.7, 3.6, 21.8, 21.4, 0.4, 6.7, 10.0);
a=0.1; b=1  #prior
n=length(x.data); y=sum(x.data)  #data

x.prior=seq(0, 0.2, length=100);
hx.prior=dgamma(x.prior, shape=a, rate=b);
x.post=seq(0, 0.2, length=100);
hx.post=dgamma(x.post, shape=a+n, rate=b+y)

png("prior.png");
plot(x.prior, hx.prior, type="l", lty=1, xlab=expression(theta), 
     ylab="Gamma(0.1, 1.0)", main="Prior Distribution");
dev.off()

png("post.png");
plot(x.post, hx.post, type="l", lty=1, xlab=expression(theta),
     ylab="Gamma(8.1, 155.5)", main="Posterior Distribution");
dev.off()