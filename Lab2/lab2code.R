size = 10000;

x = rep(0,size);
tausq = rgamma(size, shape=0.5, rate=0.5);

for(i in 1:size){
  x[i] = rnorm(1, 0, sqrt(1/(tausq[i])));
}


hist(x, breaks=35);

ks.test(x, "pt", 1);

probs = rep(0, 1000)
for (i in 1:1000) {
  y = rt(50, df = 1);
  probs[i] = ks.test(y, 'pt', 1)$p;
}
hist(probs, freq = FALSE);
