N <- 10^6
l <- 1
l_0 <- 1/10^2

harm <- function(x, n) {
  ts <- rnorm(n, 0, sqrt(l_0^-1))
  1 / (sum(1/dnorm(x, ts, sqrt(l^-1))) / n)
}
print(harm(2, N))
