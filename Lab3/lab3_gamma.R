##### Lab 4 #####

# This file will work through 2 simple implementations
# of rejection sampling. We will sample from a density
# p(x) ~ sin(pi*x)^2, x in (0,1)
# The first implementation will be using a uniform
# envelope function, the second one using a beta
# envelope (since we're on the unit interval). We
# will also call the envelope function a proposal
# function, and denote its density by g(x)

#### Plot the density ####
xseq <- seq(from=0,to=10,by=0.01)

d.values <- xseq^4*(15*exp(-(xseq/2)^5)+5/81*exp(-(xseq/6)^5))

a = 3; b = 0.8;
thresh = 500

plot(xseq,d.values,type="l", ylim=c(0,100))
lines(xseq, thresh*dgamma(xseq, shape = a, rate = b), col = "red")

#### Uniform envelope function ####

# We can see from the previous plot that the
# density seems to have its maximum at 1. Since
# the uniform density on (0,10) is constant at 10,
# we will not have to worry the envelope being
# smaller than the target density.
# In order to make things slightly more abstract
# we will write a new function which will evaluate
# above distribution for us, and use this function throughout
# the rest of the code.

density.function <- function(x){
  return(x^4*(15*exp(-(x/2)^5)+5/81*exp(-(x/6)^5)))
}


samples.gamma <- NULL
N=10000
for ( i in 1:N ) {
  proposal <- rgamma(1, shape = a, rate = b) # Here we get a proposal value
  density.ratio <- density.function(proposal)/(thresh*dgamma(proposal, a, b))  # We calculate the ratio of the densities
  if ( runif(1) < density.ratio ) samples.gamma <- c(samples.gamma,proposal) # If a random uniform is lower than our ratio, we accept our sample, otherwise we reject. Then we repeat this process
}

hist(samples.gamma,freq=FALSE)
print(paste("Acceptance Ratio: ",length(samples.gamma)/N))

# If we wanted to, we could now increase the
# sample size to make the estimate more accurate.
# Alternatively, we could try to find a better
# proposal density.


