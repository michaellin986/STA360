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
xseq <- seq(from=0,to=1,by=0.01)

d.values <- sin(pi*xseq)^2

plot(xseq,d.values,type="l")

#### Uniform envelope function ####

# We can see from the previous plot that the
# density seems to have its maximum at 1. Since
# the uniform density on (0,1) is constant at 1,
# we will not have to worry the envelope being
# smaller than the target density.
# In order to make things slightly more abstract
# we will write a new function which will evaluate
# sin(x*pi) for us, and use this function throughout
# the rest of the code.

density.function <- function(x){
  return(sin(pi*x)^2)
}


samples.uniform <- NULL
N=1000
for ( i in 1:N ) {
  proposal <- runif(1) # Here we get a proposal value
  density.ratio <- density.function(proposal)/dunif(proposal)  # We calculate the ratio of the densities
  if ( runif(1) < density.ratio ) samples.uniform <- c(samples.uniform,proposal) # If a random uniform is lower than our ratio, we accept our sample, otherwise we reject. Then we repeat this process
}

hist(samples.uniform,freq=FALSE)
print(paste("Acceptance Ratio: ",length(samples.uniform)/N))

# If we wanted to, we could now increase the
# sample size to make the estimate more accurate.
# Alternatively, we could try to find a better
# proposal density.



#### Beta proposal ####

# First let's plot the density and a beta(2,2)
# density to see whether the beta is bigger at
# all points.

plot(xseq,density.function(xseq),type="l",col="blue",ylim=c(0,1.5))
lines(xseq,dbeta(xseq,2,2),col="red")

# We see that the beta is in fact bigger at all
# points, but it is in fact *too* big. We want
# the envelope density to be as close to the
# target density as possible. Let's multiply this by
# a factor less than 1 to bring it down a bit.

M=0.7
plot(xseq,density.function(xseq),type="l",col="blue",ylim=c(0,1.5))
lines(xseq,M*dbeta(xseq,2,2),col="red")

# This looks much better. For real applications we'd
# probably have to do much better still, especially if
# we're working with more complicated densities. Note that
# when we compute the ratio of the densities for the
# acceptance rate, we will now take the factor of M=0.7
# into account as well, so it will become
# p(x)/(M*g(x))

samples.beta <- NULL

for ( i in 1:N ) {
  proposal  = rbeta(1,2,2) # We obtain a proposal from a beta(2,2) random variable
  density.ratio = density.function(proposal)/(M*dbeta(proposal,2,2)) # We calculate the density ratio, including the multiplicative constant
  if (runif(1) < density.ratio ) samples.beta <- c(samples.beta,proposal)
}

hist(samples.beta)
print(paste("Acceptance Ratio: ",length(samples.beta)/N))
  
# As we can see the acceptance ratio has gone up from
# about 50% to about 90%. In real applications the
# acceptance ratio is probably going to be much lower,
# and a lot of the work (should you choose to do rejection
# sampling) can be in finding an envelope function that
# doesn't just give you a basically 0% chance of accepting
# a point.

