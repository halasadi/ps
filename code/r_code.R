# pprime gives the frequency after selection
# w[1], w[2], w[3] are w11, w12, w22, the fitnesses of
#    the three genotypes. 
#   p is the frequency before selection.

pprime <- function(p, w){
  wbar <- w[1]*p*p + w[2]*2*p*(1-p) + w[3]*(1-p)*(1-p)
  ( w[1]*p*p + p*(1-p)*w[2] )/ wbar    
}

# similar to pprime but with mutation, A1 -> A2 at rate
# u,   A2 -> A1 at rate v.

pprimemut <- function(pin, w, u=0, v=0){
  p = pin
  wbar <- w[1]*p*p + w[2]*2*p*(1-p) + w[3]*(1-p)*(1-p)
  pp = ( w[1]*p*p + p*(1-p)*w[2] )/ wbar 
  ppp = pp*(1-u) + (1-pp)*v
  return(ppp)   
}

# evolve population with selection and drift (and 
#    mutation.)  Returns a vector of the frequency trajectory.

evolvedrift <- function( ngen, pstart, w,popsize, u=0, v=0){
  #p <- numeric(ngen)
  #p[1] <- pstart
  if (ngen <=0){
    return(pstart)
  }
  p = pstart
  for( i in 1:ngen){
    #p[i] <- rbinom(1,2*popsize, pprimemut(p[i-1], w,u,v) )/(2*popsize)
    p <- rbinom(1,2*popsize, pprimemut(p, w,u,v) )/(2*popsize)
  }
  return(p)
}

# period = 1 recapitulates the evolvedrift function above with ngen = nperiods
# Requires: period >=1
evolvedrift_periodic <- function(nperiods, pstart, w, popsize, period, u=0, v=0){
  p = pstart

  for (i in 1:nperiods){
    p = evolvedrift(1, p, w, popsize)
    p = evolvedrift(period-1, p, w = c(1,1,1), popsize)
    if (p < 1e-10 || p > (1-1e-10)){
      return(p)
    }
  }
  return(p)
}

# N diploid organisms
# s is selective advantage of allele A
kimura <- function(N, s){
  return((1-exp(-2*s))/(1-exp(-4*N*s)))
}

periods = seq(1, 100, 1)
N = 1e4
p0 = 1/(2*N)
s=0.5
w = c(1+s, 1 + (s/2), 1)

p_fixation = rep(0, length(periods))
for (i in 1:length(periods)){
  period = periods[i]
  nreps = 1000
  cnt = 0
  for(j in 1:nreps){
    nperiods = 10000 # should be enough
    p = evolvedrift_periodic(nperiods, p0, w, popsize = 2*N, period)
    if (p >0.999){
      cnt = cnt + 1
    }
  }
  p_fixation[i] = cnt/nreps
  
}
kimura_values = rep(0, length(periods))
for (i in 1:length(periods)){
  kimura_values[i] = kimura(2*N, s/periods[i])
}

plot(periods, p_fixation, xlab = "period", ylab = "probability of fixation",
       col = "black", main = "selective force occuring every period generations", lwd=4)
points(periods, kimura_values, type="l", col = "red", lwd=4)
legend("topright", c("predictions using ''effective'' selection coefficient", "simulations"),
           lty=c(1,1), lwd=c(4,4,4),col=c("red", "black"), cex=1)
