ps - periodic selection
==

### 1.1 The main question of this project: What is the probability of fixation of an allele under selection every P generations?

While there are many papers addressing randomly varying selection coefficients, all have assumed selection occuring every single generation. I argue that the dynamic of periodic selection can be quite different because the only force between periods is drift. Intuitively, if P is large then selection plays a neglible role in the dynamics of the allele.

### 1.2 Why is this interesting?
An interesting macro-evolutionary trend arises from the above intiution: longer lived organisms are more robust to periodic selection events than shorter lived organisms.

Suppose every 30 years a heatwave afflicts a region. Parrots - who have a long lifespan - see this event every generation. While, say, mice who have a lifespan of one year experience selection every 30 generations. Therefore for 29 generations, the only force affecting the dynamics of the allele conferring resistance to the heatwave is drift - resulting in a lower probability of fixation. 

Can this explain why rats cannot adapt to the flowering of bamboo: http://www.pbs.org/wgbh/nova/nature/rat-attack.html ?

### 2.1 Simulations
We write a few functions to simulate W-F model
```{r}
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
   p <- numeric(ngen)
   p[1] <- pstart
   for( i in 2:ngen){
     p[i] <- rbinom(1,2*popsize, pprimemut(p[i-1], w,u,v) )/(2*popsize)
   }
   p
}

evolvedrift_periodic <- function(ngen, pstart, w, popsize, period, u=0, v=0){
  n_w = c(1,1,1)
  p = evolvedrift(ngen, pstart, n_w, popsize, u, v)
  p = evolvedrift(ngen, p, w, popsize, u, v)
  return(p)
}

```