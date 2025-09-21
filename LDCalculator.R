########################################
####################
# Function for calculating Linkage disequilibrium (LD) 

# Input parameters: chrn, snpchr, pop, Basenum
# chrn: Number of chromosomes 
# snpchr: Number of SNPs per chromosome (could be a vector if more than one chr exists)
# pop: Matrix of marker data with two rows per each individual
# Basenum: Number of individuals in historical population

LDCalculator <- function(chrn,snpchr,pop,Basenum)
{
  totsnp <- chrn*snpchr
  ldsVector <- as.numeric()
  for ( i in 1:(totsnp-1))
  {
    Sum=0
    for ( j in 1:(Basenum*2))
    {
      if((pop[j,i] + pop[j,i+1])==2){(Sum = Sum+1)}
    }
    FreqAB = Sum/(Basenum*2)
    Sumofc <- colSums(pop)
    FreqA = Sumofc[i]/(Basenum*2)
    FreqB = Sumofc[i+1]/(Basenum*2)
    d = FreqAB - (FreqA)*(FreqB)
    r2 = d^2 / (FreqA * FreqB * (1-FreqA) * (1-FreqB))
    ldsVector[i]<-r2
  }
  LD <- mean(ldsVector , na.rm=TRUE)
  return(LD)
}

########################################
