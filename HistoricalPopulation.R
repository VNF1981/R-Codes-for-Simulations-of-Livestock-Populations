##########################################################################################################
################################################# Base Pop ###############################################
Basepop <- function(totqtl = 15 ,
                    chrn = 3 ,
                    lenchr = 1 ,
                    totsnp = 300 ,
                    Basenum = 100 ,
                    G = 50)
{
  
##################################   Definition of the map parameters
  lenchrCM <- lenchr*100
  GenomeL <- lenchrCM * chrn
  snpchr <- totsnp/chrn 
  snpchr <- floor(snpchr)
  qtlnum <- totqtl / chrn
  qtlnum <- floor(qtlnum)
  Dsnp <- snpchr/lenchrCM    # Number of SNPs per each cM
  Stp <- 1/Dsnp              # distances between snps
  Stp <- round(Stp,2)
  StpinMorgan <- Stp/100     # Distance between SNPs per each Morgan (Indicates the probability of CO)  
################################# Assigning QTL's
  qtl.ids <- seq (from = 1 , to = totsnp , by = (totsnp/totqtl))
################################# Assigning additive effects to the QTL's
  eff.add <- rnorm(length(qtl.ids),0,1)  
################################# Simulation of the founder population
  pop <- matrix(0,(Basenum * 2),totsnp)
  for (i in 1:nrow(pop))
  {
    for(j in 1:ncol(pop))
    {
      x <- runif(1,min = 0, max = 1) 
      if(x>=0.5){pop[i,j] <- 1} else{pop[i,j] <- 0}    # The probability of each SNP allele at the first generation is set to 0.5
    }
  } 
    
  random.mate.pop.temp <- matrix(0,Basenum*2,totsnp)
  LDs <- as.numeric()
  random.mate.pop <- pop
  SireNum<-Basenum/2
  for(generation in 1 : G) ## loop over generations
      {
        for(indiv in 1 : Basenum) ## loop over individuals
            {
              gameteIndex1<- (indiv-1)*2+1
              gameteIndex2<-gameteIndex1+1
              ## gamete 1 from Sire
              Sire<-sample( 1:SireNum,1) #Select one Sire Randomely
              SireGamet1<-(Sire-1)*2+1   #Finding 1st Gamet of Sire in order to Sire ID
              SireGamet2<-SireGamet1+1   #Finding 2nd Gamet of Sire in order to Sire ID
              random.mate.pop.temp[gameteIndex1,] <- Recombin(random.mate.pop[SireGamet1,],
                                                              random.mate.pop[SireGamet2,],
                                                              chrn,
                                                              totsnp,
                                                              StpinMorgan)
              ## gamete 2 From Dam
              Dams<-sample((SireNum+1): Basenum,1) #Select one Dam Randomely
              DamsGamet1<-(Dams-1)*2+1
              DamsGamet2<-DamsGamet1+1
              random.mate.pop.temp[gameteIndex2,] <- Recombin(random.mate.pop[DamsGamet1,],
                                                              random.mate.pop[DamsGamet2,],
                                                              chrn,
                                                              totsnp,
                                                              StpinMorgan)
              
            } ## End of Basenum loop
        random.mate.pop <- random.mate.pop.temp
        #cat("LD",LDs[g],"\n")
      } ## End of generation loop
  
  LDs <- LDCalculator(chrn,snpchr,random.mate.pop,Basenum)

  BaseOut<- list (random.mate.pop=random.mate.pop, LDs=LDs, eff.add=eff.add, qtl.ids=qtl.ids)
  return(BaseOut)
}  
##########################################################################################################
#################################################### LD Calculator #######################################
LDCalculator <- function(chrn,snpchr,random.mate.pop,Basenum)
{
  totsnp <- chrn*snpchr
  ldsVector <- as.numeric()
  for ( i in 1:(totsnp-1))
      {
        Sum=0
        for ( j in 1:(Basenum*2))
            {
             if((random.mate.pop[j,i] + random.mate.pop[j,i+1])==2){(Sum = Sum+1)}
            }
        FreqAB = Sum/(Basenum*2)
        Sumofc <- colSums(random.mate.pop)
        FreqA = Sumofc[i]/(Basenum*2)
        FreqB = Sumofc[i+1]/(Basenum*2)
        d = FreqAB - (FreqA)*(FreqB)
        r2 = d^2 / (FreqA * FreqB * (1-FreqA) * (1-FreqB))
        ldsVector[i]<-r2
      }
  LD <- mean(ldsVector , na.rm=TRUE)
  return(LD)
}
##########################################################################################################
#################################################### Recombination #######################################
Recombin <- function(hap1,hap2,chrn,totsnp,StpinMorgan)
{
  snpchr <- totsnp/chrn 
  recnum <- snpchr-1                                        
  matdat <- rbind(hap1,hap2) 
  recombinedhaps <- array(0,c(2,snpchr,chrn))
  for(i in 1:chrn)
      {
        firstind <- ((i-1)*snpchr)+1
        secondind <- i*snpchr
        tempmat <- matdat[,firstind:secondind]                
        for(j in 1:recnum)
            {
              x <- runif(1, min=0, max=1)
                  if(x <= StpinMorgan){x1 <- tempmat[1,(j+1):snpchr]
                  x2 <- tempmat[2,(j+1):snpchr]
                  tempmat[1,(j+1):snpchr] <- x2
                  tempmat[2,(j+1):snpchr] <- x1} 
            } # End of j loop
        recombinedhaps[,,i] <- tempmat 
      } # End of i loop
  matdat <- matrix(rbind(recombinedhaps),2,totsnp)
  recombinedhap1 <- matdat[1,]
  recombinedhap2 <- matdat[2,]
  gamete <- matrix(0,1,totsnp)
  x <- runif(1,min=0,max=1)
  if(x<0.25){gamete <- hap1}
  if(x>0.75){gamete <- hap2}
  if(x >= 0.25 && x <=0.5){gamete <- recombinedhap1}
  if(x > 0.5 && x <=0.75){gamete <- recombinedhap2}
  return(gamete)
}# End of function
##########################################################################################################
##########################################################################################################
