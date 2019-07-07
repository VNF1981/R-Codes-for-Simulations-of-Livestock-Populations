########################################
####################
# Random Mating algorithm used in our research

# Input parameters: selectedmales,selectedfemales,totsnp,generation,Lines,nRef,chrn,StpinMorgan

# selectedmales: Vector of identification numbers of selected males.              
# selectedfemales: Vector of identification numbers of selected females.
# totsnp: Total number of SNPs.
# generation: Indicates the index of the current generation in the generation loop. 
#             Used for construction of pedigree. 
#             If there is only one generation in evaluations, this parameter can be ommited. 
# Lines: Matrix of marker information (2n*m, n = number of individuals, m= number of markers).
# nRef: Number of individuals.
# chrn: Number of chromosomes.
# StpinMorgan: Distance between SNPs in Morgan.

RandommatN <- function(selectedmales,selectedfemales,totsnp,generation,Lines,nRef,chrn,StpinMorgan)
{
  
  N <- (nrow(selectedmales) + nrow(selectedfemales)) 
  SireNum <- nrow(selectedmales)    
  DamNum <- N-SireNum  
  Sirevectemp <- selectedmales[,1]
  Sirevectemp.2 <- rep(Sirevectemp,each=100)
  Damvectemp <- selectedfemales[,1]
  Damvectemp.2 <- rep(Damvectemp,each=10)
  Sirevec <- sample(Sirevectemp.2,nRef,replace=FALSE) 
  Damvec <- sample(Damvectemp.2,nRef,replace=FALSE)          
  valpop <- matrix(0,nRef*2,totsnp)
  ped <- matrix(0,nRef,3)
  colnames(ped) <- c("ID","Sire","Dam")
  
  for (indiv in 1:nRef)
      {
        gameteIndex1 <- (indiv-1)*2+1
        gameteIndex2 <- gameteIndex1+1
        ## Sires gamete
        SireINDEX <- Sirevec[indiv]                        
        SireGamet1 <- (((Sirevec[indiv])-1)*2)+1                 
        SireGamet2 <- (SireGamet1)+1                  
        valpop[gameteIndex1,] <- Recombin(Lines[SireGamet1,],
                                          Lines[SireGamet2,],
                                          chrn,
                                          totsnp,
                                          StpinMorgan)
        ## Dams gamete
        DamINDEX <- Damvec[indiv]                           
        DamsGamet1 <- (((Damvec[indiv])-1)*2)+1                   
        DamsGamet2 <- (DamsGamet1)+1                   
        valpop[gameteIndex2, ] <- Recombin(Lines[DamsGamet1,],
                                           Lines[DamsGamet2,],
                                           chrn,
                                           totsnp,
                                           StpinMorgan)
          
        ped[,"ID"] <- (1+(nRef*(generation-1))):(nRef*generation)
        ped[indiv,"Sire"] <- (Sirevec[indiv]+((generation-2)*nRef))
        ped[indiv,"Dam"] <- (Damvec[indiv]+((generation-2)*nRef))
      } 
  
  matout <- list(valpop=valpop , ped=ped)
  return(matout)
} # End of function