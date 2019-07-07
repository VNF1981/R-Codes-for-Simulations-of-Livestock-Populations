########################################
####################
# Mate allocation algorithm used in our research

# Input parameters: selectedmales,selectedfemales,totsnp,generation,Lines,nRef,Arel,chrn,StpinMorgan

# selectedmales: Vector of identification numbers of selected males.              
# selectedfemales: Vector of identification numbers of selected females.
# totsnp: Total number of SNPs.
# generation: Indicates the index of the current generation in the generation loop. 
#             Used for construction of pedigree. 
#             If there is only one generation in evaluations, this parameter can be ommited. 
# Lines: Matrix of marker information (2n*m, n = number of individuals, m= number of markers).
# nRef: Number of individuals.
# Arel: Could be pedigree-based, genomic, or combined H relationship matrices.
# chrn: Number of chromosomes.
# StpinMorgan: Distance between SNPs in Morgan. 

mateselN <- function(selectedmales,selectedfemales,totsnp,generation,Lines,nRef,Arel,chrn,StpinMorgan)
    {
      N <- (nrow(selectedmales) + nrow(selectedfemales)) 
      SireNum <- nrow(selectedmales)    
      DamNum <- N-SireNum  
      sexratio <- DamNum/SireNum 
      malesID <- selectedmales[,1]
      femalesID <- selectedfemales[,1]
      mID <- malesID
      fID <- femalesID
      valpop <- matrix(0,nRef*2,totsnp)
      ped <- matrix(0,nRef,3)
      colnames(ped) <- c("ID","Sire","Dam")
      MAmat <- matrix(0,SireNum,DamNum)
      rownames(MAmat) <- malesID
      colnames(MAmat) <- femalesID
      for(i in 1:nrow(MAmat))                   # MAmat with the corresponding elements in A or H
      {                                            # Sires in rows
        for(j in 1:ncol(MAmat))
        {                                            # Dams in columns
          A <- selectedmales[i,1]
          B <- selectedfemales[j,1]
          MAmat[i,j] <- (Arel[A,B])
          MAmat[i,j] <- MAmat[i,j]+(runif(1)/100)
        }
      }# 
      
      MM <- matrix(0,DamNum,2)
      colnames(MM) <- c("male","female")
      
      for(i in 1:SireNum) # for all males
      {
        relvec <- MAmat[i,] # i index
        for(j in 1:sexratio) # for total number of matings per male = DamNum/SireNum
        {
          maletemind <- ((i-1)*sexratio) + j     
          MM[maletemind,"male"] <- mID[i]                  
          yindex <- which(relvec[]==min(relvec[]),arr.ind=TRUE) 
          MM[maletemind,"female"] <- fID[yindex] 
          relvec <- relvec[-yindex]
          if(i*j < N-SireNum){MAmat <- MAmat[,-yindex]}
          fID <- fID[-yindex]
        } # end of matings for each male 
      } # end for all males
      
      Sirevectemp <- MM[,1]
      Sirevec <- rep(Sirevectemp,each=10) 
      Damvectemp <- MM[,2]
      Damvec <- rep(Damvectemp,each=10)  
      
      for (indiv in 1:nRef)
            {
              gameteIndex1 <- (indiv-1)*2+1
              gameteIndex2 <- gameteIndex1+1
              ## Sires gamete
              SireINDEX <- Sirevec[indiv]                        
              SireGamet1 <- (((Sirevec[indiv])-1)*2)+1                ## 1st gammete of sire
              SireGamet2 <- (SireGamet1)+1                              ## 2nd gammete of sire
              valpop[gameteIndex1,] <- Recombin(Lines[SireGamet1,],
                                                Lines[SireGamet2,],
                                                chrn,
                                                totsnp,
                                                StpinMorgan)
              ## Dams gamete
              DamINDEX <- Damvec[indiv]                           
              DamsGamet1 <- (((Damvec[indiv])-1)*2)+1                 ## 1st Gamete of dam
              DamsGamet2 <- (DamsGamet1)+1                              ## 2nd Gamete of dam
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
########################################