########################################
####################
# Function for simulating meiosis & recombination

# Input parameters: 
# hap: Vectors of marker information of an individual
# chrn: Number of chromosomes (Must have equal lengths)
# totsnp: Total number of SNP markers (equal number of SNPs per chromosome are considered)
# StpinMorgan: Distance between SNPs in Morgan (equal distances are considered)

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
########################################