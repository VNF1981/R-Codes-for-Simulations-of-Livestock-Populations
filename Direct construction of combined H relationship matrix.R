########################################
####################
# Function for making combined H relationship matrix following Legarra et al, (2009) 

# Input parameters: 
# A: Pedigree-based relationship matrix for all individuals
# A22: Pedigree-based relationship matrix for genotyped individuals
# G22: Genomic relationship matrix for genotyped individuals
# A12: A matrix including pedigree-based relationships between un-genotyped and genotyped individuals
# Alphadiag: Difference between the average diagonal elements of A22 and G22 (used for mitigating bias in predictions)
# Alphaofdiag: Difference between the average off-diagonal values of A22 and G22 (used for mitigating bias in predictions) 
# GA.ID: Id of genotyped animals  
# UGA.ID: Id of un-genotyped animals

HCreator <- function(A,A22,G22,A12,Alphadiag,Alphaofdiag,GA.ID,UGA.ID)
{
  A21 <- t(A12)
  Mix1 <- A12 %*% (solve(A22))
  #Mix1 <- A12 %*% (chol2inv(chol(A22)))                        # Another way to inverse the matrix
  Mix1 <- round(Mix1,4)
  G22t <- (0.95*G22) + (0.05*A22)                               # To avoid singularity following VanRaden (2008)                     
  G22temp <- G22t                                               
  diag(G22temp) <- 0
  G22temp <- G22temp+Alphaofdiag
  diag(G22temp) <- 0                                            
  diag(G22temp) <- diag(G22t)+Alphadiag
  G22temp <- round(G22temp,4)                                   # Gw22, to mitigate bias following Vitezica et al, (2011)
  Mix2 <- (G22temp-A22)
  Mix2 <- round(Mix2,4)                 
  Mix3 <- (solve(A22)) %*% A21
  #Mix3 <- (chol2inv(chol(A22))) %*% A21
  Mix3 <- round(Mix3,4)
  mIden <- matrix(0,length(GA.ID),length(GA.ID))                # Identity matrix with dimensions equal to number of genotyped individuals
  diag(mIden) <- 1 

  vv <- matrix(0,nrow(A),ncol(A))                               
  vv[GA.ID,UGA.ID] <- Mix3
  vv[GA.ID,GA.ID] <- mIden
  ss <- matrix(0,length(GA.ID),nrow(A))                                
  ss[,GA.ID] <- mIden
  qq <- matrix(0,nrow(A),length(GA.ID))                                                                  
  qq[GA.ID,] <- mIden
  nn <- matrix(0,nrow(A),ncol(A)) 
  nn[UGA.ID,GA.ID] <- Mix1
  nn[GA.ID,GA.ID] <- mIden   
  
  Adelta <- (nn %*% (qq %*% (Mix2 %*% (ss %*% vv))))
  Adelta <- round(Adelta,4)
  H <- A + Adelta
  H <- round(H,4)                  
  return(H)                        
}
########################################