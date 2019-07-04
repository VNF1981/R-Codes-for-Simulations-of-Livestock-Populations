########################################
####################
# Function for making pedigree-based relationship matrix 

# Input parameters: Vector of sires, Vector of dams 

ACreator <- function(Sire,Dam)
      {
        N <- length(Sire)
        A <- matrix(0, N, N)
        Sire <- (Sire == 0)*(N) + Sire
        Dam <- (Dam == 0)*N + Dam
        
        for(i in 1:N)
            {
              A[i,i] <- 1 + (A[Sire[i], Dam[i]]/2 )
              for(j in (i+1):N)
                  {
                   if (j > N) break
                   A[i,j] <- ( A[i, Sire[j]] + A[i,Dam[j]] )/2
                   A[j,i] <- A[i,j]
                  }
            }
        return(A[1:N,1:N])
      }
########################################