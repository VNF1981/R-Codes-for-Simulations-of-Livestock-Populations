########################################
####################
# Function for solving MME

# Input parameters: 
# Z: Incidence matrix (r*n), r = number of records, n = number of animals
# y: Vector of phenotypes
# rel: Relationship matrix
# alpha: (1 - h2)/h2

BLUP <- function(Z,rel,y,alpha)
{
  AInv <- solve(rel)
  #AInv <- chol2inv(chol(rel))       # Another way to inverse the matrix 
  k <- AInv*alpha
  ZZ <- crossprod(Z)
  w <- ZZ+k
  Zy <- crossprod(Z,y)
  bhat <- solve(w,Zy)
  bhat <- round(bhat,4)             
  return(bhat)
}
########################################