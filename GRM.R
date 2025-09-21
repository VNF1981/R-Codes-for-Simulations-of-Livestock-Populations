########################################
####################
# Function for making Genomic Relationship Matrix using VanRaden method

# Input parameters: Z matrix 
# Z: Matrix of genotypes(n*m), n = number of individuals, m = number of markers.                  
# elements of the matrix are 0,1, and 2, representing the number of major allele

GCreator <- function(Z)
     {
        for(v in 1:ncol(Z))
            {frq <-(sum(Z[,v])/(2*nrow(Z)))
             if(frq < 0.05 | frq > 0.95 ){Z[,v]<- NA}}
             snpremove <- colSums(is.na(Z)) != nrow(Z)
             Z <- Z[,snpremove]  # Omitting Markers with MAF < 0.05      
        p <- apply (Z,2,function(x) sum(x)/(length(x)*2))
        pt2 <- 2*p
        M <- t(apply(Z,1,function(x) x-pt2)) # centered matrix
        k <- 2*(sum(p*(1-p)))
        ZZ <- tcrossprod(M)
        GRM <- ZZ/k
        GRM <- round(GRM,4)
        return(GRM)
     }

########################################
