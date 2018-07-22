posProb.dnOneTrait <- function(dnData, muAll, gamma.mean.dn,
                                Ndn, prob0, beta.dn){
    nGdn <- dim(dnData)[2]
    m <- dim(dnData)[1]
    bF <- array(0, dim = c(m, 2))
    ###############Calculate BF for all categories
    for (j in 1:nGdn){
        message("j = ", j)
###NULL
        marg.lik0 <- dpois(dnData[, j], 2*Ndn[j]*muAll[, j])
###Both traits, first trait and then second trait
        marg.lik1 <- dnbinom(dnData[, j], gamma.mean.dn[j]*beta.dn[j],
                                  beta.dn[j]/(beta.dn[j] + 2*Ndn[j]*muAll[, j]))
        bF[, j] <- log(marg.lik1) - log(marg.lik0)
        
        }        
    ##########################
    bF <- exp(bF)
    bF <- apply(bF, 1, prod)
     PP <- bF1*prob0/(1 - prob0 + prob0*bF1)
    
    return(list(BF = bF, PP = PP))

}
