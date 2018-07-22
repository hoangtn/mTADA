#An example of gamma.mean.dn:
#1                   24.81551                         82.82128
#2                   24.81551                          1.00000
#3                    1.00000                         82.82128

#dataDN, muAll, Ndn should be across categories and 2 traits
#probs = NULL, BOTH, FIRST, SECOND

posProb.dn  <- function(dnData, muAll, gamma.mean.dn, Ndn, prob0, beta.dn){
    nGdn <- dim(dnData)[2]
    m <- dim(dnData)[1]
    bF <- array(0, dim = c(m, 4)) ##The first column is for NULL model
    ###############Calculate BF for all categories
    for (j in 1:nGdn){
        message("j = ", j)
###NULL
        marg.lik0 <- dpois(dnData[, j], 2*Ndn[j]*muAll[, j])
#        marg.lik0 <- dnbinom(dnData[, j], beta.dn[i, j], beta.dn[i, j]/(beta.dn[i, j] + 2*Ndn[j]*muAll[, j]))
        ###Both traits, first trait and then second trait
        marg.lik1 <- array(1, dim = c(m, 4))
#        bF[, 1] <-  bF[, 1] + log(marg.lik0) #For NULL model
        for (i in 1:4){ ##Add NULL
            if (gamma.mean.dn[i, j] == 1)
                marg.lik1[, i] <- dpois(dnData[, j], 2*Ndn[j]*muAll[, j])
            else
                marg.lik1[, i] <- dnbinom(dnData[, j], gamma.mean.dn[i, j]*beta.dn[i, j], beta.dn[i, j]/(beta.dn[i, j] + 2*Ndn[j]*muAll[, j]))
            bF[, i] <- bF[, i ] + log(marg.lik1[, i])# - log(marg.lik0)
            ##bF = 1:4, but calculate for 2:4
        }
        print(j)
        }
    ##########################
    bF <- exp(bF)
   
    PP <- t(apply(bF, 1, function(x)        x*prob0/sum(x*prob0)))
    colnames(PP) <- c("NO", "BOTH", "FIRST", "SECOND")
    return(list(BF = bF, PP = PP))

}
