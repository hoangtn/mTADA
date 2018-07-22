estimatePars <- function(pars, mcmcResult, nThin = NULL){
    mcmcDataFrame <- as.data.frame(mcmcResult)
        if (!is.null(nThin))
        mcmcDataFrame <- mcmcDataFrame[seq(1, dim(mcmcDataFrame)[1], by = nThin),]
    allPars <- colnames(mcmcDataFrame)

    if (is.null(pars)){
        message("\nYou are estimating all parameters, you can input a specific parameter(s) for pars\n")
        pars = allPars
    }

#    pars <- pars[grep("hyper|pi|alpha", pars)]
    message("====\nOnly pi, alpha and hyper parameters are estimated in this step\n",
            "gTADA does not calculate HPDs for hyper betas, just their medians\n===\n")

    if (length(pars[!is.element(pars, colnames(mcmcDataFrame))]) > 0)
        warning((pars[!is.element(pars, colnames(mcmcDataFrame))]), " is/are not in mcmc chains")
    pars <- pars[is.element(pars, colnames(mcmcDataFrame))]

    hpdList <- NULL
    for (iPar in pars) {
        message("Estimating for ", iPar)
        xData <-  mcmcDataFrame[, iPar]
        if ((sd(xData) > 10^-10) & (length(grep("Beta", iPar)) == 0))
            hpdList <- rbind(hpdList, loc1stats(xData)[1:3])
        else
            hpdList <- rbind(hpdList, rep(median(xData), 3))
    }
    rownames(hpdList) <- pars

    colnames(hpdList) <- c("Mode", "lCI", "uCI")

return(hpdList)
}
