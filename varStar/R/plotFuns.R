
plotStarObs = function(obsJD, obsMag, obsSigma, 
                       xlab = "Julian Date - 2450000", 
                       ylab = expression(paste(italic("I"), "  [mag]")), ...){
    ylim = c(max(obsMag), min(obsMag))
    plot(obsJD, obsMag, pch = 20, ylim = ylim, 
         xlab = xlab, ylab = ylab, type = "n", ...)
    arrows(obsJD, obsMag + obsSigma,
           obsJD, obsMag - obsSigma,
           col = "grey", angle = 90, length = 0.02, code = 3)
    points(obsJD, obsMag, pch = 20)
}

gpModel_JD_seq = function(modelObj){
    tmin = min(modelObj$JD)
    tmax = max(modelObj$JD)
    ttSeq = seq(tmin, tmax, length.out = 500)
    return(ttSeq)
}


plotAllComponents.complex = function(modelObj){
    par0 = par(mfrow=c(4,1))
    par(mar=c(2,4,1,1))
    xaxt = "n"
    tt = gpModel_JD_seq(modelObj)
    for(componentI in 0:3){
        if(componentI == 3) xaxt = "s"
        plotStarObs(modelObj$JD, modelObj$mag, modelObj$sigma, xaxt = xaxt)
        predy = modelObj$gp_predict(tt, componentI)
        lines(tt, predy)
    }
    ##mtext("Julian Date - 2450000", side = 1, line = 1, cex = 0.9)
    par(par0)
}

plotAllComponents = function(modelObj){
    plotAllComponents.complex(modelObj)
}
