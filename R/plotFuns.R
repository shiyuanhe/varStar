
plotStarObs = function(starObs, 
                       xlab = "Julian Date - 2450000", 
                       ylab = NULL, ...){
    if(is.null(ylab)){
        ylab = expression(paste(italic("I"), "  [mag]"))
    }
    ylim = c(max(starObs$V2), min(starObs$V2))
    plot(starObs$V1, starObs$V2, pch = 20, ylim = ylim, 
         xlab = xlab, ylab = ylab, type = "n", ...)
    arrows(starObs$V1, starObs$V2 + starObs$V3,
           starObs$V1, starObs$V2 - starObs$V3,
           col = "grey", angle = 90, length = 0.02, code = 3)
    points(starObs$V1, starObs$V2, pch = 20, ylim = ylim)
}


