if(T){
folder = "~/Dropbox/project/Astronomy/M33Var/CoreData/"
rdfiles = paste0(folder, "all_star_table.csv")

## The information table from OGLE
starTable = read.csv(rdfiles, header=TRUE)
starTable = subset(starTable, Type=="Mira")

## Get the k-th variable star
k=35 # 33 86 92 99 151 1128 1385
folder = "~/Dropbox/project/Astronomy/M33Var/CoreData/lightcurves/"

fdNgt = "~/Dropbox/project/Astronomy/M33Var/inputfiles/"
fdOutput = "~/Moutput/M33"

sBegin = 100 ## beginning of the search region
sEnd = 1500  ## end of the search region
###first 25
rmC1 = c(6,145, 95,161,170,232)
rmC2 = c(165,156,122,165,179,181,196,202,204,209,221,236)
rmC = c(rmC1, rmC2)
for (k in 1:10) {
    if (k %in% rmC) next
    print(k)
    starID = starTable$ID[k]
    rdfiles = paste0(folder, starID, ".dat")
    
    ## actual observations
    starObs = read.table(rdfiles)
    ##True period
    p0 = starTable$P_1[k]
    
    ## define new object of class "varStar"
    ## following three inputs should be MJD, mag, error
    x = new(sinusoidModel,starObs$V1, starObs$V2 , starObs$V3)
    x$set_period(p0)
    x$model_fit()
    starModel = x$get_model()
    
    
    for(j in 1:25){
        fields=c(0:9,letters[1:19])
        oneField=sample(fields,1)
        ngtMJD=read.table(paste0(fdNgt,"m33i-",oneField,"-ngt.dat"))
        nObsWeight=read.table(paste0(fdNgt, "m0",oneField,"-obs-hst.dat"),
                              header=TRUE)
        nObs=with(nObsWeight, sample(Fraction,1,prob=frequency))
        if(nObs < 20 ) next
        nights=sample(ngtMJD[,1],nObs,replace=FALSE)
        nights=sort(nights)
        sData = x$get_fake(nights, 4.7)
        sData=round(sData,4)
         write.csv(sData,file=paste0(fdOutput,starID,"-",oneField,"-",j,".csv"),
                  row.names=FALSE,col.names=FALSE)
#          deltaMJD = max(starModel$MJD)-min(starModel$MJD)
#          deltaMJD = ceiling(deltaMJD/p0) - 1
#          fTime = (sData[,1]-min(sData[,1])+1666.43) %% (deltaMJD * p0) + min(starModel$MJD)
#          plotVS(starModel,FALSE)
#          points(fTime,sData[,2],col="red",pch=20)
#         
        
        #
        # Period estimation

        y = new(sinusoidModel,sData[,1], sData[,2], sData[,3])
        residSeq = y$period_est(sBegin, sEnd, 1000)
        y$set_period(p0)
        y$model_fit()
        starModel2 = y$get_model()
        ## residuals plot
        
        png(paste0(fdOutput,starID,"-",oneField,"-",j,".png"),width=500,height=700)
        par(mfrow=c(2,1))
        plotVS(starModel2,FALSE)
        luB = apply(residSeq[,-1], 1, quantile, probs=c(0.05,0.95))
        freqC = residSeq[,1]
        freqC = c(freqC, rev(freqC))
        perc = c(luB[1,], rev(luB[2,]))
        par0 = par(mar=c(3,4,1,1))
        plot(residSeq[,1], rowMeans(residSeq[,-1]), 
             xlab = "Frequency", ylab = "RSS",
             type="l",main=nObs,ylim=range(perc))
        polygon(freqC,perc,col="grey",border=NA)
        lines(residSeq[,1], rowMeans(residSeq[,-1]),col="red",lwd=2)
        abline(v=1/p0,col="blue",lwd = 3)
        dev.off()
    }
}

#stopCluster(cl)




