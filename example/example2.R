## Use one fake light curve and its spectrum as example

library(varStar)
starObs = read.table("./mira_73762_wli31709_131.5.flc")
obsJD = starObs$V1
obsMag = starObs$V2
obsSigma = starObs$V3
head(starObs)


plotStarObs(obsJD,obsMag,obsSigma)


vsObj = new(simple_gpModel,
            obsJD, obsMag,obsSigma)
spc = vsObj$freq_est()


plot(spc[,1],spc[,2],type="n", main="",
     xlab=expression(Frequency(day^-1)),
    ylab="Log-Likelihood")
f0 = 0.006236
abline(v = f0, col = "blue")
abline(v = f0+1/365, col = "red", lty = 2)
abline(v = f0-1/365, col = "red", lty = 2)
lines(spc[,1],spc[,2])



