library(VarStar)
starObs = read.table("./OGLE-LMC-LPV-00055.dat")
obsJD = starObs$V1
obsMag = starObs$V2
obsSigma = starObs$V3
head(starObs)


plotStarObs(obsJD,obsMag,obsSigma)


p0 = 290.9 ## true period of the light curve
modelObj = new(gpModel, obsJD, obsMag, obsSigma)
modelObj$set_freq(1/p0)

set.seed(100)
modelObj = fitComplexGP(modelObj)

plotAllComponents(modelObj)




