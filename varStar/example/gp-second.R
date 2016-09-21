library(VarStar)
rm(list=ls())
folderTable = "~/M33Project/"
folder = "~/M33Project/lightcurves/"

## The information table from OGLE
rdfiles = paste0(folderTable, "all_star_table.csv")
starTable = read.csv(rdfiles, header=TRUE)
starTable = subset(starTable, Type=="Mira")

## Get the k-th variable star
k=1 # 33 86 92 99 151 1128 1385
starID = as.character(starTable$ID[k])
rdfiles = paste0(folder, starID, ".dat")

## actual observations
starObs = read.table(rdfiles)
nObs = dim(starObs)[1]
if(nObs>500){
    smSub = sample(1:nObs, 500,replace=FALSE)
    smSub = sort(smSub)
    starObs = starObs[smSub, ]
}
##True period
p0 = starTable$P_1[k]

## define new object of class "simple_gpModel"
## following three inputs should be MJD, mag, error
## prior: mean level
meanMag = mean(starObs$V2)
priorVec = c(meanMag, 10^2,
             3^2,
             0,0, #mu1 mu2
             100^2,100^2) #Sigma1 Sigma2
x = new(simple_gpModel, starObs$V1, 
        starObs$V2 , starObs$V3, priorVec)


## MUST set theta before MLE
x$set_period(p0,0)
theta0 = c(3.1,200.1)
x$gp_setTheta(theta0)
currL = x$gp_mLoglik(theta0);currL
dcurrL = x$gp_DmLoglik(theta0)
x$gp_mLoglik(theta0- 0.0001*dcurrL)

optRes = simple_gp_fit(x,jIter = 6)
optTheta = optRes[2:3]

for(j in 1:5){
    print(j)
    theta0[1]= runif(1,0.2,3.5)
    theta0[2] = runif(1,70,320)
    opsR = optim(theta0,x$gp_mLoglik, x$gp_DmLoglik, method="BFGS")
    print(opsR$value)
    if(opsR$value<currL){
        currL=opsR$value
        print(currL)
        optTheta = opsR$par
    }
}

##General flexible framework
##Draw back: multiple starting value
## not always get the desired shape
x$gp_setTheta(optTheta)
newX = seq(min(starObs$V1),max(starObs$V1),by=5)
newY = x$gp_predict(newX)
simple_gp_plot(starObs,newX,newY)
