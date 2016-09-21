library(VarStar)
library(scales)

folder = "~/Dropbox/project/Astronomy/M33Var/CoreData/"
rdfiles = paste0(folder, "all_star_table.csv")

## The information table from OGLE
starTable = read.csv(rdfiles, header=TRUE)
starTable = subset(starTable, Type=="Mira")

## Get the k-th variable star
k=1 # 33 86 92 99 151 1128 1385
folder = "~/Dropbox/project/Astronomy/M33Var/CoreData/lightcurves/"
starID = as.character(starTable$ID[k])
rdfiles = paste0(folder, starID, ".dat")

## actual observations
starObs = read.table(rdfiles)

##True period
p0 = starTable$P_1[k]

## define new object of class "varStar"
## following three inputs should be MJD, mag, error
## prior mean level
meanMag = mean(starObs$V2)
x = new(gpModel,starObs$V1, starObs$V2 , starObs$V3, meanMag)
theta0 = rnorm(7,10,3)

## MUST set theta before MLE
x$set_period(p0,0)
x$gp_setTheta(theta0)
currL = x$gp_mLoglik(theta0)


for(j in 1:10){
    print(j)
    theta0[c(2,4)] = runif(2,300,1300)
    theta0[7] = rnorm(1,100,300)
    theta0[c(1,3,5)] = rnorm(3,5,2.2)
    theta0[6] = rnorm(1,2.5,1.1)
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
newX = seq(min(starObs$V1),max(starObs$V1),by=1)
plot(starObs[,1:2],pch=20)
newY = x$gp_predict(newX, 0)
lines(newX,newY$predy,pch=20,col="red",lwd=2)
newY = x$gp_predict(newX, 1)
newY2 = x$gp_predict(newX, 3)$predy-meanMag
lines(newX,newY$predy+newY2,
      col="blue",lwd = 2)
newY = x$gp_predict(newX, 2)
lines(newX,newY$predy,col="black",lty=2)


