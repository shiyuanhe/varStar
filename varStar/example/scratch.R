folder = "~/Dropbox/project/Astronomy/M33Var/CoreData/"
rdfiles = paste0(folder, "all_star_table.csv")

## The information table from OGLE
starTable = read.csv(rdfiles, header=TRUE)
starTable = subset(starTable, Type=="SRV")


## Get the k-th variable star
k=1 # 33 86 92 99 151 1128 1385
folder = "~/Dropbox/project/Astronomy/M33Var/CoreData/lightcurves/"


p0 = 1
mag = 1
nObs = 1
for (k in 1:1200) {
    print(k)
    starID = starTable$ID[k]
    rdfiles = paste0(folder, starID, ".dat")
    
    ## actual observations
    starObs = read.table(rdfiles)
    nObs[k] = dim(starObs)[1]
    ##True period
    p0[k] = starTable$P_1[k]
    mag[k] = median(starObs$V2)
}

plot(log(p0), mag, pch=19)
