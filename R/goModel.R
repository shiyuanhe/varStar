fitComplexGP <- function(starObs, p0){
    x = new(gpModel, starObs$V1,
            starObs$V2,
            starObs$V3, mean(starObs$V2))
    x$set_freq(1/p0, 0)
    theta0 = exp(rnorm(7, 1, 4))
    x$gp_setTheta(theta0)
    currL = x$gp_mLoglik(theta0)
    for(i in 1:20){
        theta0[c(2,4)] = runif(2,400,1000)
        theta0[1] = rnorm(1, 2.5, 1.2)
        theta0[3] = rnorm(1, 5.5, 2)
        theta0[5] =  runif(1, 0.1, 2)
        theta0[6] = rnorm(1,0.5,0.2)
        theta0[7] = rnorm(1,60,20)
        
        opsR = optim(theta0,
                     x$gp_mLoglik, 
                     x$gp_DmLoglik, 
                     method="BFGS")
        if(opsR$value<currL){
            currL = opsR$value
            optTheta = opsR$par
        }
    }
    x$gp_setTheta(optTheta)
    return(x)
}