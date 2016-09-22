fitComplexGP.single <- function(modelObj, theta0){
    modelObj$set_theta(theta0)
    opsR = optim(theta0,
                 modelObj$gp_mLoglik, 
                 modelObj$gp_DmLoglik, 
                 method="BFGS")
    return(opsR)
}

fitComplexGP.RandInit = function(){
    theta0 = rep(0, 7)
    theta0[c(2,4)] = runif(2,400,1000)
    theta0[1] = rnorm(1, 2.5, 1.2)
    theta0[3] = rnorm(1, 5.5, 2)
    theta0[5] =  runif(1, 0.1, 2)
    theta0[6] = rnorm(1,0.5,0.2)
    theta0[7] = rnorm(1,60,20)
    return(theta0)
}


fitComplexGP <- function(modelObj, initialTheta0 = NULL){
    if(is.null(initialTheta0)){
        currL = 1e100
        for(i in 1:5){
            theta0 = fitComplexGP.RandInit()
            opsR = fitComplexGP.single(modelObj, theta0)
            if(opsR$value<currL){
                currL = opsR$value
                optTheta = opsR$par
            }
        }
    }else{
        opsR = fitComplexGP.single(modelObj, theta0)
        optTheta = opsR$par
    }
    modelObj$set_theta(optTheta)
    return(modelObj)
}