# ## find the optimal value of theta
# ## by maximizing CV-loglik
# 
# simple_gp_fit=function(simple_x, theta2_init,
#                        trial_init = NULL){
#     currL = 1e50
#     theta0 = c(0,0)
#     optTheta = c(0,0)
#     
#     t0 = 10^seq(-0.8,1.6,length.out=3)
#     t0 = c(t0,-t0)
#     t1 = theta2_init
#     for(j in 1:length(t0)){
#         for(k in 1:length(t1)){
#             theta0[1]= t0[j]
#             theta0[2] = t1[k]
#             simple_x$gp_setTheta(theta0)
#             tmpL = simple_x$gp_mLoglik(theta0);
#             try({
#                 opsR = optim(theta0, fn = simple_x$gp_mLoglik,
#                              gr = simple_x$gp_DmLoglik, method="BFGS")
#                 if(opsR$value<currL){
#                     currL=opsR$value
#                     optTheta = opsR$par
#                 }
#             })
#         }
#     }
#     
#     if(!is.null(trial_init)){
#         for(j in 1:dim(trial_init)[2]){
#             theta0[1] = trial_init[1,j]
#             theta0[2] = trial_init[2,j]
#             try({
#                 opsR = optim(theta0, 
#                              fn = simple_x$gp_mLoglik,
#                              gr = simple_x$gp_DmLoglik, 
#                              method="BFGS")
#                 if(opsR$value<currL){
#                     currL=opsR$value
#                     optTheta = opsR$par
#                 }
#             })
#         }
#     }
#     return(list(currL = currL,optTheta = optTheta))
# }


simple_gp_plot=function(starObs, newX, newY){
    plot(starObs[,1:2],pch=20,
         ylim=rev(range(starObs[,2])),
         ylab="Mag",xlab="MJD",type="n")
    bandX = c(newX,rev(newX))
    diagT = abs(diag(newY$predy_cov))
    newY_upper = newY$predy+1.643*sqrt(diagT) 
    newY_lower = newY$predy-1.643*sqrt(diagT) 
    bandY = c(newY_upper,rev(newY_lower))
    polygon(bandX,bandY,col="grey",border = NA)
    lines(newX,newY$predy,pch=20,col="red",lwd=3)
    sinPart= (newY$Ht_star)%*%(newY$gamma_bar)
    lines(newX,sinPart,col="blue",lwd = 2)
    points(starObs[,1:2],pch=20)
}


