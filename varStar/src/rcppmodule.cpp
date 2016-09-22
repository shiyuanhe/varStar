#include "gpmodel.h"
#include "simple_gpmodel.h"


RCPP_MODULE(varStar_mod){
    class_<varStar>("varStar")
    .constructor<NumericVector, NumericVector, NumericVector>()
    .field("JD", &varStar::MJD)
    .field("mag", &varStar::mag)
    .field("sigma", &varStar::mag_sigma)
    ;
    
    class_ <gpModel>( "gpModel")
        .derives<varStar>("varStar")
        .constructor<NumericVector, NumericVector, NumericVector>()
        .method("set_freq", &gpModel::set_freq)
        .method("set_theta", &gpModel::gp_setTheta)
        .method("gp_mLoglik", &gpModel::gp_mLoglik)
        .method("gp_DmLoglik", &gpModel::gp_DmLoglik)
        .method("gp_predict", &gpModel::gp_predict)
    //.method("get_fake", &gpModel::get_fake)
    ;
    
    class_ <simple_gpModel>( "simple_gpModel")
        .derives<varStar>("varStar")
        .constructor<NumericVector, NumericVector, NumericVector>()
        .method("set_theta", &simple_gpModel::set_theta)
        .method("set_freq", &simple_gpModel::set_freq)
        .method("gp_mLoglik", &simple_gpModel::gp_mLoglik)
        .method("predict", &simple_gpModel::predict)
        .method("BFGSopt", &simple_gpModel::BFGSopt)
        .method("freq_est", &simple_gpModel::freq_est)
    ;
    
}
