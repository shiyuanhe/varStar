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
}
