
#include <RcppArmadillo.h>
using namespace Rcpp;

#include "commCpp.cpp"

RCPP_MODULE(commCpp) {
    class_<metaCommCpp>("metaCommCpp")
        .constructor<NumericVector,NumericMatrix,int>()
    
        .field("abundance", &metaCommCpp::abundance)
        .field("traits", &metaCommCpp::traits)
        .field("Smax", &metaCommCpp::Smax)
    ;
    class_<localCommCpp>("localCommCpp")
        .constructor<NumericVector,NumericVector,int,NumericVector>()
    
        .field("abundance_indv", &localCommCpp::abundance_indv)
        .field("species_ids", &localCommCpp::species_ids)
        .field("traits", &localCommCpp::traits)
        .field("Imax", &localCommCpp::Imax)
        .field("Smax", &localCommCpp::Smax)
        .field("pi", &localCommCpp::pi)
        .field("traitdiffs", &localCommCpp::traitdiffs)
        .field("abundance_sp", &localCommCpp::abundance_sp)
        .field("traits_sp", &localCommCpp::traits_sp)
        .field("print", &localCommCpp::print)
    
        .method("birth", &localCommCpp::birth)
        .method("death", &localCommCpp::death)
        .method("speciation", &localCommCpp::speciation)
        .method("immigration", &localCommCpp::immigration)
    ;
}

#include "rolePhyloCpp.cpp"

RCPP_MODULE(phyloCpp) {
    class_<rolePhyloCpp>("rolePhyloCpp")
    .constructor<int,NumericMatrix,NumericVector,LogicalVector,StringVector,long>()
    .field("alive", &rolePhyloCpp::alive)
    .field("n", &rolePhyloCpp::n)
    .field("e", &rolePhyloCpp::e)
    .field("l", &rolePhyloCpp::l)
    .method("birth", &rolePhyloCpp::birth)
    .method("death", &rolePhyloCpp::death)
    .method("speciation", &rolePhyloCpp::speciation)
    .method("immigration", &rolePhyloCpp::immigration)
    ;
}

#include "roleModelCpp.cpp"

RCPP_MODULE(modelCpp) {
    class_<roleModelCpp>("roleModelCpp")
    .constructor<localCommCpp,metaCommCpp,rolePhyloCpp,roleParamsCpp>()
    .field("local", &roleModelCpp::localComm)
    .field("meta", &roleModelCpp::metaComm)
    .field("phylo", &roleModelCpp::phylo)
    .field("params", &roleModelCpp::params)
    .method("birth", &roleModelCpp::birth)
    .method("death", &roleModelCpp::death)
    .method("speciation", &roleModelCpp::speciation)
    .method("immigration", &roleModelCpp::immigration)
    .field("print", &roleModelCpp::print)
    ;
}

#include "roleParamsCpp.cpp"

RCPP_MODULE(paramsCpp) {
    class_<roleParamsCpp>("roleParamsCpp")
    .constructor<paramValuesCpp,std::string,int>()
    .field("values", &roleParamsCpp::values)
    ;
}

#include "paramValuesCpp.cpp"

RCPP_MODULE(paramValuesCpp) {
    class_<paramValuesCpp>("paramValuesCpp")
    .constructor()
    .field("speciation_meta", &paramValuesCpp::speciation_meta)
    .field("extinction_meta", &paramValuesCpp::extinction_meta)
    .field("trait_sigma", &paramValuesCpp::extinction_meta)
    .field("species_meta", &paramValuesCpp::species_meta)
    .field("individuals_meta", &paramValuesCpp::individuals_meta)
    .field("individuals_local", &paramValuesCpp::individuals_local)
    .field("dispersal_prob", &paramValuesCpp::dispersal_prob)
    .field("speciation_local", &paramValuesCpp::speciation_local)
    .field("dispersal_prob", &paramValuesCpp::dispersal_prob)
    //.constructor<double, double, double, double,double, double, double, double, double, double, double, double>()
    ;
}
