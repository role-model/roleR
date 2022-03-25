
#include <RcppArmadillo.h>
using namespace Rcpp;

#include "commCpp.cpp"

RCPP_MODULE(commCpp) {
    class_<metaCommCpp>("metaCommCpp")
        .constructor<NumericVector,NumericMatrix>()
        .field("abundance", &metaCommCpp::abundance)
        .field("traits", &metaCommCpp::traits)
    ;
    class_<localCommCpp>("localCommCpp")
        .constructor<NumericVector,NumericVector,NumericVector,int>()
    
        .field("abundance_indv", &localCommCpp::abundance_indv)
        .field("species_ids", &localCommCpp::species_ids)
        .field("traits", &localCommCpp::traits)
        .field("J", &localCommCpp::J)
        .field("S_index", &localCommCpp::S_index)
        .field("pi_sp", &localCommCpp::pi_sp)
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
    .constructor<int,NumericMatrix,NumericVector,LogicalVector,StringVector,double>()
    .field("alive", &rolePhyloCpp::alive)
    .field("n", &rolePhyloCpp::n)
    .field("e", &rolePhyloCpp::e)
    .field("l", &rolePhyloCpp::l)
    .field("tipNames", &rolePhyloCpp::tipNames)
    .field("scale", &rolePhyloCpp::scale)
    .method("birth", &rolePhyloCpp::birth)
    .method("death", &rolePhyloCpp::death)
    .method("speciation", &rolePhyloCpp::speciation)
    .method("immigration", &rolePhyloCpp::immigration)
    ;
}

#include "roleModelCpp.cpp"

RCPP_MODULE(modelCpp) {
    class_<roleModelCpp>("roleModelCpp")
    .constructor<localCommCpp,metaCommCpp,rolePhyloCpp,List>()
    .field("local", &roleModelCpp::local)
    .field("meta", &roleModelCpp::meta)
    .field("phylo", &roleModelCpp::phylo)
    .field("params", &roleModelCpp::params)
    .field("print", &roleModelCpp::print)
    .field("print_vectors", &roleModelCpp::print_vectors)
    .field("timeseries", &roleModelCpp::timeseries)
    .method("birth", &roleModelCpp::birth)
    .method("death", &roleModelCpp::death)
    .method("speciation", &roleModelCpp::speciation)
    .method("immigration", &roleModelCpp::immigration)
    .method("copyData", &roleModelCpp::copyData)
    ;
}


#include "roleDataCpp.cpp"

RCPP_MODULE(dataCpp) {
  class_<roleDataCpp>("roleDataCpp")
  .constructor<localCommCpp,metaCommCpp,rolePhyloCpp>()
  .field("local", &roleDataCpp::local)
  .field("meta", &roleDataCpp::meta)
  .field("phylo", &roleDataCpp::phylo)
  .field("stats", &roleDataCpp::stats)
  .field("iter_num", &roleDataCpp::iter_num)
  ;
}

