#include <RcppArmadillo.h>
#include "commCpp.cpp"
#include "rolePhyloCpp.cpp"
#include "roleParamsCpp.cpp"
#include <string>
#pragma once

using namespace Rcpp;

//' @name roleDataCpp
//' @title a class to specify a frozen roleModel and its data - one index of the model timeseries
//' @field new Constructor
//' @field localComm the local community, an object of class localCommCpp
//' @field metaComm the metacommunity, an object of class metaCommCpp
//' @field phylo the phylogeny of the metacommunity, an object of class rolePhyloCpp
//' @field iter_num the iteration of the model held in this object
//' 
class roleDataCpp {
    public:
        localCommCpp localComm;
        metaCommCpp metaComm;
        rolePhyloCpp phylo;
        int iter_num; 
        
        // constructor
        roleDataCpp(localCommCpp local_, metaCommCpp meta_, rolePhyloCpp phy_,
                     int iter_num_) : localComm(local_), metaComm(meta_),
                     phylo(phy_), iter_num(iter_num_)
        {
        }
};

RCPP_EXPOSED_CLASS(roleDataCpp)
