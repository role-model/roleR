#pragma once

#include <RcppArmadillo.h>
#include "commCpp.cpp"
#include "rolePhyloCpp.cpp"
#include <string>

using namespace Rcpp;

//' @name roleDataCpp
//' @title a class to specify a frozen roleModel and its data - one index of the model timeseries
//' @field new Constructor
//' @field localComm the local community, an object of class localCommCpp
//' @field metaComm the metacommunity, an object of class metaCommCpp
//' @field phylo the phylogeny of the metacommunity, an object of class rolePhyloCpp

class roleDataCpp {
    public:
        localCommCpp local;
        metaCommCpp meta;
        rolePhyloCpp phylo;
        NumericVector stats; 
        int iter_num;
        
        // constructor
        roleDataCpp(localCommCpp local_, metaCommCpp meta_, 
                    rolePhyloCpp phy_) : 
                     local(local_), meta(meta_), 
                     phylo(phy_)
        {
        }
        
        roleDataCpp(){
        }
};

RCPP_EXPOSED_CLASS(roleDataCpp)

  