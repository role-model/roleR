#include <Rcpp.h>
#include <commCpp.cpp>
#include <rolePhyloCpp.cpp>
#include <roleParamsCpp.cpp>

using namespace Rcpp;

//' @name roleModelCpp
//' @title a C++ class to specify the entire RoLE model
//' @field new Constructor
//' @field localComm the local community, an object of class localCommCpp
//' @field metaComm the metacommunity, an object of class metaCommCpp

class roleModelCpp {
    public:
        localCommCpp localComm;
        metaCommCpp metaComm;
        rolePhyloCpp phylo;
        roleParamsCpp params;

        //constructor
        roleModelCpp(localCommCpp local_, metaCommCpp meta_, rolePhyloCpp phy_,
                     roleParamsCpp params_) : localComm(local_), metaComm(meta_),
                     phylo(phy_), params(params_)
        {
        }
};

// expose C++ class to R through Rcpp
RCPP_EXPOSED_CLASS(roleModelCpp)

// the Rcpp module required to load using loadModule() to use the class in R
RCPP_MODULE(modelCpp) {
    class_<roleModelCpp>("roleModelCpp")
        .constructor<localCommCpp,metaCommCpp,rolePhyloCpp,roleParamsCpp>()
        .field("local", &roleModelCpp::localComm)
        .field("meta", &roleModelCpp::metaComm)
    ;
}
