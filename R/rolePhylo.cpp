#include <Rcpp.h>
#include <commCpp2.cpp>

using namespace Rcpp;

class rolePhylo {
    public:
        int n;
        NumericMatrix e;
        NumericVector l;
        LogicalVector alive;
        StringVector tipNames;
        long scale;

        //constructor
        rolePhylo(int n_, NumericMatrix e_,NumericVector l_,
                  LogicalVector alive_,StringVector tipNames_,long scale_)
                : n(n_), e(e_), l(l_), alive(alive_), tipNames(tipNames_), scale(scale_)
        {
        }
};

RCPP_EXPOSED_CLASS(roleModelCpp)

    RCPP_MODULE(modelCpp) {
        class_<roleModelCpp>("roleModelCpp")
        .constructor<localCommCpp,metaCommCpp>()
        .field("local", &roleModelCpp::localComm)
        .field("meta", &roleModelCpp::metaComm)
        ;
    }
