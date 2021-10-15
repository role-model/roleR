#include <Rcpp.h>
#pragma once

using namespace Rcpp;

// localComm class
class localCommCpp {
    private:

    public:
        NumericVector abundance;
        NumericMatrix traits;
        int Smax;
        NumericVector pi;

        //constructor
        localCommCpp(NumericVector abundance_, NumericMatrix traits_, int Smax_,
                     NumericVector pi_) : abundance(abundance_),
                     traits(traits_), Smax(Smax_), pi(pi_)
        {
        }

        //test birth function
        void birth(int i) {
            abundance[i] = abundance[i] + 1;
        }
};

// metaComm class
class metaCommCpp {
    private:

    public:
        NumericVector abundance;
        NumericMatrix traits;
        int Smax;

        //constructor
        metaCommCpp(NumericVector abundance_, NumericMatrix traits_,int Smax_) :
            abundance(abundance_), traits(traits_), Smax(Smax_)
        {
        }
};

RCPP_EXPOSED_CLASS(localCommCpp)
RCPP_EXPOSED_CLASS(metaCommCpp)

RCPP_MODULE(commCpp) {
    class_<metaCommCpp>("metaCommCpp")
        .constructor<NumericVector,NumericMatrix,int>()
        .field("abundance", &metaCommCpp::abundance)
        .field("traits", &metaCommCpp::traits)
        .field("Smax", &metaCommCpp::Smax)
    ;
    class_<localCommCpp>("localCommCpp")
        .constructor<NumericVector,NumericMatrix,int,NumericVector>()
        .field("abundance", &localCommCpp::abundance)
        .field("traits", &localCommCpp::traits)
        .field("Smax", &localCommCpp::Smax)
        .field("pi", &localCommCpp::pi)
    ;
}
