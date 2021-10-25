#include <Rcpp.h>
#pragma once

using namespace Rcpp;

class commCpp {
private:

public:
    NumericVector abundance;
    NumericMatrix traits;
    int Smax;

    //constructor
    commCpp(NumericVector abundance_, NumericMatrix traits_, int Smax_)
            : abundance(abundance_),traits(traits_), Smax(Smax_)
    {
    }
};

//' @name localCommCpp
//' @title a C++ class to specify the local community
//' @field new Constructor
//' @param abundance a numeric vector of abundances for each species
//' @param traits matrix of traits; the first column specifies the species index
//' (i.e. the index of that species in the \code{abundance} vector) and the
//' subsequent columns specify the trait values
//' @param Smax a single integer specifying the total number of species ever
//' recorded in the local community (both locally extinct and extant)
//' @details Smax is used for internal bookkeeping.  The dimensions of
//' \code{abundance} and \code{traits} can be greater than
//' \code{Smax}.  In such cases, \code{Smax} is used to index where to begin
//' adding new information (e.g. a new species can be added at the index
//' \code{Smax + 1}).

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
