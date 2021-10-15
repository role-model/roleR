#include <Rcpp.h>

using namespace Rcpp;

// Base comm class
class commCpp {
    private:

    //note - private variables will not be inherited
    public:
        NumericVector abundance;
        NumericMatrix traits;
        int Smax;

        //constructor
        commCpp(NumericVector abundance_, NumericMatrix traits_, int Smax_) : abundance(abundance_), traits(traits_), Smax(Smax_)
        {
        }
};

// localComm subclass, which includes pi vector of genetic diversities
class localCommCpp : public commCpp {
    private:

    public:
        NumericVector pi;

        //constructor
        localCommCpp(NumericVector abundance_, NumericMatrix traits_, int Smax_,
                     NumericVector pi_) : commCpp(abundance_,traits_,Smax_) //: abundance(abundance_), traits(traits_), Smax(Smax_)
        {
            pi = pi_;
        }

        //birth function
        void birth(int i) {
            abundance[i] = abundance[i] + 1;
        }
};

// metaComm subclass
class metaCommCpp : public commCpp {
    private:

    public:

    //constructor
    metaCommCpp(NumericVector abundance_, NumericMatrix traits_,int Smax_) :
        commCpp(abundance_, traits_, Smax_)
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
