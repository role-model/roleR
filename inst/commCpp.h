//#ifndef COMM
//#define COMM
#include <Rcpp.h>

using namespace Rcpp;

class commCpp{
public:
    NumericVector abundance;
    NumericMatrix traits;
    int Smax;
    commCpp(NumericVector abundance_, NumericMatrix traits_, int Smax_);
};

class localCommCpp{
public:
    NumericVector abundance;
    NumericMatrix traits;
    int Smax;
    NumericVector pi;
    localCommCpp(NumericVector abundance_, NumericMatrix traits_, int Smax_, NumericVector pi_);
};

//#endif
