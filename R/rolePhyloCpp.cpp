#include <Rcpp.h>
#pragma once

using namespace Rcpp;

class rolePhyloCpp {
    public:
        int n;
        NumericMatrix e;
        NumericVector l;
        LogicalVector alive;
        StringVector tipNames;
        long scale;

        //constructor
        rolePhyloCpp(int n_, NumericMatrix e_,NumericVector l_,
                  LogicalVector alive_,StringVector tipNames_,long scale_)
                : n(n_), e(e_), l(l_), alive(alive_), tipNames(tipNames_), scale(scale_)
        {
        }
};

RCPP_EXPOSED_CLASS(rolePhyloCpp)

RCPP_MODULE(phyloCpp) {
    class_<rolePhyloCpp>("rolePhyloCpp")
        .constructor<int,NumericMatrix,NumericVector,LogicalVector,StringVector,long>()
        ;
}
