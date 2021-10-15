#include <Rcpp.h>
#include <commCpp.cpp>
#include <roleModelCpp.cpp>

using namespace Rcpp;

// [[Rcpp::export]]

localCommCpp birthLocal(localCommCpp x, int i)
{
    x.abundance[i] += 1;
    return x;
}

roleModelCpp birthRole(roleModelCpp x)
{
    //sample a species for birth relative to local abundance
    //0 vs 1 start indices may cause problems
    NumericVector probs = x.localComm.abundance[Rcpp::Range(1,x.localComm.Smax)];
    IntegerVector i = sample(x.localComm.Smax, 1, false, probs);

    x.localComm = birthLocal(x.localComm, i[0]);
    return x;
}

//todo - permit calling birth() rather than birthL/R, however this may not be
//possible thru modules
RCPP_MODULE(birthCpp) {
    function("birthL" , &birthLocal);
    function("birthR" , &birthRole);
}

