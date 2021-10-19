#include <Rcpp.h>
#include <commCpp.cpp>
#include <roleModelCpp.cpp>

using namespace Rcpp;

// [[Rcpp::export]]

localCommCpp immigrationLocal(localCommCpp x, int i)
{
    x.abundance[i] += 1;
    return x;
}

roleModelCpp immigrationRole(roleModelCpp x)
{
    //sample a species for birth relative to local abundance
    //0 vs 1 start indices may cause problems
    NumericVector probs = x.metaComm.abundance[Rcpp::Range(1,x.metaComm.Smax)];
    IntegerVector i = sample(x.metaComm.Smax, 1, false, probs);

    x.localComm = birthLocal(x.localComm, i[0]);
    return x;
}

//todo - permit calling birth() rather than birthL/R, however this may not be
//possible thru modules
RCPP_MODULE(immigrationCpp) {
    function("immigrationL" , &immigrationLocal);
    function("immigrationR" , &immigrationRole);
}

