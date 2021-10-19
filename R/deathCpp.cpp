#include <Rcpp.h>
#include <commCpp.cpp>
#include <roleModelCpp.cpp>
#include <rolePhyloCpp.cpp>

using namespace Rcpp;

// [[Rcpp::export]]

localCommCpp deathLocal(localCommCpp x, int i)
{
    x.abundance[i] -= 1;
    return(x);
}

rolePhyloCpp deathPhylo(rolePhyloCpp x, int i)
{
    // set tip to dead
    x.alive[i] = false;

    return(x);
}

roleModelCpp deathRole(roleModelCpp x)
{

    // sample a species for death proportional to species abundance
    NumericVector probs = x.localComm.abundance[Rcpp::Range(1,x.localComm.Smax)];
    IntegerVector i = sample(x.localComm.Smax, 1, false, probs);

    x.localComm = deathLocal(x.localComm, i[0]);

    // if death led to extinction, call death on rolePhylo
    if(x.localComm.abundance[i] <= 0)
    {
        x.phylo = deathPhylo(x.phylo, i[0]);
    }

    return(x);
}

//todo - permit calling birth() rather than birthL/R, however this may not be
//possible thru modules
RCPP_MODULE(deathCpp) {
    function("deathL" , &deathLocal);
    function("deathR" , &deathRole);
}
