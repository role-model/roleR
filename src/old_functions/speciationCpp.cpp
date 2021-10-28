#pragma once

#include <Rcpp.h>
#include "commCpp.cpp"
#include "roleModelCpp.cpp"
#include "roleParamsCpp.cpp"

using namespace Rcpp;

//' @name specLocal
//' @title function to implement speciation for \code{localComm} class objects
//' @param x an object of class \code{localComm}
//' @param i the index to undergo speciation
//' @param p the roleModel params

localCommCpp specLocal(localCommCpp x, int i, roleParamsCpp p)
{
    // update Smax and initialize abundance
    // update number of species
    x.Smax += 1;

    // initialize abundance for new species
    x.abundance[x.Smax] = 1;

    // index of where unrealized traits begin
    int j = x.traitMax + 1;

    // add trait
    x.traits(j, 1) = x.Smax;
    x.traits(j, 2) = x.traits(i, 2) + R::rnorm(0, p.values.trait_sigma);

    //increase traitMax
    x.traitMax += 1;

    return(x);
}

roleModelCpp specRole(roleModelCpp x)
{

    // note: `Smax` from `@localComm` and `@metaComm` and `n` from `@phylo` are
    // all enforced to be equal, so we can sample from any but we have to
    // weight the probabilities by abundances and immigration

    // dispersal prob
    double dp = x.params.values.dispersal_prob;

    // normalized abundances at meta and local levels
    NumericVector mp = x.metaComm.abundance[Rcpp::Range(1,x.localComm.Smax)];
    mp = mp / sum(mp);
    NumericVector lp = x.localComm.abundance[Rcpp::Range(1,x.localComm.Smax)];
    lp = lp / sum(lp);

    // prob of selecting a parent
    // metacomm abundance weighted by dispersal prob + local comm abundance weighted by birth
    NumericVector pp = dp * mp + (1 - dp) * lp;

    // index of parent
    NumericVector i = sample(x.phylo.e, 1, false, pp);

    // update slots of the role model object
    x.localComm = specLocal(x.localComm, i[0], x.params);
    //x.phylo = speciation(x.phylo, i);

    return(x);
}

RCPP_MODULE(speciationCpp) {
    function("speciationL" , &specLocal);
    function("speciationR" , &specRole);
}

