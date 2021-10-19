#include <Rcpp.h>
#include <commCpp.cpp>
#include <roleModelCpp.cpp>
#include <roleParamsCpp.cpp>

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
    x.abundance[x.Smax] = 1

    // index of where unrealized traits begin
    j = min(which(is.na(x.traits[,1])))

    // add trait
    x.traits[j, 1] = x.Smax
    x.traits[j, 2] = x.traits[i, 2] + rnorm(1, 0, p.values.trait_sigma)

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
    mp = mp / R::sum(mp);
    lp = x.localComm.abundance[1:x.localComm.Smax];
    lp = lp / R::sum(lp);

    // prob of selecting a parent
    // metacomm abundance weighted by dispersal prob + local comm abundance weighted by birth
    pp = dp * mp + (1 - dp) * lp;

    // index of parent
    i = R::sample(x.phylo.n, 1, pp);

    // update slots of the role model object
    x.localComm = speciation(x.localComm, i, x.params);
    x.metaComm = speciation(x.metaComm);
    x@phylo = speciation(x.phylo, i);

    return(x);
}
