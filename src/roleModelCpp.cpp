#include <Rcpp.h>
#include "commCpp.cpp"
#include "rolePhyloCpp.cpp"
#include "roleParamsCpp.cpp"
#pragma once

using namespace Rcpp;

//' @name roleModelCpp
//' @title a C++ class to specify the entire RoLE model
//' @field new Constructor
//' @field localComm the local community, an object of class localCommCpp
//' @field metaComm the metacommunity, an object of class metaCommCpp
//' @field phylo the phylogeny of the metacommunity, an object of class rolePhyloCpp
//' @field params the parameters for the model, an object of class roleParamsCpp
//'

class roleModelCpp {
    public:
        localCommCpp localComm;
        metaCommCpp metaComm;
        rolePhyloCpp phylo;
        roleParamsCpp params;

        //constructor
        roleModelCpp(localCommCpp local_, metaCommCpp meta_, rolePhyloCpp phy_,
                     roleParamsCpp params_) : localComm(local_), metaComm(meta_),
                     phylo(phy_), params(params_)
        {
        }

        void birth()
        {
            //sample a species for birth relative to local abundance
            //0 vs 1 start indices may cause problems
            NumericVector probs = localComm.abundance[Rcpp::Range(0,localComm.Smax-1)];

            //for(int i=0; i<probs.length(); ++i){
            //    Rprintf("the value of v[%i] : %f \n", i, probs[i]);
            //}

            IntegerVector i = sample(localComm.Smax, 1, false, probs);

            // make i from 0 to Smax - 1 (previously 1 to Smax)
            i[0] -= 1;

            // call birth
            localComm.birth(i[0]);
        }

        void death()
        {
            // sample a species for death proportional to species abundance
            NumericVector probs = localComm.abundance[Rcpp::Range(0,localComm.Smax-1)]; //localComm.Smax
            IntegerVector i = sample(localComm.Smax, 1, false, probs);

            localComm.death(i[0]);

            // if death led to extinction, call death on rolePhylo
            if(localComm.abundance[i[0]] <= 0)
            {
                phylo.death(i[0]);
            }
        }

        void speciation()
        {
            // note: `Smax` from `@localComm` and `@metaComm` and `n` from `@phylo` are
            // all enforced to be equal, so we can sample from any but we have to
            // weight the probabilities by abundances and immigration

            // dispersal prob
            double dp = params.values.dispersal_prob;

            // normalized abundances at meta and local levels
            // mp has length equal to the number of species in the metacomm
            NumericVector mp = metaComm.abundance[Rcpp::Range(0,localComm.Smax)]; //Smax - 1
            mp = mp / sum(mp);
            // lp has length equal to the number of species in the localcom
            NumericVector lp = localComm.abundance[Rcpp::Range(0,localComm.Smax)]; //Smax - 1
            lp = lp / sum(lp);
            
            // prob of selecting a parent for speciation depends of abundance 
            // metacomm abundance weighted by dispersal prob + local comm abundance weighted by birth
            NumericVector pp = dp * mp + (1 - dp) * lp;
            
            // vector of phylo parents
            // this is not the case because includes extinct edges 
            //NumericVector v = phylo.e(_, 0); 
            
            //Rcout << "vect: " << v << "\n";
            Rcout << "probs size: " << pp.size() << "\n";
            Rcout << "phylo n size : " << phylo.n << "\n";
            
            // index of parent
            IntegerVector i = sample(phylo.n, 1, false, pp);
            
            // update slots of the role model object
            localComm.speciation(i[0], params);
            phylo.speciation(i[0]);
        }

        void immigration()
        {
            //for(int i=0; i<probs.length(); ++i){
            //    Rprintf("the value of v[%i] : %f \n", i, probs[i]);
            //}

            //sample a species for birth relative to local abundance
            //0 vs 1 start indices may cause problems
            NumericVector probs = metaComm.abundance[Rcpp::Range(0,metaComm.Smax)];//metaComm.Smax-1
            IntegerVector i = sample(metaComm.Smax, 1, false, probs);

            localComm.immigration(i[0]);
        }
};

RCPP_EXPOSED_CLASS(roleModelCpp)
