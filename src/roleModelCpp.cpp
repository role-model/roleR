#include <RcppArmadillo.h>
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

        // constructor
        roleModelCpp(localCommCpp local_, metaCommCpp meta_, rolePhyloCpp phy_,
                     roleParamsCpp params_) : localComm(local_), metaComm(meta_),
                     phylo(phy_), params(params_)
        {
        }
        
        // samples an individual and calls localComm.birth(individual), replacing the indv at dead_index
        void birth(int dead_index) 
        {
            // sample an individual for birth randomly
            NumericVector probs = localComm.abundance_indv;
            IntegerVector i = Rcpp::sample(localComm.J, 1, false, probs);
            
            // print contents of NumericVector
            //for(int i=0; i<probs.length(); i++){
            //    Rprintf("the value of v[%i] : %f \n", i, probs[i]);
            //}

            // make i from 0 to J - 1 (previously 1 to J)
            i[0] -= 1;

            // call birth on an individual
            //localComm.birth(i[0]);
            localComm.birth(i[0], dead_index);
        }

        int death()
        {
            // trait_z is the optimal trait for an environment 
            // sigma_e determines how quickly fitness decays with the distance to the optimum
            // NOTE - allow trait_z to vary as a time series later
            
            // probs of death due to environmental filtering
            NumericVector f_probs = 1 - exp(-1/params.values.sigma_e * pow(localComm.traits - params.values.trait_z, 2));
            //calculate vector once and each time individual is replaced then update index 
              
            // init c_probs
            //NumericVector c_probs = NumericVector(localComm.J);
            
            // WIP for each individual
            // does not use traitdiffs matrix currently
            // for(int i = 0; i < localComm.Imax; i++)
            // {
            //   float sum = 0; 
            //   
            //   for(int j = 0; j < localComm.Imax; j++)
            //   {
            //     sum += exp(-1/params.values.sigma_e * pow(localComm.traits[i] - localComm.traits[j], 2));
            //   }
            //   
            //   c_probs[i] = (1/localComm.Imax) * sum; 
            // }

            // calculate probs of death due to competitive filtering
            NumericVector c_probs = as<NumericVector>(wrap(1/localComm.J * arma::sum(exp((-1/params.values.sigma_c) * arma::pow(localComm.traitdiffs, 2)),0)));
            
            
            Rcout << "fprobs size: " << f_probs.size() << "\n";
            
            // print probs
            for(int i=0; i<f_probs.length(); i++){
              Rprintf("the value of v[%i] : %f \n", i, f_probs[i]);
            }
            
            Rcout << "cprobs size: " << c_probs.size() << "\n";
            
            // print probs
            for(int i=0; i<c_probs.length(); i++){
              Rprintf("the value of v[%i] : %f \n", i, c_probs[i]);
            }
            
            // probs is sum of f_probs and c_probs 
            NumericVector probs = f_probs + c_probs; 
            
            Rcout << "probs size: " << probs.size() << "\n";
            Rcout << "J size : " << localComm.J << "\n";
            
            // print probs
            for(int i=0; i<probs.length(); i++){
               Rprintf("the value of v[%i] : %f \n", i, probs[i]);
            }
            
            IntegerVector i = Rcpp::sample(localComm.J, 1, false, probs);
            
            // make i from 0 to J - 1 (previously 1 to J)
            i[0] -= 1;
            
            // call death on individual
            localComm.death(i[0]);

            // if death of indv led to extinction of species, call death on rolePhylo
            if(localComm.abundance_sp[localComm.species_ids[i[0]]] <= 0)
            {
                phylo.death(i[0]);
            }
            
            // return the index thqt death was called on 
            return i[0];
        }

        void speciation(int dead_index)
        {
            // note: `Smax` from `@localComm` and `@metaComm` and `n` from `@phylo` are
            // all enforced to be equal, so we can sample from any but we have to
            // weight the probabilities by abundances and immigration

            // dispersal prob
            double dp = params.values.dispersal_prob;

            // normalized abundances at meta and local levels
            // mp has length equal to the number of species in the metacomm
            NumericVector mp = metaComm.abundance[Rcpp::Range(0,localComm.Smax-1)]; //Smax - 1
            mp = mp / sum(mp);
            
            // prob of selecting a parent for speciation depends on abundance 
            // metacomm abundance weighted by dispersal prob + local comm abundance weighted by birth
            
            //CHANGE - how to do this? decollapse metacomm abundance? cant really do that 
            // collapse individual vector?
            // lp has length equal to the number of species in the localcom
            // EZ - save species abundances, should be easy enough 
            NumericVector sa = NumericVector(localComm.Smax); 
   
            //for(int s = 0; s < localComm.Smax; s++)
            //{
            // NumericVector indices = match(localComm.species_ids, s);
            
            //  sa[s] = match(localComm.species_ids
            //}
            
            NumericVector lp = localComm.abundance_sp[Rcpp::Range(0,localComm.Smax-1)]; //Smax - 1
            lp = lp / sum(lp);
            
            NumericVector pp = dp * mp + (1 - dp) * lp;
            
            // vector of phylo parents
            // this is not the case because includes extinct edges 
            //NumericVector v = phylo.e(_, 0); 
            
            //Rcout << "vect: " << v << "\n";
            Rcout << "probs size: " << pp.size() << "\n";
            Rcout << "phylo n size : " << phylo.n << "\n";
            
            // index of parent
            IntegerVector i = sample(phylo.n, 1, false, pp);
            
            // make i from 0 to phylo.n - 1 (previously 1 to phylo.n)
            i[0] -= 1;
            
            // call speciation on individual index to replace dead indv with
            // first member of new species
            localComm.speciation(i[0], dead_index, params);
            
            // add a tip to the phylogeny
            phylo.speciation(i[0]);
        }

        void immigration(int dead_index)
        {
            // Smax was not changing - I think I fixed this 
            Rcout << "Smax: " << metaComm.Smax << "\n";
        
            //sample a species for birth relative to local abundance
            //0 vs 1 start indices may cause problems
            
            NumericVector probs = metaComm.abundance[Rcpp::Range(0,params.values.species_meta-1)];//metaComm.Smax-1
          
            Rcout << "probs size: " << probs.size() << "\n";
            Rcout << "Smax n size : " << metaComm.Smax << "\n";
            
            IntegerVector i = sample(params.values.species_meta, 1, false, probs);
            
            // make i from 0 to Smax - 1 (previously 1 to Smax)
            i[0] -= 1;
            
            //Rcout << "i: " << i[0] << "\n";
            
            // call immigration on species i
            localComm.immigration(i[0], dead_index, metaComm);
        }
};

RCPP_EXPOSED_CLASS(roleModelCpp)
