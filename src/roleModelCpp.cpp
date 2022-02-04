#include <RcppArmadillo.h>
#include "commCpp.cpp"
#include "rolePhyloCpp.cpp"
#include "roleParamsCpp.cpp"
#include "roleDataCpp.cpp" 
#include <string>
#pragma once

using namespace Rcpp;

//' @name roleModelCpp
//' @title a class to specify the entire RoLE model
//' @field new Constructor
//' @field localComm the local community, an object of class localCommCpp
//' @field metaComm the metacommunity, an object of class metaCommCpp
//' @field phylo the phylogeny of the metacommunity, an object of class rolePhyloCpp
//' @field params the parameters for the model, an object of class roleParamsCpp
//' @field timeseries a list of past simulation states 

class roleModelCpp {
    public:
        localCommCpp localComm;
        metaCommCpp metaComm;
        rolePhyloCpp phylo;
        List params; 
        //std::map<std::string, NumericVector> params; 
        std::list<roleDataCpp> timeseries;
        NumericVector stats; 
        int iter;  
        //int niter_save; 
        bool print; 
        
        // constructor
        roleModelCpp(localCommCpp local_, metaCommCpp meta_, rolePhyloCpp phy_, List params_) : localComm(local_), metaComm(meta_),
                     phylo(phy_), params(params_)
        {
        }
        
        // samples an individual and calls localComm.birth(individual), replacing the indv at dead_index
        void birth(int dead_index) 
        {
            // sample an individual for birth randomly
            NumericVector probs = localComm.abundance_indv;
            IntegerVector i = Rcpp::sample(localComm.J, 1, false, probs);

            // make i from 0 to J - 1 (previously 1 to J)
            i[0] -= 1;

            // call birth on an individual
            localComm.birth(i[0], dead_index);
        }

        int death()
        {
            // NOTE - calculate vector once and each time individual is replaced then update index 
            
            // get params from params List
            NumericVector env_sigma = params["env_sigma"];
            NumericVector trait_z = params["trait_z"];
            
            // compute probs of death due to environmental filtering
            NumericVector f_probs = 1 - exp(-1/env_sigma[iter] * pow(localComm.traits - trait_z[iter], 2));

            // WIP alternate version for each individual without traitdiffs 
            // init c_probs
            //NumericVector c_probs = NumericVector(localComm.J);
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
            
            // get params from params List
            NumericVector comp_sigma = params["comp_sigma"];
            
            // compute probs of death due to competitive filtering
            NumericVector c_probs = as<NumericVector>(wrap(1/localComm.J * 
            arma::sum(exp((-1/comp_sigma[iter]) * arma::pow(localComm.traitdiffs, 2)),0)));
            
            // prints size of vector, contents, and local J value to compare 
            if(print){printVector(f_probs, "f_probs");}
            if(print){printVector(c_probs, "c_probs");}
            
            // probs is sum of f_probs and c_probs 
            NumericVector probs = f_probs + c_probs; 
            
            if(print){printVector(probs, "sum probs");}
            
            IntegerVector i = Rcpp::sample(localComm.J, 1, false, probs);
            
            // make i from 0 to J - 1 (previously 1 to J)
            i[0] -= 1;
            
            if(print){Rprintf("chosen index : ", i[0]);}
            
            // call death on individual
            localComm.death(i[0]);

            // if death of indv led to extinction of species, call death on rolePhylo
            if(localComm.abundance_sp[localComm.species_ids[i[0]]] <= 0)
            {
                if(print){Rprintf("extinction occured, calling phylodeath");}
                phylo.death(i[0]);
            }
            
            // return the index thqt death was called on 
            return i[0];
        }

        void speciation(int dead_index)
        {
            // NOTE - `Smax` from `@localComm` and `@metaComm` and `n` from `@phylo` are
            // all enforced to be equal, so we can sample from any but we have to
            // weight the probabilities by abundances and immigration
          
            // get param from params List
            NumericVector dispersal_prob = params["dispersal_prob"];
            double dp = dispersal_prob[iter];

            // compute normalized abundances at meta and local levels - mp has length equal to the number of species in the metacomm
            NumericVector mp = metaComm.abundance[Rcpp::Range(0,localComm.Smax-1)];
            mp = mp / sum(mp);
            
            //NumericVector sa = NumericVector(localComm.Smax); 
            //for(int s = 0; s < localComm.Smax; s++)
            //{
            // NumericVector indices = match(localComm.species_ids, s);
            
            //  sa[s] = match(localComm.species_ids
            //}
          
            // prob of selecting a parent for speciation is metacomm abundance weighted by dispersal prob + local comm abundance weighted by birth
            // NOTE -  can come up with negative probs now, is this intended behavior? 
            NumericVector lp = localComm.abundance_sp[Rcpp::Range(0,localComm.Smax-1)]; 
            lp = lp / sum(lp);
            
            NumericVector probs = dp * mp + (1 - dp) * lp;
            
            // remove negative probabilities
            probs = (abs(probs)+probs)/2;
            
            if(print){printVector(probs, "speciation probs");}
            
            // index of parent
            IntegerVector i = sample(phylo.n, 1, false, probs);
        
            // make i from 0 to phylo.n - 1 (previously 1 to phylo.n)
            i[0] -= 1;
            
            if(print){Rprintf("chosen index : ", i[0]);}
            
            // get param from params List to feed into localComm speciation
            NumericVector trait_sigma = params["trait_sigma"];
            
            // call speciation on individual index to replace dead indv with first member of new species
            if(print){Rprintf("calling localComm speciation");}
            localComm.speciation(i[0], dead_index, trait_sigma[iter]);
            
            // add a tip to the phylogeny
            if(print){Rprintf("calling phylo speciation");}
            phylo.speciation(i[0]);
        }

        void immigration(int dead_index)
        {
            if(print){Rcout << "metacomm Smax size: " << metaComm.Smax << "\n";}
            // Smax was not changing - I think I fixed this 
            
            //sample a species for birth relative to local abundance
            NumericVector probs = metaComm.abundance[Rcpp::Range(0,metaComm.abundance.length() - 1)];
          
            if(print){Rcout << "imm from meta probs size: " << probs.size() << "\n";}
            
            IntegerVector i = sample(metaComm.abundance.length() - 1, 1, false, probs);
            
            // make i from 0 to Smax - 1 (previously 1 to Smax)
            i[0] -= 1;
            
            if(print){Rcout << "chosen index : " << i[0] << "\n";}
            
            // call immigration on species i
            localComm.immigration(i[0], dead_index, metaComm);
        }
        
        //NOTE - unused for now, we are doing this in R 
        // save summary statistics
        void computeStatistics()
        {
          
          double shannon_entropy = 0;
          NumericVector abundance_sp_relative = localComm.abundance_sp / max(localComm.abundance_sp);
          for(int i = 0; i < abundance_sp_relative.length(); i++)
          {
            shannon_entropy += abundance_sp_relative[i] * log(abundance_sp_relative[i]);
          }
          shannon_entropy = -1 * shannon_entropy; 
          
          //get only branch lengths of alive in meta
          double faith_pd = sum(phylo.l);
          stats = NumericVector::create(Named("faith_pd",faith_pd), Named("y")=2 , _["z"]=3);
        }
        
        // create a deep copy of this object's data for timeseries
        roleDataCpp copyData(int i)
        {
          localCommCpp l = localCommCpp(localComm.abundance_sp,localComm.traits_sp,localComm.Smax,localComm.pi_sp);
          metaCommCpp m = metaCommCpp(metaComm.abundance,metaComm.traits,metaComm.Smax);
          rolePhyloCpp ph = rolePhyloCpp(phylo.n,phylo.e,phylo.l,phylo.alive,phylo.tipNames,phylo.scale);
          // compute stats
          roleDataCpp out = roleDataCpp(l,m,ph,i);
          return out; 
        }
        
        //NOTE - unused for now, we are doing this in R 
        // return a timeseries
        void getTimeseries(NumericVector v, std::string name)
        {
          Rcout << name << " size: " << v.size() << "\n";
          Rcout << "local J size : " << localComm.J << "\n";
          
          for(int i=0; i<v.length(); i++){
            Rprintf("v[%i] : %f ", i, v[i],",");
          }
        }
        
        void printVector(NumericVector v, std::string name)
        {
          Rcout << name << " size: " << v.size() << "\n";
          Rcout << "local J size : " << localComm.J << "\n";
          
          for(int i=0; i<v.length(); i++){
              Rprintf("v[%i] : %f ", i, v[i],",");
          }
        }
};

RCPP_EXPOSED_CLASS(roleModelCpp)
