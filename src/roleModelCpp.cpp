#include <RcppArmadillo.h>
#include "commCpp.cpp"
#include "rolePhyloCpp.cpp"
#include "roleDataCpp.cpp" 
#include <string>
#include <array>
#pragma once

using namespace Rcpp;

//' @name roleModelCpp
//' @title a class to specify the entire RoLE model
//' @field new Constructor
//' @field local the local community, an object of class localCommCpp
//' @field meta the metacommunity, an object of class metaCommCpp
//' @field phylo the phylogeny of the metacommunity, an object of class rolePhyloCpp
//' @field params the parameters for the model, an object of class roleParamsCpp
//' @field timeseries a list of past simulation states 

class roleModelCpp {
    public:
        localCommCpp local;
        metaCommCpp meta;
        rolePhyloCpp phylo;
        List params; 
        //std::vector<roleDataCpp> timeseries;
        List timeseries;
        // std::list<roleDataCpp> timeseries;
        //std::array<roleDataCpp,10> timeseries;
        NumericVector stats; 
        int iter;  
        //int niter_save; 
        bool print; 
        bool print_vectors; 
        //roleDataCpp timeseries[];
        
        // constructor
        roleModelCpp(localCommCpp local_, metaCommCpp meta_, rolePhyloCpp phy_, List params_) : local(local_), meta(meta_),
                     phylo(phy_), params(params_)
        {
          timeseries = List(10);
          //timeseries = std::vector<roleDataCpp>(10);
          //timeseries = roleDataCpp[10];
          //timeseries = std::array<roleDataCpp,10>(); 
        }
        
        // samples an individual and calls localComm.birth(individual), replacing the indv at dead_index
        void birth(int dead_index) 
        {
            // sample an individual for birth randomly
            NumericVector probs = local.abundance_indv;
            IntegerVector i = Rcpp::sample(local.J, 1, false, probs);

            // make i from 0 to J - 1 (previously 1 to J)
            i[0] -= 1;

            // call birth on an individual
            local.birth(i[0], dead_index);
        }

        int death()
        {
            // NOTE - calculate vector once and each time individual is replaced then update index 
            
            // get params from params List
            NumericVector env_sigma = params["env_sigma"];
            NumericVector trait_z = params["trait_z"];
            
            if(print){Rcout << "calculating probs of death due to env" << "\n";}

            // compute probs of death due to environmental filtering
            NumericVector f_probs = 1 - exp(-1/env_sigma[iter] * pow(local.traits - trait_z[iter], 2));
            
            // get params from params List
            NumericVector comp_sigma = params["comp_sigma"];
            
            if(print){Rcout << "calculating probs of death due to comp" << "\n";}
            
            // compute probs of death due to competitive filtering
            NumericVector c_probs = as<NumericVector>(wrap(1/local.J * arma::sum(exp((-1/comp_sigma(iter)) * arma::pow(local.traitdiffs, 2)),0)));
            
            // prints size of vector, contents, and local J value to compare 
            if(print){printVector(f_probs, "f_probs");}
            if(print){printVector(c_probs, "c_probs");}
            
            // probs is sum of f_probs and c_probs 
            NumericVector probs = f_probs + c_probs; 
            
            if(print){printVector(probs, "sum probs");}
            
            // remove negative probabilities
            probs = (abs(probs)+probs)/2;
            
            IntegerVector i = Rcpp::sample(local.J, 1, false, probs);
            
            // make i from 0 to J - 1 (previously 1 to J)
            i(0) -= 1;
            
            if(print){Rcout << "chosen index : " << i(0) << "\n";}
            
            // call death on individual
            local.death(i(0));

            // if death of indv led to extinction of species, call death on rolePhylo
            if(local.abundance_sp[local.species_ids[i(0)]] <= 0)
            {
                if(print){Rcout << "extinction occured, calling phylodeath" << "\n";}
                
                phylo.death(i(0));
            }
            
            // return the index thqt death was called on 
            return i(0);
        }

        void speciation(int dead_index)
        {
            // NOTE - `Smax` from `@local` and `@metaComm` and `n` from `@phylo` are
            // all enforced to be equal, so we can sample from any but we have to
            // weight the probabilities by abundances and immigration
          
            // get param from params List
            NumericVector dispersal_prob = params("dispersal_prob");
            double dp = dispersal_prob(iter);
              
            if(print){Rcout << "computing probs of speciation as sum of weighted meta and local abundances" << "\n";}
              
            // prob of selecting a parent for speciation is metacomm abundance weighted by dispersal prob + local comm abundance weighted by birth
            // NOTE -  can come up with negative probs now, is this intended behavior? 
            //NumericVector lp = local.abundance_sp; //[Rcpp::Range(0,local.Smax - 1)]; //local.Smax-1
            //lp = lp / sum(lp);
            
            NumericVector probs = dp * (meta.abundance / sum(meta.abundance)) + 
                                  (1 - dp) * (local.abundance_sp / sum(local.abundance_sp));
            
            // remove negative probabilities
            probs = (abs(probs)+probs)/2;
            
            if(print){printVector(probs, "speciation probs");}
            //if(print){Rcout << "phylo n" << phylo.n << "\n";}
            
            // index of parent
            IntegerVector i = sample(meta.abundance.length(), 1, false, probs);
        
            // make i from 0 to phylo.n - 1 (previously 1 to phylo.n)
            i(0) -= 1;
            
            if(print){Rprintf("chosen index : ", i(0));}
            
            // get param from params List to feed into localComm speciation
            NumericVector trait_sigma = params("trait_sigma");
            
            // call speciation on individual index to replace dead indv with first member of new species
            if(print){Rcout << "calling local speciation" << "\n";}
            local.speciation(i(0), dead_index, trait_sigma(iter));
            
            // add a tip to the phylogeny
            if(print){Rcout << "calling phylo speciation" << "\n";}
            // phylo.speciation(i(0));
        }

        void immigration(int dead_index)
        {
            if(print){Rcout << "metacomm Smax size: " << meta.abundance.length() << "\n";}
            // Smax was not changing - I think I fixed this 
            
            //sample a species for birth relative to local abundance
            NumericVector probs = meta.abundance[Rcpp::Range(0,meta.abundance.length() - 1)];
          
            if(print){Rcout << "imm from meta probs size: " << probs.size() << "\n";}
            
            IntegerVector i = sample(meta.abundance.length(), 1, false, probs); //metaComm.abundance.length() - 1
            
            // make i from 0 to Smax - 1 (previously 1 to Smax)
            i(0) -= 1;
            
            if(print){Rcout << "chosen index : " << i(0) << "\n";}
            
            // if species i was dead in the phylo, make it alive again
            phylo.alive(i(0)) = true; 
            
            // call immigration on species i
            local.immigration(i(0), dead_index, meta);
        }
        
        //NOTE - unused for now, we are doing this in R 
        // save summary statistics
        void computeStatistics()
        {
          
          double shannon_entropy = 0;
          NumericVector abundance_sp_relative = local.abundance_sp / max(local.abundance_sp);
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
        roleDataCpp copyData()
        {
          if(print){Rcout << "copying local" << "\n";}
          localCommCpp l = localCommCpp(clone(local.abundance_sp),clone(local.traits_sp),clone(local.pi_sp),
                                        clone(local.abundance_indv),clone(local.traits));
          if(print){Rcout << "copying meta" << "\n";}
          metaCommCpp m = metaCommCpp(clone(meta.abundance),clone(meta.traits));
          if(print){Rcout << "copying phylo" << "\n";}
          rolePhyloCpp ph = rolePhyloCpp(int(phylo.n),clone(phylo.e),clone(phylo.l),
                                         clone(phylo.alive),clone(phylo.tipNames),double(phylo.scale));
          // compute stats
          if(print){Rcout << "creating data" << "\n";}
          roleDataCpp out = roleDataCpp(l,m,ph);
          out.iter_num = iter;
          return out; 
        }
        
        void addTimeseriesStep(int i)
        {
          if(print){Rcout << "adding timeseries step" << "\n";}
          timeseries.push_back(copyData());
          //timeseries[i] = copyData();
        }
        
        void printTimeseries()
        {
          roleDataCpp temp = timeseries[0];
          Rcout << temp.local.abundance_sp << "\n";
          temp = timeseries[1];
          Rcout << temp.local.abundance_sp << "\n";
          temp = timeseries[2];
          Rcout << temp.local.abundance_sp << "\n";
          temp = timeseries[9];
          Rcout << temp.local.abundance_sp << "\n";
        }
        //NOTE - unused for now, we are doing this in R 
        // return a timeseries
        void getTimeseries(NumericVector v, std::string name)
        {
        }
        
        void printVector(NumericVector v, std::string name)
        {
          Rcout << name << " size: " << v.size() << "\n";
          Rcout << "local J size : " << local.J << "\n";
          if(print_vectors)
          {
            for(int i=0; i<v.length(); i++){
                if(true){
                Rprintf("v[%i] : %f ", i, v[i],",");}
            }
          }
        }
};

RCPP_EXPOSED_CLASS(roleModelCpp)
