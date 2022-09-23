//#pragma once

#include <string>
#include <Rcpp.h>
using namespace Rcpp;

class roleParamsCpp {
public:
    NumericVector individuals_local;
    //NumericVector individuals_meta;
    //NumericVector species_meta;
    NumericVector speciation_local;
    NumericVector speciation_meta;
    NumericVector extinction_meta;
    NumericVector trait_sigma; 
    NumericVector env_sigma;
    NumericVector comp_sigma; 
    NumericVector neut_delta; 
    NumericVector dispersal_prob; 
    //NumericVector mutation_rate; 
    NumericVector equilib_escape;
    //NumericVector num_basepairs;
    
   roleParamsCpp(RObject params, int niter){
        individuals_local = Rcpp::as<NumericVector>(params.slot("individuals_local"));
        speciation_local = Rcpp::as<NumericVector>(params.slot("speciation_local"));
        speciation_meta = Rcpp::as<NumericVector>(params.slot("speciation_meta"));
        extinction_meta = Rcpp::as<NumericVector>(params.slot("extinction_meta"));
        trait_sigma = Rcpp::as<NumericVector>(params.slot("trait_sigma"));
        env_sigma = Rcpp::as<NumericVector>(params.slot("env_sigma"));
        comp_sigma = Rcpp::as<NumericVector>(params.slot("comp_sigma"));
        neut_delta = Rcpp::as<NumericVector>(params.slot("neut_delta"));
        dispersal_prob = Rcpp::as<NumericVector>(params.slot("dispersal_prob"));
        
        if(individuals_local.length() == 1){
            individuals_local = rep(individuals_local, niter);
        }
        
        if(speciation_local.length() == 1){
            speciation_local = rep(speciation_local, niter);
        }
        
        if(speciation_meta.length() == 1){
            speciation_meta = rep(speciation_meta, niter);
        }
        if(extinction_meta.length() == 1){
            extinction_meta = rep(extinction_meta, niter);
        }
        if(trait_sigma.length() == 1){
            trait_sigma = rep(trait_sigma, niter);
        }
        if(env_sigma.length() == 1){
            env_sigma = rep(env_sigma, niter);
        }
        if(comp_sigma.length() == 1){
            comp_sigma = rep(comp_sigma, niter);
        }
        if(neut_delta.length() == 1){
            neut_delta = rep(neut_delta, niter);
        }
        if(dispersal_prob.length() == 1){
            dispersal_prob = rep(dispersal_prob, niter);
        }
    }
};

RCPP_EXPOSED_CLASS(roleParamsCpp)