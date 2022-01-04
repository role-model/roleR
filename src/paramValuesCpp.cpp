#pragma once

#include <RcppArmadillo.h>

using namespace Rcpp;

//' @name paramValuesCpp
//' @title an object containing the values of parameters, a field of roleParamsCpp
//' @field new Constructor initialize an empty paramValuesCpp containing defaults
//' @field speciation_meta the rate of speciation in the metacommunity
//' @field extinction_meta the rate of extinction in the metacommunity
//' @field trait_sigma the rate of trait change
//' @field species_meta the number of species in the metacommunity
//' @field individuals_meta the number of individuals in the metacommunity
//' @field individuals_local the number of individuals in the local community
//' @field dispersal_prob the probability of dispersal from the metacommunity to the local community
//' @field speciation_local the rate of speciation in the localcommunity
//' @field env_sigma
//' @field comp_sigma
//' @field mu
//' @field alpha
//'
class paramValuesCpp {
public:
    double speciation_meta;
    double extinction_meta;
    double trait_sigma;
    double species_meta;
    double individuals_meta;
    double individuals_local;
    double dispersal_prob;
    double speciation_local;
    double env_sigma;
    double comp_sigma;
    double mu;
    double alpha;
    double trait_z; // the optimal trait for an environment 
    double sigma_e; // determines how quickly fitness decays with the distance to the optimum
    double sigma_c; 

    //empty constructor
    paramValuesCpp(){
        //empty constructor assigns defaults
        species_meta = 15;
        individuals_meta = 500;
        individuals_local = 100;
        dispersal_prob = 0.5;
        speciation_local = 0.1;
        extinction_meta = 0.8;
        speciation_meta = 1;
        trait_sigma = 0.1;
        sigma_e = 0.1;
        sigma_c = 0.1; 
        trait_z = 0.5;
    }

    // //full constructor - may want to use this in the future
    // paramValuesCpp(double sm, double em, double ts, double spm,double im, double il, double dp, double sl, double es, double cs, double m, double a)
    //     : speciation_meta(sm), extinction_meta(em), trait_sigma(ts), species_meta(spm),
    //       individuals_meta(im), individuals_local(il), dispersal_prob(dp),
    //       speciation_local(sl), env_sigma(es), comp_sigma(cs), mu(m), alpha(a)
    // {
    // }
};

RCPP_EXPOSED_CLASS(paramValuesCpp)
