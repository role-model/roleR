#include <Rcpp.h>
#pragma once

using namespace Rcpp;

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
    }

    // //full constructor
    // paramValuesCpp(double sm, double em, double ts, double spm,double im, double il, double dp, double sl, double es, double cs, double m, double a)
    //     : speciation_meta(sm), extinction_meta(em), trait_sigma(ts), species_meta(spm),
    //       individuals_meta(im), individuals_local(il), dispersal_prob(dp),
    //       speciation_local(sl), env_sigma(es), comp_sigma(cs), mu(m), alpha(a)
    // {
    // }
};

RCPP_EXPOSED_CLASS(paramValuesCpp)

RCPP_MODULE(paramValsCpp) {
    class_<paramValuesCpp>("paramValuesCpp")
        .constructor()
        .field("speciation_meta", &paramValuesCpp::speciation_meta)
        .field("extinction_meta", &paramValuesCpp::extinction_meta)
        .field("trait_sigma", &paramValuesCpp::extinction_meta)
        .field("species_meta", &paramValuesCpp::species_meta)
        .field("individuals_meta", &paramValuesCpp::individuals_meta)
        .field("individuals_local", &paramValuesCpp::individuals_local)
        .field("dispersal_prob", &paramValuesCpp::dispersal_prob)
        .field("speciation_local", &paramValuesCpp::speciation_local)
        .field("dispersal_prob", &paramValuesCpp::dispersal_prob)
        //.constructor<double, double, double, double,double, double, double, double, double, double, double, double>()
        ;
}
