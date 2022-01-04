#pragma once

#include <RcppArmadillo.h>
#include "roleModelCpp.cpp"
//#include "birthCpp.cpp"
//#include "deathCpp.cpp"
//#include "speciationCpp.cpp"
//#include "immigrationCpp.cpp"

using namespace Rcpp;

//' @name iterSimCpp
//' @title a C++ function to iterate the RoLE model through nstep steps of simulation
//' @param model the entire C++ model, an object of class roleModelCpp
//' @param nstep the number of steps to iterate

roleModelCpp iterSim(roleModelCpp model,int nstep) {

    roleModelCpp out = model;
    for(int i = 0; i < nstep; i++) {

        // death always occurs
        Rcout << "death triggered" << "\n";
        int dead_index = out.death();
        Rcout << "death has occurred" << "\n";
        
        //if speciation occurs based on local speciation param
        if(R::runif(0,1) <= model.params.values.speciation_local){
            Rcout << "speciation triggered" << "\n";
            out.speciation(dead_index);
            Rcout << "speciation occurred" << "\n";
        }

        // else if immigration occurs based on dispersal chance param
        if(R::runif(0,1) <= model.params.values.dispersal_prob) { //else
            Rcout << "immigration triggered" << "\n";
            out.immigration(dead_index);
            Rcout << "immigration occurred" << "\n";
        }

        // else a birth event occurs
        else{
            Rcout << "birth triggered" << "\n";
            out.birth(dead_index);
            Rcout << "birth has occurred" << "\n";
        }
    }
    return(out);
}

RCPP_MODULE(iterSimCpp) {
  function("iterSim", &iterSim);
}