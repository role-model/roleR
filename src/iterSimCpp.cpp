#pragma once

#include <RcppArmadillo.h>
#include "roleModelCpp.cpp"

using namespace Rcpp;

//' @name iterSim
//' @title a function to iterate the RoLE model through \code{nstep} steps of simulation
//' @param an object of class \code{roleModelCpp} to iterate 
//' @param nstep an \code{integer} of the number of steps to iterate
//' @param print a \code{bool} specifying whether to print the step outcomes
//' 
// NOTE - doesn't appear to make a deep copy, keep this in mind 
roleModelCpp iterSim(roleModelCpp model,int nstep, bool print) {
  
    roleModelCpp out = model;
    for(int i = 0; i < nstep; i++) {

        // death always occurs
        if(print){Rcout << "start death event" << "\n";}
        int dead_index = out.death();
        if(print){Rcout << "end death event" << "\n";}
        
        //if speciation occurs based on local speciation param
        if(R::runif(0,1) <= model.params.values.speciation_local){
            if(print){Rcout << "start speciation event: chance was" <<
              model.params.values.speciation_local << "\n";}
            out.speciation(dead_index);
            if(print){Rcout << "end speciation event" << "\n";}
        }

        // else if immigration occurs based on dispersal chance param
        else if(R::runif(0,1) <= model.params.values.dispersal_prob) { //else
          if(print){Rcout << "start immigration event: chance was" <<
            model.params.values.dispersal_prob << "\n";}
            out.immigration(dead_index);
          if(print){Rcout << "end immigration event" << "\n";}
        }

        // else a birth event occurs
        else{
            if(print){Rcout << "start birth event" << "\n";}
            out.birth(dead_index);
            if(print){Rcout << "end birth event" << "\n";}
        }
    }
    return(out);
}

RCPP_MODULE(iterSimCpp) {
  function("iterSim", &iterSim);
}