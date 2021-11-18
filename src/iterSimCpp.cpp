#pragma once

#include <Rcpp.h>
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
        out.death();

        //if speciation occurs based on local speciation param
        if(R::runif(1,1) <= model.params.values.speciation_local){
            //out.speciation();
        }

        // else if immigration occurs based on dispersal chance param
        if(R::runif(1,1) <= model.params.values.dispersal_prob) { //else
            out.immigration();
        }

        // else a birth event occurs
        else{
            out.birth();
        }
    }
    return(out);
}

RCPP_MODULE(iterSimCpp) {
    function("iterSim" , &iterSim);
}
