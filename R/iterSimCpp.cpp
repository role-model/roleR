#include <Rcpp.h>
#include <roleModelCpp.cpp>
#include <birthCpp.cpp>
#include <deathCpp.cpp>
//#include <speciationCpp.cpp>
#include <immigrationCpp.cpp>

using namespace Rcpp;

// [[Rcpp::export]]

//' @name iterSimCpp
//' @title a C++ function to iterate the RoLE model
//' @param model the entire C++ model, an object of class roleModelCpp
//' @param nstep the number of steps to iterate

roleModelCpp iterSim(roleModelCpp model,int nstep) {

    roleModelCpp out = model;
    for(int i = 0; i < nstep; i++) {

        // death always occurs
        roleModelCpp out = deathRole(model);

        // if speciation occurs based on local speciation param
        //if(runif(1) <= model.params.values.speciation_local){
        //    model = speciationR(model)
        //}

        // else if immigration occurs based on dispersal chance param
        if(R::runif(1,1) <= model.params.values.dispersal_prob) { //else
            out = immigrationRole(out);
        }

        // else a birth event occurs
        else{
            out = birthRole(out);
        }
    }
    return(out);
}


RCPP_MODULE(iterSim) {
    function("iterSim" , &iterSim);
}
