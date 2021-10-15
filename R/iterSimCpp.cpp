#include <Rcpp.h>
#include <roleModelCpp.cpp>

using namespace Rcpp;

// [[Rcpp::export]]

void iterSim(roleModelCpp model,int nstep) {
    for(int i = 0; i < nstep; i++) {

        //death
        model.death();

        if(runif(1) <= model.params.values.speciation_local){
            model = speciation(model)
        }

        // if immigration occurs based on dispersal chance param
        else if(runif(1) <= model.params.values.dispersal_prob) {

            model = immigration(model)
        }

        else{
            model = birth(model)
        }
    }
}


RCPP_MODULE(iterSim) {
    function("iterSim" , &iterSim);
}
