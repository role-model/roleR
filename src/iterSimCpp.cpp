#pragma once

#include <RcppArmadillo.h>
#include "roleModelCpp.cpp"

using namespace Rcpp;

//' @name iterSim
//' @title a function to iterate the RoLE model through \code{nstep} steps of simulation
//' @param an object of class \code{roleModelCpp} to iterate 
//' @param nstep an \code{integer} of the number of steps to iterate
//' @param print a \code{bool} specifying whether to print the step outcomes

roleModelCpp iterSim(roleModelCpp model,int niter, int niter_timestep, bool print) {

    if(print){Rcout << "iter loop started" << "\n";}
    
    // must do a check somewhere to determine if niter is a multiple of niter_timestep
    //model.timeseries = List(niter/niter_timestep);
    
    //model.timeseries = roleDataCpp[int(niter/niter_timestep)];
    //model.timeseries = std::vector<roleDataCpp>(niter/niter_timestep);
    int timeseries_index = 0;
    
    //roleModelCpp model = model.copy(); NOTE - deprecated iterSim returning a roleModel, now acts on specified roleModel 
    for(int i = 0; i < niter - 1; i++) {
      
        // set the iter of the model
        model.iter = i; 

        // if reached steps per save
        if((i % niter_timestep) == 0)
        {
          if(print){Rcout << "copying timestep data " << timeseries_index <<"\n";}
          // add a deep copy of the simulation as it stands before the next step to the time series
          //if(print){Rcout << "copying timestep data" << "\n";}
          //model.timeseries[timeseries_index] = model.copyData();
          //roleDataCpp temp = model.timeseries[timeseries_index];
          //Rcout << "copied abundance sp " << temp.local.abundance_sp <<"\n";
          model.addTimeseriesStep(timeseries_index);
          timeseries_index++;
        }
        
        // get params
        NumericVector dispersal_prob = model.params["dispersal_prob"];
        NumericVector speciation_local = model.params["speciation_local"];
        
        // death always occurs
        if(print){Rcout << "start death event" << "\n";}
        int dead_index = model.death();
        if(print){Rcout << "end death event" << "\n";}
        
        //if speciation occurs based on local speciation param
        Rcout << "i" << i << "\n";
        double ru = R::runif(0,1);
        if(ru <= speciation_local[i]){
            //if(print){Rcout << "start speciation event: chance was" <<
            //model.params.values.speciation_local << "\n";}
            model.speciation(dead_index);
            if(print){Rcout << "end speciation event" << "\n";}
        }
        // else if immigration occurs based on dispersal chance param
        else if(R::runif(0,1) <= dispersal_prob[i]) { //else
          if(print){Rcout << "start immigration event: chance was " <<
            dispersal_prob[i] << "\n";}
            model.immigration(dead_index);
          if(print){Rcout << "end immigration event" << "\n";}
        }

        // else a birth event occurs
        else{
            if(print){Rcout << "start birth event" << "\n";}
              model.birth(dead_index);
            if(print){Rcout << "end birth event" << "\n";}
        }
    }
    
    return(model);
}

RCPP_MODULE(iterSimCpp) {
  function("iterSim", &iterSim);
}
