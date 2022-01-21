#pragma once

#include <RcppArmadillo.h>
#include "paramValuesCpp.cpp"
#include <map> 

using namespace Rcpp;

//' @name roleParamsCpp
//' @title a C++ class holding the model parameter values as a paramValuesCpp object plus some additional metadata
//' @param values a paramValuesCpp object containing the parameter values proper 
//' @param runType the type of model run, either "sim" or "fit"

class roleParamsCpp {
public:
    paramValuesCpp values; // rename as current values
    std::map <std::string,double> vals; 
    std::string runType;
    
    //placeholder for priors, which contain distributions from which to draw parameters
    int pr;

    //constructor
    roleParamsCpp(paramValuesCpp values_, std::string runType_, int pr_) :
        values(values_), runType(runType_), pr(pr_)
    {
    }
};

RCPP_EXPOSED_CLASS(roleParamsCpp)
