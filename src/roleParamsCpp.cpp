#pragma once

#include <Rcpp.h>
#include "paramValuesCpp.cpp"

using namespace Rcpp;

class roleParamsCpp {
public:
    paramValuesCpp values;
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
