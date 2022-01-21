#pragma once

#include <RcppArmadillo.h>
#include "paramValuesCpp.cpp"
#include <map> 

using namespace Rcpp;

//' @name paramDistributionCpp 
//' @title a class holding the values to create distribution from which to draw parameters
//' @param distrType a string naming the type of distribution

class paramDistributionCpp {
public:
  std::string distrType;
  double val1;
  double val2; 
  
  //constructor
  paramDistributionCpp(std::string distrType_, double val1_, double val2_) :
    distrType(distrType_), val1(val1_), val2(val2_)
  {
  }
  
  double sample()
  {
    if(distrType == "uniform")
    {
      return Rcpp::runif(1,val1,val2)[1];
    }
    
    if(distrType == "normal")
    {
      return Rcpp::rnorm(1,val1,val2)[1];
    }
    
    return -1;
  }
};

RCPP_EXPOSED_CLASS(roleParamsCpp)
  