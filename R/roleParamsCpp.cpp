#include <Rcpp.h>
#include <map>
#include <paramValuesCpp.cpp>

using namespace Rcpp;

class roleParamsCpp {
public:
    paramValuesCpp values;
    int n;
    std::string runType;
    // todo - add priors

    //constructor
    roleParamsCpp(paramValuesCpp values_, int n_, std::string runType_) :
        values(values_), n(n_), runType(runType_)
    {
    }
};

RCPP_EXPOSED_CLASS(roleParamsCpp)

RCPP_MODULE(paramsCpp) {
    class_<roleParamsCpp>("roleParamsCpp")
        .constructor<paramValuesCpp,int,std::string>()
        .field("values", &roleParamsCpp::values)
        ;
}
