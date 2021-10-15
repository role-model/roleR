#include <Rcpp.h>
#include <commCpp2.R>
#include <roleModelCpp.>

using namespace Rcpp;

localCommCpp death(localCommCpp x)
{
    x.abundance[i] += 1;
    return x;
}

roleModelCpp birth(roleModelCpp x)
{
    i <- sample(x@localComm@Smax, 1,
                prob = x@localComm@abundance[1:x@localComm@Smax])

    x@localComm <- birth(x@localComm, i)

    return(x)
}
