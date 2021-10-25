#pragma once
#include <Rcpp.h>
#include "roleParamsCpp.cpp"

using namespace Rcpp;

class commCpp {
private:

public:
    NumericVector abundance;
    NumericMatrix traits;
    int Smax;

    //constructor
    commCpp(NumericVector abundance_, NumericMatrix traits_, int Smax_)
            : abundance(abundance_), traits(traits_), Smax(Smax_)
    {
    }
};

//' @name localCommCpp
//' @title a C++ class to specify the local community
//' @field new Constructor
//' @param abundance a numeric vector of abundances for each species
//' @param traits matrix of traits; the first column specifies the species index
//' (i.e. the index of that species in the \code{abundance} vector) and the
//' subsequent columns specify the trait values
//' @param Smax a single integer specifying the total number of species ever
//' recorded in the local community (both locally extinct and extant)
//' @details Smax is used for internal bookkeeping.  The dimensions of
//' \code{abundance} and \code{traits} can be greater than
//' \code{Smax}.  In such cases, \code{Smax} is used to index where to begin
//' adding new information (e.g. a new species can be added at the index
//' \code{Smax + 1}).
//' @param pi a numeric vector of genetic diversities for each species

class localCommCpp : public commCpp {
    private:
    public:
        NumericVector pi;
        int traitMax;

        //constructor
        localCommCpp(NumericVector abundance_, NumericMatrix traits_, int Smax_,
                     NumericVector pi_)
                    : commCpp(abundance_, traits_, Smax_)
        {
            pi = pi_;
            //I believe traitMax can start as Smax_ - check this
            traitMax = Smax_;
        }

        void birth(int i)
        {
            abundance[i] += 1;
        }

        void death(int i)
        {
            abundance[i] -= 1;
        }

        void speciation(int i, roleParamsCpp p)
        {
            // update Smax and initialize abundance
            // update number of species
            Smax += 1;

            // initialize abundance for new species
            abundance[Smax] = 1;

            // index of where unrealized traits begin
            int j = traitMax + 1;

            // add trait
            traits(j, 1) = Smax;
            traits(j, 2) = traits(i, 2) + R::rnorm(0, p.values.trait_sigma);

            //increase traitMax
            traitMax += 1;
        }

        void immigration(int i)
        {
            birth(i);
        }

};

//' @name metaCommCpp
//' @title a C++ class to specify the meta community
//' @field new Constructor
//' @param abundance a numeric vector of abundances for each species
//' @param traits matrix of traits; the first column specifies the species index
//' (i.e. the index of that species in the \code{abundance} vector) and the
//' subsequent columns specify the trait values
//' @param Smax a single integer specifying the total number of species ever
//' recorded in the local community (both locally extinct and extant)

class metaCommCpp : public commCpp {
    private:

    public:
        NumericVector abundance;
        NumericMatrix traits;
        int Smax;

        //constructor
        metaCommCpp(NumericVector abundance_, NumericMatrix traits_,int Smax_)
                    : commCpp(abundance_, traits_, Smax_)
        {
        }
};

RCPP_EXPOSED_CLASS(localCommCpp)
RCPP_EXPOSED_CLASS(metaCommCpp)

