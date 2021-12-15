#pragma once

#include <RcppArmadillo.h>
#include "roleParamsCpp.cpp"

using namespace Rcpp;
using namespace arma; 

// [[Rcpp::depends(RcppArmadillo)]]

//' @name metaCommCpp
//' @title a C++ class to specify the meta community
//' @field new Constructor
//' @param abundance a numeric vector of relative abundances for each species
//' @param traits matrix of traits; the first column specifies the species index
//' (i.e. the index of that species in the \code{abundance} vector) and the
//' subsequent columns specify the trait values
//' @param Smax a single integer specifying the total number of species ever
//' recorded in the local community (both locally extinct and extant)

class metaCommCpp {
    private:
    public:
        NumericVector abundance;
        NumericMatrix traits;
        int Smax;
  
        //constructor
        metaCommCpp(NumericVector abundance_, NumericMatrix traits_,int Smax_)
        {
          abundance = abundance_; 
          traits = traits_;
          Smax = Smax_; 
        }
};

//' @name localCommCpp
//' @title a C++ class to specify the local community
//' @field new Constructor
//' @param abundance_binary a numeric vector specifying the binary abundance (alive or dead) 
//' of each individual
//' @param traits a numeric vector of trait values for each individual
//' @param Imax a float specifying the total number of individuals ever
//' recorded in the local community (both locally extinct and extant)
//' @param pi a numeric vector of genetic diversities for each species

class localCommCpp {
    private:
    public:
        NumericVector abundance_indv; // a vector of 0s and 1s specifying alive or dead individuals
        NumericVector species_ids; // the species num that each individual belongs to
        NumericVector traits; // trait values 
        int Imax; // max number of individuals in the community, index used often
        int Smax; // max number of species in community, index used to create new species ids 
        NumericVector pi; // genetic diversities of species - unsure how this plays into new structure
  
        arma::sp_mat traitdiffs; // a sparse matrix that is the outer product of traits*traits 
                                  // calculated at object creation then altered rather than recalculated
        //int traitMax; don't think we need this anymore? 
        NumericVector abundance_sp; // abundances of every species i at index i, kept to save on computation
        NumericVector traits_sp; // traits of every species i at index i, kept to save on computation
        
        // constructor takes species abundances and species traits and decollapses
        // them into a vector of individuals
        localCommCpp(NumericVector abundance_, NumericVector traits_, int Smax_,
                     NumericVector pi_): abundance_sp(abundance_), traits_sp(traits_)
        {
            // init abundance, species, traits vectors
            abundance_indv = NumericVector(10000);  
            species_ids = NumericVector(10000);
            traits = NumericVector(10000); 
            
            // init Imax
            Imax = 0; 
            // init Smax
            Smax = Smax_; 
            
            Rcout << "started sp decollapse" << "\n";

            // for every species
            for(int s = 0; s < abundance_.length(); s++)
            {
              // for every individual in species
              for(int i = 0; i < abundance_[s]; i++)
              {
                // add that individual to vector of individuals
                abundance_indv[Imax] = 1;

                // add species id
                species_ids[Imax] = s;

                // add trait value
                traits[Imax] = traits_[i]; //maybe deviate indv traits from sp randomly later?

                // increment Imax
                Imax += 1;
              }
            }
            pi = pi_;
            
            Rcout << "init traitdiffs" << "\n";
            
            // init traitdiffs as outer product of traits 
            //arma::mat tr = as<arma::mat>(traits_);
            //traitdiffs = arma::sp_mat(tr*tr);
            traitdiffs = arma::sp_mat(10000,10000);
            
            Rcout << "diffs dims " << arma::size(traitdiffs) << "\n";
            Rcout << "started traitdiffs calc" << "\n";
            Rcout << "traits length = " << traits_.length() << "\n";

            // add all current traitdiffs
            for(int i = 0; i < traits_.length(); i++)
            {
              for(int j = 0; j < traits_.length(); j++)
              {
                traitdiffs(i,j) = traits_[i] - traits_[j];
              }
            }

            Rcout << "diffs dims " << arma::size(traitdiffs) << "\n";

            //traitdiffs = traitdiffs_;
            //abundance_sp = abundance_;
        }
        
        // birth is called on individual i 
        void birth(int i)
        {
            // add individual 
            abundance_indv[Imax] = 1;
            // indv gets the species and trait of the indv that gave birth to it 
            species_ids[Imax] = species_ids[i];
            traits[Imax] = traits[i]; // should traits vary randomly when birth is called? 
            
            // update traitdiffs
            arma::vec v = traits * -1 + traits[Imax];
            
            traitdiffs.row(Imax) = v.t(); 
            
            //traitdiffs(Imax,) = traits[Imax] - traits;
            
            Imax += 1; 
        }
        
        // death is called on individual i 
        void death(int i)
        {
            abundance_indv[i] = 0;
            // decrement species abundances 
            abundance_sp[species_ids[i]] -= 1; 
        }
        
        // speciation is called on species s
        void speciation(int s, roleParamsCpp p)
        {
            // calculate random trait deviation by sigma
            float tdev = R::rnorm(0, p.values.trait_sigma);
            
            // add abundance and trait for new individual of species
            abundance_indv[Imax] = 1;
            traits[Imax] = traits_sp[Smax] + tdev; 
            
            // update traitdiffs
            //traitdiffs(Imax,_) = traits[Imax] - traits;
            arma::rowvec v = traits * -1 + traits[Imax];
            traitdiffs.row(Imax) = v; 
            
            // add abundance and trait for new species
            abundance_sp[Smax] = 1;
            traits_sp[Smax] = traits_sp[s] + tdev; 
            // NOTE - may have to update species trait means whenever birth occurs if keeping this 
            
            // increment num indv and sp
            Imax += 1; 
            Smax += 1;
            
            // NOTE - should speciation dev be of mean trait? or of a random individual in the species instead? 
        }
        
        // immigration is called on species s 
        void immigration(int s, metaCommCpp m)
        {
          // add indv of species s
          abundance_indv[Imax] = 1;
          species_ids[Imax] = s;
          traits[Imax] = m.traits[s];
          
          // update traitdiffs
          //traitdiffs(Imax,_) = traits[Imax] - traits;
          arma::rowvec v = traits * -1 + traits[Imax];
          traitdiffs.row(Imax) = v; 
          
          Imax += 1; 
        }

};

RCPP_EXPOSED_CLASS(localCommCpp)
RCPP_EXPOSED_CLASS(metaCommCpp)

