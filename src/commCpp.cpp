#pragma once

#include <RcppArmadillo.h>
#include "roleParamsCpp.cpp"

using namespace Rcpp;
using namespace arma; 

// [[Rcpp::depends(RcppArmadillo)]]

//' @name metaCommCpp
//' @title a class to specify the meta community
//' @param abundance a \code{numeric vector} of relative abundances for each species
//' @param traits a \code{numeric vector} of traits for each species
//' @param Smax a single \code{integer} specifying the total number of species ever
//' recorded in the local community (both locally extinct and extant)

class metaCommCpp {
    private:
    public:
        NumericVector abundance;
        NumericVector traits;
        int Smax;
  
        //constructor
        metaCommCpp(NumericVector abundance_, NumericVector traits_,int Smax_)
        {
          abundance = abundance_; 
          traits = traits_;
          Smax = Smax_; 
        }
};

//' @name localCommCpp
//' @title a class to specify the local community
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
        int J; // constant number of individuals in the community
        int Imax; // max number of individuals in the community
        int Smax; // max number of species in community, index used to create new species ids 
        NumericVector pi; // genetic diversities of species - unsure how this plays into new structure
  
        arma::mat traitdiffs; // a matrix that is the outer product of traits*traits 
                                  // calculated at object creation then altered rather than recalculated
        //int traitMax; don't think we need this anymore? 
        NumericVector abundance_sp; // abundances of every species i at index i, kept to save on computation
        NumericVector traits_sp; // traits of every species i at index i, kept to save on computation
        
        bool print;
        
        // constructor takes species abundances and species traits and decollapses
        // them into a vector of individuals
        localCommCpp(NumericVector abundance_, NumericVector traits_, int Smax_,
                     NumericVector pi_): abundance_sp(abundance_), traits_sp(traits_), Smax(Smax_), pi(pi_)
        {
            // init abundance, species, traits vectors
            abundance_indv = NumericVector(10000);  
            species_ids = NumericVector(10000);
            traits = NumericVector(10000); 
            
            // init Imax
            Imax = 0; 
            
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
            
            // init J, which is the same as Imax at init 
            J = Imax;
            
            // trim abundance, species, traits vectors to J 
            abundance_indv = abundance_indv[Rcpp::Range(0,J-1)];
            species_ids = species_ids[Rcpp::Range(0,J-1)];
            traits = traits[Rcpp::Range(0,J-1)];
            
            Rcout << "init traitdiffs" << "\n";
            // init traitdiffs as outer product of traits 
            //arma::mat tr = as<arma::mat>(traits_);
            //traitdiffs = arma::sp_mat(tr*tr);
            traitdiffs = arma::mat(J,J);
            
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
        
        // birth is called on indv i, the new indv REPLACING the index at r
        void birth(int i, int r)
        {
            // add individual 
            abundance_indv[r] = 1;
          
            // indv gets the species and trait of the indv that gave birth to it 
            species_ids[r] = species_ids[i];
            traits[r] = traits[i]; // soon traits should vary randomly when birth is called
            
            // increment the abundance of the species holding the indv that gave birth
            abundance_sp[species_ids[i]] = abundance_sp[species_ids[i]] + 1;
            
            if(print){Rcout << "incrementing species " << species_ids[i] << "\n";}
            
            // update traitdiffs
            arma::vec v = (traits * -1) + traits[r];
            
            traitdiffs.row(r) = v.t(); 
            
            // increment Imax
            Imax += 1; 
        }
        
        // death is called on individual i 
        void death(int i)
        {
            // change binary abundance from 1 to 0
            abundance_indv[i] = 0;
          
            // decrement the abundance of the species holding the indv that gave birth
            abundance_sp[species_ids[i]] = abundance_sp[species_ids[i]] - 1;
        }
        
        // speciation is called on species s, the new individual REPLACING the index at r
        void speciation(int s, int r, roleParamsCpp p)
        {
            // calculate random trait deviation by sigma
            float tdev = R::rnorm(0, p.values.trait_sigma);
            
            // add abundance, species id and trait for new individual of species
            abundance_indv[r] = 1;
            species_ids[r] = s; 
            traits[r] = traits_sp[Smax] + tdev; 
            
            // update traitdiffs
            //traitdiffs(Imax,_) = traits[Imax] - traits;
            arma::rowvec v = traits * -1 + traits[r];
            traitdiffs.row(r) = v; 
            
            // add abundance and trait for new species
            abundance_sp[Smax] = 1;
            
            // TODO calculate new species trait here as a function of individual traits
            traits_sp[Smax] = traits_sp[s] + tdev; 
            
            // increment num indv and sp
            Imax += 1; 
            Smax += 1;
        }
        
        // immigration is called on species s 
        void immigration(int s, int r, metaCommCpp m)
        {
          // add indv of species s
          abundance_indv[r] = 1;
          species_ids[r] = s;
          traits[r] = m.traits[s];
          
          // increment the abundance of the species holding the indv that gave birth
          abundance_sp[s] = abundance_sp[s] + 1;
          
          if(print){Rcout << "incrementing species " << s << "\n";}
          
          // update traitdiffs
          arma::rowvec v = traits * -1 + traits[s];
          traitdiffs.row(r) = v; 
          
          // increment Imax
          Imax += 1; 
        }

};

RCPP_EXPOSED_CLASS(localCommCpp)
RCPP_EXPOSED_CLASS(metaCommCpp)

