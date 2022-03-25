#pragma once

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma; 

// [[Rcpp::depends(RcppArmadillo)]]

//' @name metaCommCpp
//' @title a class to specify the meta community
//' @field abundance a \code{numeric vector} of relative abundances for each species - stays static 
//' @field traits a \code{numeric vector} of traits for each species
//' @field Smax a single \code{integer} specifying the total number of species ever
//' recorded in the local community (both locally extinct and extant)

class metaCommCpp {
    private:
    public:
        NumericVector abundance;
        NumericVector traits;
  
        //constructor
        metaCommCpp(NumericVector abundance_, NumericVector traits_)
        {
          abundance = abundance_; 
          traits = traits_;
        }
};

//' @name localCommCpp
//' @title a class to specify the local community
//' @field abundance_indv a numeric vector of 0s and 1s specifying the binary abundance (alive or dead) 
//' of each individual
//' @field species_ids a numeric vector specifying species of each individual as
//' an index found in the species vectors 
//' @field traits a numeric vector of trait values for each individual
//' @field J constant number of individuals in the community
//' @field Smax max number of species in community, index used to create new species ids 
//' @field traitdiffs a matrix that is the outer product of traits*traits - 
//' calculated at object creation then altered rather than recalculated
//' @field abundance_sp abundances of every species i at index i, kept to save on computation
//' @field traits_sp traits of every species i at index i, kept to save on computation
//' @field pi_sp genetic diversities of species - unsure how this plays into new structure
//' @field print a bool specifying whether to print outputs related to class functions 

class localCommCpp {
    private:
    public:
        NumericVector abundance_indv; 
        NumericVector species_ids;
        NumericVector traits; 
        int J;
  
        arma::mat traitdiffs;
                                  
        NumericVector abundance_sp;
        NumericVector traits_sp;
        NumericVector pi_sp;
        int S_index;
        
        bool print;
        
        // constructor takes species abundances and species traits and decollapses
        // them into a vector of individuals
        localCommCpp(NumericVector abundance_, NumericVector traits_, NumericVector pi_, int aug_length): 
          abundance_sp(abundance_), traits_sp(traits_), pi_sp(pi_)
        {
            // init S_index
            S_index = abundance_.length() - 1;
  
            // init J as the sum of species abundances
            J = sum(abundance_sp);
            abundance_indv = NumericVector(J);
            species_ids = NumericVector(J);
            traits = NumericVector(J);
            
            if(print){Rcout << "started sp decollapse" << "\n";}
            
            // for every species
            int k = 0;
            for(int s = 0; s < abundance_sp.length(); s++)
            {
              // for every individual in species
              for(int i = 0; i < abundance_sp(s); i++)
              {
                // add that individual to vector of individuals
                abundance_indv(k) = 1;

                // add species id
                species_ids(k) = s;

                // add trait value
                traits(k) = traits_(s); //maybe deviate indv traits from sp randomly later?

                // increment counter
                k += 1; 
              }
            }
            
            if(print){Rcout << "init traitdiffs" << "\n";}
            
            // init traitdiffs, which will be equal to the outer product of traits 
            traitdiffs = arma::mat(J,J);
            
            if(print){Rcout << "diffs dims " << arma::size(traitdiffs) << "\n";}
            if(print){Rcout << "started traitdiffs calc" << "\n";}
            if(print){Rcout << "traits length = " << traits_.length() << "\n";}

            // add all current traitdiffs
            for(int i = 0; i < traits_.length(); i++)
            {
              for(int j = 0; j < traits_.length(); j++)
              {
                traitdiffs(i,j) = traits_(i) - traits_(j);
              }
            }
            
            // copy abundance_sp, traits_sp, and pi_species to aug_length vectors
            abundance_sp = augmentVector(abundance_sp,aug_length);
            traits_sp = augmentVector(traits_sp, aug_length);
            pi_sp = augmentVector(pi_sp, aug_length);
        }
        
        // constructor used to copy without augmenting or flipping species to indv 
        localCommCpp(NumericVector abundance_sp_, NumericVector traits_sp_, NumericVector pi_sp_, 
                     NumericVector abundance_indv_, NumericVector traits_): 
                     abundance_sp(abundance_sp_), traits_sp(traits_sp_), pi_sp(pi_sp_),
                     abundance_indv(abundance_indv_), traits(traits_)
        {
        }
        
        // birth is called on indv i, the new indv REPLACING the index at r
        void birth(int i, int r)
        {
            // add individual 
            abundance_indv(r) = 1;
          
            // indv gets the species and trait of the indv that gave birth to it 
            species_ids(r) = species_ids(i);
            //traits[r] = traits[i] + rnorm(0,1)
            //rnorm(0, sqrt(sigma_bm / (1 / speciation + 1/extinction) / J)
                    
            traits(r) = traits(i); // soon traits should vary randomly when birth is called
            
            // increment the abundance of the species holding the indv that gave birth
            abundance_sp(species_ids(i)) = abundance_sp(species_ids(i)) + 1;
            
            if(print){Rcout << "incrementing species " << species_ids(i) << "\n";}
            //CRASHES HERE
            
            // update traitdiffs
            arma::vec v = (traits * -1) + traits(r);
            
            traitdiffs.row(r) = v.t(); 
            //traitdiffs.col(r) = v.t(); 
            
            // update traits_sp 
            updateSpeciesTraitMean(species_ids[i]);
        }
        
        // death is called on individual i 
        void death(int i)
        {
            // change binary abundance from 1 to 0
            abundance_indv(i) = 0;
          
            // decrement the abundance of the species holding the indv that died
            abundance_sp(species_ids(i)) = abundance_sp(species_ids(i)) - 1;
            
            if(print){Rcout << "updating traits_sp" << "\n";}

            // update traits_sp 
            updateSpeciesTraitMean(species_ids[i]);
        }
        
        // speciation is called on species s (the species UNDERGOING speciation), 
        //   the new individual REPLACING the index at r
        // new species is placed at Smax 
        void speciation(int s, int r, double trait_sigma)
        {
          
            // increment species index 
            S_index += 1;
          
            // calculate random trait deviation by sigma
            float tdev = R::rnorm(0, trait_sigma);
            //tdev is also scaled by branch length, same scaling but just not divided by J 
            
            if(print){Rcout << "adding info for new indv" << "\n";}
            
            // add abundance, species id and trait for new individual of species
            abundance_indv(r) = 1;
            species_ids(r) = S_index; 
            traits(r) = traits_sp(s) + tdev; 
            
            
            if(print){Rcout << "updating traitdiffs" << "\n";}
            
            // update traitdiffs
            arma::rowvec v = traits * -1 + traits(r);
            traitdiffs.row(r) = v; 
            
            if(print){Rcout << "adding info for new species" << "\n";}
            
            // add abundance and trait for new species
            abundance_sp[S_index] = 1; 
            
            traits_sp[S_index] = traits(r);
        }
        
        // immigration is called on species s 
        void immigration(int s, int r, metaCommCpp m)
        {
          // add indv of species s
          abundance_indv(r) = 1;
          species_ids(r) = s;
          traits(r) = m.traits(s);
          
          // increment the abundance of the species holding the indv that gave birth
          abundance_sp(s) = abundance_sp(s) + 1;
          
          if(print){Rcout << "incrementing species " << s << "\n";}
          
          // update traitdiffs
          arma::rowvec v = traits * -1 + traits(s);
          traitdiffs.row(r) = v; 
          
          // update traits_sp 
          updateSpeciesTraitMean(s);
        }
        
        // update the traits of a given species id 
        void updateSpeciesTraitMean(int id)
        {
          
          // subtract old trait / N, add new trait / N 
          if(print){Rcout << "updating species trait mean" << "\n";}
          if(print){Rcout << "id " << id << "\n";}
          LogicalVector matches = LogicalVector(traits.size());
          
          
          for(int i = 0; i < matches.size(); i++)
          {
            if(species_ids[i] == id){
              matches[i] = true;
            }
            else{
              matches[i] = false;
            }
          }
          if(print){Rcout << "matches " << matches << "\n";}
          NumericVector trait_matches = traits[matches];
          if(print){Rcout << "trait_matches " << trait_matches << "\n";}
          double trait_mean = Rcpp::mean(trait_matches);
          if(print){Rcout << "trait_mean " << trait_mean << "\n";}
          traits_sp[id] = trait_mean;
        }
        
        // return vector v augmented in size by aug_length (used on abundance_sp, traits_sp, pi_sp)
        NumericVector augmentVector(NumericVector v, int aug_length)
        {
          NumericVector out = NumericVector(aug_length);
          for(int i = 0; i < v.length(); i++)
          {
            out[i] = v[i];
          }
          return(out);
        }

};

RCPP_EXPOSED_CLASS(localCommCpp)
RCPP_EXPOSED_CLASS(metaCommCpp)

