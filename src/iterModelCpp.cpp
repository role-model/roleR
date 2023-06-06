// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "roleDataCpp.cpp"
#include "roleParamsCpp.cpp"
#include <iostream>

using namespace Rcpp;

// iterModelCpp.cpp contains:
//      + Cpp functions used within the core loop of iterModelCpp
//          1. sample_zero_to_x
//          2. sample_index_using_probs
//          3. call_birth
//          4. call_dispersal
//          5. call_speciation
//
//      + iterModelCpp, a function exported to R that iterates a model
//      + Wrappers around the loop C++ functions exported to R for testing - this structure avoids have to redundantly wrap every function
//          1. intFunCpp
//          2. dataFunCpp
//          3. vectFunCpp

// sample a random integer from 0 to x
// replicates sample(1:x, 1)
int sample_zero_to_x(int x)
{
    return((int) (R::runif(0,1) * (double) x));
}

// sample an index using a vector of relative probabilities
int sample_index_using_probs(NumericVector probs){
    IntegerVector v = Rcpp::sample(probs.length(), 1, false, probs, false);
    return(v(0));
}

// get prob of death for each individual due to env filtering
// raising traits - optimum to power of 2
// multiply by gaussian kernel of env_sigma
// ANDY  - could save this to avoid recomputing
//  this would be a new function similar to update_trait_diffs_sq
// JACOB NOTE - this actually happens currently - every iter the prob is recomputed ONLY at the dead_indv location
NumericVector get_filtering_death_probs(int i, roleDataCpp &d, roleParamsCpp &p){
    return(1 - exp(-1/p.env_sigma[0] * pow(d.indTraitL - 0, 2)));
}

// get prob of death for each individual due to competition
NumericVector get_comp_death_probs(int i, roleDataCpp &d, roleParamsCpp &p){
    arma::colvec comp_probs_a = arma::sum(arma::exp((-1/p.comp_sigma(i)) * d.traitDiffsSq),1) / p.individuals_local[0];
    return(Rcpp::wrap(comp_probs_a.begin(),comp_probs_a.end()));
}

// get combined probs of death from filtering and competition, weighted by degree of neutrality
NumericVector get_combined_death_probs(int i, roleParamsCpp &p, NumericVector f_probs, NumericVector c_probs){
    return(p.neut_delta[0] + (1-p.neut_delta[0]) * ((p.env_comp_delta[0] * f_probs) + ((1-p.env_comp_delta[0]) * c_probs)));
}

// check if extinction happens and if it does update the saved species extinction steps & the alive vector in the phylo
void update_extinct(int i, roleDataCpp &d, roleParamsCpp &p, int dead_index){ // - uses dead_index as dead species index
    
    int dead_species = dead_index;
    // check for extinction
    // extinct if species abundance of dead individual is now <0 
    bool extinct_in_local = d.spAbundL(dead_species) <= 0;
    // not in meta if species of the dead individual is not in meta
    bool not_in_meta = dead_species > d.spAbundM.length();
    
    // check for local extinction
    if(extinct_in_local){
        d.spExtinctionStepL(dead_species) = i;
    }
    
    // check for local and meta extinction
    if(extinct_in_local & not_in_meta){
        d.aliveP(dead_species) = false;
    }
}

// call death on a roleDataCpp object using the params at iteration i
// returns the individual chosen for death which is needed in the rest of the loop
int call_death(int i, roleDataCpp &d, roleParamsCpp &p){
    
    int n_indv = p.individuals_local[0];
    // start with neutral death probs where all probs are equal: rep(1,n_indv)
    NumericVector death_probs(n_indv,1.0);
    
    // if not neutral, calculate initial probs of death from environmental filtering
    if(p.neut_delta[0] != 1)
    {
        // ANDY  - figure out how ifs should be structured
        // JACOB NOTE - done!
        
        // envFilterProbs is updated at the END of every loop iter to be current for the next iter
        // comp probs are NOT updated, but the most expensive part of it which ALWAYS changes every iter (traitDiffsSquared),
        //  is updated at the end of every loop iter
        
        // if env_comp_delta is 1, only env filtering matters, so don't worry about combining
        //  and death probs is just the current envFilterProbs
        if(p.env_comp_delta[0] == 1){
            death_probs = d.envFilterProbs;
        }
        else if(p.env_comp_delta[0] == 0){ // if is 0, only competition matters
            // METHOD - calculate probs of death due to competition
            death_probs = get_comp_death_probs(i,d,p);
        }
        else{ // otherwise, must combine
            // METHOD - calculate combined death probs
            death_probs = get_combined_death_probs(i,p,d.envFilterProbs,get_comp_death_probs(i,d,p));
        }
        
        // save the used env sigma as the prevEnvSigma for the next loop
        // JACOB NOTE - should neutrality or env comp delta ever be time-varying? currently they cannot be but 
        //   but minor refactoring would make it possible
        // ANDY NOTE - based on conceptual principles decide which things can time vary and which can't
        //   maybe hackathon topic - but keep track of which params are invariant
        d.prevEnvSigma = p.env_sigma[0];
    }
    
    // get dead index and dead species
    int dead_index = sample_index_using_probs(death_probs);
    int dead_species = d.indSpeciesL(dead_index);

    // subtract from abundance of dead species
    d.spAbundL(dead_species) = d.spAbundL(dead_species) - 1;
    
    // check if extinct and update phylo
    update_extinct(i,d,p,dead_species);
    
    return(dead_index);
}

// passing explicitly by reference using '&' makes things as fast or almost as fast as if code is pasted into loop
// given a presampled parent, sets the sp of the new individual, adds to the abundance of the sp matrix, 
// and calculate and add new trait value
void call_birth(int i, int dead_index, int parent_indv, roleDataCpp &d, roleParamsCpp &p, bool print){
    
    // get the birthed species 
    int birthed_species = d.indSpeciesL(parent_indv); 
    
    // set the species of the new indv to that of the parent 
    d.indSpeciesL(dead_index) = birthed_species;

    // add to the abundance of the species matrix
    d.spAbundL(birthed_species) = d.spAbundL(birthed_species) + 1; 
    
    // calculate trait change from parent
    // ANDY NOTE - change rnorm and unif to use std C++
    // JACOB NOTE - done! with different distr tnorm & sptnorm
    // but note that if params affecting the norm for sampling traits vary over iters,
    //  the distr object will have to be remade and reassigned as the stdev cannot be edited
    //  however, now that I finished this I wonder if we should just stick with Rcpp sugar i.e. Rcpp:rnorm
    //  note that Rcpp sugar DOES NOT call R, rather it is a C++ implementation
    //  we may be reinventing the wheel in a more ugly way rather than making things faster - can test if unsure!
    
    float trait_change = d.norm(d.rng) * p.trait_sigma(i) / (p.speciation_meta(i) + p.extinction_meta(i));

    // add new trait value
    d.indTraitL(dead_index) = d.indTraitL(parent_indv) + trait_change;
}

// set the species of the new indv, update last time of origin, add to sp abundance matrix, and calculate and add new trait
void call_dispersal(int i, int dead_index, int parent_indv, roleDataCpp &d, roleParamsCpp &p, bool print){
    
    // sample a parent index using species abundances as probs
    int parent_index = parent_indv;

    // set the species to the parent species from meta
    d.indSpeciesL(dead_index) = parent_index;
    
    // update time of last origin if the species' last local abundance was 0 
    if(d.spAbundL(d.indSpeciesL(dead_index)) <= 0){
        d.spLastOriginStepL(d.indSpeciesL(dead_index)) = i;
    }
    
    // add to the abundance of the species matrix
    d.spAbundL(d.indSpeciesL(dead_index)) = d.spAbundL(d.indSpeciesL(dead_index)) + 1;
    
    // calculate trait change from parent
    float trait_change = d.norm(d.rng) * p.trait_sigma(i) / (p.speciation_meta(i) + p.extinction_meta(i));
    
    // add new trait value
    d.indTraitL(dead_index) = d.spTraitM(parent_index) + trait_change;
    
    // set founderFlag of dead indv to 1 
    d.founderFlagL(dead_index) = 1; 
}

// get the probs of speciation, determines by weighted local and meta abundance
NumericVector get_speciation_probs(int i, roleDataCpp &d, roleParamsCpp &p){
    
    // get dispersal prob
    double dp = p.dispersal_prob[i];
    
    // calculate speciation probs as sum of local and meta abundances weighted by dispersal prob
    NumericVector probs = dp * (d.spAbundM / sum(d.spAbundM)) + 
        (1 - dp) * (d.spAbundL / sum(d.spAbundL));
    
    // normalize 
    probs = (abs(probs)+probs)/2;
    
    return(probs);
}

// // helper functions to augment vectors if an index out of range is about to happen despite the buffering
// Rcpp::NumericVector bufferNumericVector(Rcpp::NumericVector v, int value, int n) {
//     int vsize = v.size(); // get the size of vec1
//     Rcpp::NumericVector result(vsize + n, value); // create a new vector of the appropriate size
//     std::copy(v.begin(), v.end(), result.begin()); // copy the elements of vec1 to the beginning of result
//     std::copy(v.begin(), v.end(), result.begin() + vsize); // copy the elements of vec2 to the end of result
//     return result;
// }
Rcpp::NumericVector buffer_numeric_vector(Rcpp::NumericVector vec, int value, int n) {
    int len = vec.length() + n; // calculate the new length of the vector
    Rcpp::NumericVector result(len); // create a new vector with the new length
    std::fill(result.begin(), result.end(), value); // fill the new vector with zeroes
    std::copy(vec.begin(), vec.end(), result.begin() + n); // copy the original vector to the new vector, starting at the nth position
    return result; // return the new vector
}
Rcpp::LogicalVector buffer_logical_vector(Rcpp::LogicalVector vec, bool value, int n) {
    int len = vec.length() + n; // calculate the new length of the vector
    Rcpp::LogicalVector result(len); // create a new vector with the new length
    std::fill(result.begin(), result.end(), value); // fill the new vector with zeroes
    std::copy(vec.begin(), vec.end(), result.begin() + n); // copy the original vector to the new vector, starting at the nth position
    return result; // return the new vector
}
Rcpp::CharacterVector buffer_character_vector(Rcpp::CharacterVector vec, std::string value, int n) {
    int len = vec.length() + n; // calculate the new length of the vector
    Rcpp::CharacterVector result(len); // create a new vector with the new length
    std::fill(result.begin(), result.end(), value); // fill the new vector with zeroes
    std::copy(vec.begin(), vec.end(), result.begin() + n); // copy the original vector to the new vector, starting at the nth position
    return result; // return the new vector
}
arma::imat buffer_arma_mat(arma::imat imat, int value, int n) {
    //arma::mat mat = arma::conv_to<arma::mat>::from(imat);
    //int n_rows = mat.n_rows + n; // calculate the new number of rows of the matrix
    //int n_cols = mat.n_cols; // keep the same number of columns
    //arma::mat result(n_rows, n_cols, arma::fill::zeros); // create a new matrix with the new number of rows, filled with zeroes
    //result.submat(n, 0, n_rows - 1, n_cols - 1) = mat; // copy the original matrix to the new matrix, starting at the (n+1)th row
    //result.submat(0, 0, n - 1, n_cols - 1).fill(value); // fill the top n rows of the new matrix with -1
    arma::mat result = arma::join_cols(arma::conv_to<arma::mat>::from(imat), arma::ones<arma::mat>(n, 2) * value);
    
    return(arma::conv_to<arma::imat>::from(result)); // return the new matrix
}
arma::vec buffer_arma_vec(arma::vec v, int value, int n) {
    v.resize(v.n_elem + n);
    v.tail(n).fill(value);
    return(v); 
}
// Rcpp::LogicalVector bufferLogicalVector(Rcpp::LogicalVector v, bool value, int n) {
//     int vsize = v.size(); // get the size of vec1
//     Rcpp::LogicalVector result(vsize + n, value); // create a new vector of the appropriate size
//     std::copy(v.begin(), v.end(), result.begin()); // copy the elements of vec1 to the beginning of result
//     std::copy(v.begin(), v.end(), result.begin() + vsize); // copy the elements of vec2 to the end of result
//     return result;
// }
// Rcpp::CharacterVector bufferCharacterVector(Rcpp::CharacterVector v,  std::string value, int n) {
//     int vsize = v.size(); // get the size of vec1
//     Rcpp::CharacterVector result(vsize + n,value); // create a new vector of the appropriate size
//     std::copy(v.begin(), v.end(), result.begin()); // copy the elements of vec1 to the beginning of result
//     std::copy(v.begin(), v.end(), result.begin() + vsize); // copy the elements of vec2 to the end of result
//     return result;
// }

// call the local and meta aspects of speciation
void update_speciation_local_meta(int i, int dead_index, roleDataCpp &d, roleParamsCpp &p, bool dispersed_this_iter, int speciation_sp, bool print){
    
    //if(print){ Rcout << "updating speciation local meta"<< "\n";}
    
    // augment vectors by a set amount if the id of the new species (the n tips in the phylo)
    //  exceeds the size of the local sp vectors
    if(d.nTipsP(0) > d.spAbundL.length() - 1){
        if(print){ Rcout << "buffering in local meta" << "\n";}
        d.spAbundL = buffer_numeric_vector(d.spAbundL,0,100);
        d.spTraitL = buffer_numeric_vector(d.spTraitL,0,100);
        d.spLastOriginStepL = buffer_numeric_vector(d.spLastOriginStepL,0,100);
        d.founderFlagL = buffer_numeric_vector(d.founderFlagL,0,100);
    }
    
    //if(print){ Rcout << "indTraitL" << "\n";}
    
    // calculate deviation of trait from previous species trait and add new value
    double trait_dev = d.norm(d.rng) * p.trait_sigma(i);
    double parent_trait_val = d.indTraitL(dead_index);
    d.indTraitL(dead_index) = parent_trait_val + trait_dev;
    
    //if(print){ Rcout << "indSpeciesL" << "\n";}
    // the indv that was birthed or immigrated becomes the new species
    d.indSpeciesL(dead_index) = d.nTipsP(0); 
    
    //if(print){ Rcout << "spLastOriginStepL" << "\n";}
    // update time of last origin
    d.spLastOriginStepL(d.indSpeciesL(dead_index)) = i;
    
    //if(print){ Rcout << "founderFlagL" << "\n";}
    // set founderFlag of dead indv to 1 
    d.founderFlagL(dead_index) = 1; 
}

// call the phylo aspect of speciation
void update_speciation_phylo(int i, roleDataCpp &d, roleParamsCpp &p, int speciation_sp){
    
    // augment vectors by a set amount if the id of the new species (the n tips in the phylo)
    //  exceeds the size of the local sp vectors
    if(d.nTipsP(0) > d.aliveP.length() - 1){
        //if(print){ Rcout << "buffering in phylo" << "\n";}
        d.aliveP = buffer_logical_vector(d.aliveP,"FALSE",100); // buffer alive
        d.tipNamesP = buffer_character_vector(d.tipNamesP,"",100); // buffer tip names
        d.lengthsP = buffer_arma_vec(d.lengthsP,-2,100); // buffer edge lengths
        d.edgesP = buffer_arma_mat(d.edgesP,-2,100); // buffer edges
    }
    
    // set easy named variables for phylogeny speciation
    //arma::imat e = d.edgesP;
    //arma::vec l = d.lengthsP;
    int n = d.nTipsP(0);
    
    LogicalVector alive = d.aliveP;
    
    // nrows of the edge matrix
    //int eMax = e.nrow();
    int eMax = arma::size(d.edgesP)[0];
    
    // find index of where unrealized edges in edge matrix start
    // equivalent to eNew <- min(which(e[, 1] == -1))
    int eNew = -1;
    
    for (int k = 0; k < eMax; k++) {
        if (d.edgesP(k, 0) == -2) {
            eNew = k;
            break;
       }
    }
    
    // additional safety check
    if(eNew == -1){
        //if(print){ Rcout << "do safety buffer" << "\n";}
        
        //if(print){ Rcout << "lengthsP before" << "\n";}
        //if(print){ Rcout << d.lengthsP << "\n";}
        d.lengthsP = buffer_arma_vec(d.lengthsP,-2,100); // buffer edge lengths
        //if(print){ Rcout << "lengthsP after" << "\n";}
        //if(print){ Rcout << d.lengthsP << "\n";}
        
        //if(print){ Rcout << "edgesP before" << "\n";}
        //if(print){ Rcout << d.edgesP << "\n";}
        d.edgesP = buffer_arma_mat(d.edgesP,-2,100); // buffer edges
        //if(print){ Rcout << "edgesP after" << "\n";}
        //if(print){ Rcout << d.edgesP << "\n";}
        
        int eMax = arma::size(d.edgesP)[0];
        
        for (int k = 0; k < eMax; k++) {
            if (d.edgesP(k, 0) == -2) {
                eNew = k;
                break;
            }
        }
    }
    
    // index of the edge matrix of where to add new edge
    //arma::uvec inds = find(e.col(1) == i);
    //int j = inds(0);
    int j = -1;
    for (int k = 0; k < eMax; k++) {
       if (d.edgesP(k, 1) == speciation_sp) {
            j = k;
            break;
        }
    }
    
    // add one to internal node indices
    for (int r = 0; r < eNew; r++) {
        for (int c = 0; c < 2; c++){
            if (d.edgesP(r, c) >= n) {
                d.edgesP(r, c) ++;
            }
        }
    }
    
    //if(print){ Rcout << "edgesP" << "\n";}
    //if(print){ Rcout << "J" << j << "\n";}
    //if(print){ Rcout << "eNew" << eNew << "\n";}
    
    // add new internal node
    //int newNode = 2 * eMax + 1; // index of new node n+1
    int newNode = 2 * n; // this worked with prev approach
    d.edgesP(eNew, 0) = newNode;
    d.edgesP(1 + eNew, 0) = newNode;
    
    // add tips
    d.edgesP(eNew, 1) = d.edgesP(j, 1); // replace old tip
    d.edgesP(eNew + 1, 1) = n; // add new tip
    
    // update ancestry of internal nodes
    d.edgesP(j, 1) = newNode;
    
    //if(print){ Rcout << "lengthsP" << "\n";}
    
    // set new edge lengths to 0
    d.lengthsP[eNew] = 0;
    d.lengthsP[1 + eNew] = 0;
    
    // increase all tip edge lengths by 1 time step
    for (int r = 0; r <= eNew + 1; r++) {
        if (d.edgesP(r, 1) <= n + 1) { //n+1
            d.lengthsP(r) ++;
        }
    }
    
    //if(print){ Rcout << "alive" << "\n";}
    // update alive vector
    alive(n) = TRUE; //  - double check that this updates properly
    
    // increment nTipsP
    d.nTipsP(0) = d.nTipsP(0) + 1;
}

// call the phylo aspect of speciation
void update_speciation_phylo_new(int i, roleDataCpp &d, roleParamsCpp &p, int speciation_sp){
    
    // set easy named variables for phylogeny speciation
    arma::imat e = d.edgesP;
    arma::vec l = d.lengthsP;
    int n = d.nTipsP(0);
    LogicalVector alive = d.aliveP;
    
    // nrows of the edge matrix
    int eMax = arma::size(e)[0];
    
    // index of where unrealized edges in edge matrix start
    // an unrealized edge is the next un-utilized edge, used to avoid augmenting
    int eNew = 2* eMax - 2; // 2* eMax -2
    
    // find index of the edge matrix of where to add new edge
    arma::uvec inds = find(e.col(1) == i);
    int j = inds(0);
    
    // add one to internal nodes
    arma::uvec internalNode = find(e > eMax); // should it be > or >=??????
    e.elem(internalNode) += 1;
    
    //  - possible strategy for augmentation
    // add if statement to catch whether eNew >= eMax
    // if eNew >= eMax, augment e with addition eMax rows
    // else leave as is
    
    // add new internal node
    //int newNode = 2 * eMax + 1; // index of new node n+1
    int newNode = 2 * n; // this worked with prev approach
    e(eNew, 0) = newNode;
    e(1 + eNew, 0) = newNode;
    
    // add tips
    e(eNew, 1) = e(j, 1); // replace old tip
    e(eNew + 1, 1) = n; // add new tip
    
    // update ancestry of internal nodes
    e(j, 1) = newNode;
    
    // augment edge lengths
    l[eNew] = 0;
    l[1 + eNew] = 0;
    
    // increase all tip edge lengths by 1 time step
    l(find(e.col(1) <= eNew + 1)) += 1;
    
    // update alive vector
    alive(n) = TRUE; //  - double check that this updates properly
    
    // increment nTipsP
    d.nTipsP(0) = d.nTipsP(0) + 1;
}

// call speciation
void call_speciation(int i, int dead_index, roleDataCpp &d, roleParamsCpp &p, bool dispersed_this_iter, bool print){
    
    // computing probs of speciation as sum of weighted meta and local abundances
    // NumericVector probs = get_speciation_probs(i,d,p);
    
    // sample using probs (may be able to replace with sample_index_using_probs)
    // IntegerVector s = sample(d.spAbundM.length(), 1, false, probs);
    
    // make i from 0 to phylo.n - 1 (previously 1 to phylo.n)
    //int speciation_sp = s(0) - 1;
    int speciation_sp = d.indSpeciesL(dead_index);
    
    // update the local and meta traits, sp vector, and time of last origin
    update_speciation_local_meta(i,dead_index,d,p,dispersed_this_iter,speciation_sp,print);
    
    // update the phylogeny
    update_speciation_phylo(i,d,p,speciation_sp);
}

// update the squared differences between the traits of every individual
void update_trait_diffs_sq(int dead_index, roleDataCpp &d, roleParamsCpp &p){
    
    // update traitDiffsSq using two for loops 
    // updates ONLY the row and column of the dead_index
    for(int r = 0; r < p.individuals_local[0]; r++) 
    {
        d.traitDiffsSq(r,dead_index) = pow(d.indTraitL(r) - d.indTraitL(dead_index),2);
    }
    for(int c = 0; c < p.individuals_local[0]; c++)
    {
        d.traitDiffsSq(dead_index,c) = pow(d.indTraitL(dead_index) - d.indTraitL(c),2);
    }
}

// update a vector containing every local trait raised to the power of 2
// JACOB  - don't think we need this as we are updating the index specifically in update_env_filter_probs
void update_trait_pow(int dead_index, roleDataCpp &d, roleParamsCpp &p){
    
    // update traitPow using two for loops 
    // updates ONLY the row and column of the dead_index
    d.traitPow(dead_index) = pow(d.indTraitL(dead_index)-0,2);
}

// if non-neutral, save on computation by updating env_filtering probs fully if env_sigma changes and only for the new indv if it doesn't
// make sure env_filter_probs initial null state is used properly in death
void update_env_filter_probs(int i, int dead_index,roleDataCpp &d,roleParamsCpp &p){
    
    // check if new env_sigma, if so must recalculate envFilterProbs entirely
    if(d.prevEnvSigma != p.env_sigma[i]){
        d.envFilterProbs = 1 - exp((-1/p.env_sigma(i)) * pow(d.indTraitL - 0, 2));
    }
    //otherwise just update the probs of the new individual
    else{
        d.envFilterProbs(dead_index) = 1 - exp((-1/p.env_sigma(i)) * pow(d.indTraitL(dead_index), 2));
    }
    d.prevEnvSigma = p.env_sigma[i]; // NOTE - double check this is right
}

// update the local species sum of reciprocals
// used in R in the genetic simulation 
void update_local_sp_sum_recipr(int i, int dead_index, roleDataCpp &d, roleParamsCpp &p){

    // for each species...
    for(int s = 0; s < d.nTipsP(0); s++) {
        // if species is currently alive in local
        if(d.spAbundL(s) > 0){
            
            // add new abundance to species reciprocal sum 
            d.spReciprSumL(s) = d.spReciprSumL(s) + (1/d.spAbundL(s));
            
            // get n, the number of steps in this emergence period, as
            // the current iteration - the iteration of origin 
            int n = i - d.spLastOriginStepL(s);
            
            // harmonic mean is then n / the current reciprocal sum 
            d.spAbundHarmMeanL(s) = n / d.spReciprSumL(s);
        }
        // else species is dead in local
        else{
            d.spReciprSumL(s) = 0; 
            d.spAbundHarmMeanL(s) = 0;
        }
    }
}

// copy a roleDataCpp object to it's R side S4 equivalent
// also calculate equilibProp before the copy 
S4 role_data_from_cpp(roleDataCpp &d){
    
    // ANDY NOTE - clone each element of d or clone d and extract elements
    // JACOB NOTE - tried this, but pure C++ objects cannot be cloned
    //      trying to clone a roleDataCpp fails because it is not a matching R/Cpp type like a NumericVector
    //      instead tried, and none worked:
    // 1. copy constructor 
    // 2. copy assignment operator
    // 3. using * to get values of pointers rather than pointers themselves
    // All failed to deep copy, possibly some quirk of Rcpp where it would normally work in C++
    // take a bit of a closer look at clone
    
    // calculate equilib prop as the sum of the 0/1 numeric founderFlagL / J 
    double equilib_prop = sum(d.founderFlagL) / d.indSpeciesL.length();

    // construct S4 object, cloning each member of roleDataCpp directly to slots
    S4 out_l("localComm");
    out_l.slot("indSpecies") = Rcpp::clone(d.indSpeciesL);
    out_l.slot("indTrait") = Rcpp::clone(d.indTraitL);
    out_l.slot("spAbund") = Rcpp::clone(d.spAbundL);
    out_l.slot("spTrait") = Rcpp::clone(d.spTraitL);
    out_l.slot("spAbundHarmMean") = Rcpp::clone(d.spAbundHarmMeanL);
    out_l.slot("spLastOriginStep") = Rcpp::clone(d.spLastOriginStepL);
    out_l.slot("spExtinctionStep") = Rcpp::clone(d.spExtinctionStepL);
    out_l.slot("equilibProp") = equilib_prop;

    S4 out_m("metaComm");
    out_m.slot("spAbund") = Rcpp::clone(d.spAbundM);
    out_m.slot("spTrait") = Rcpp::clone(d.spTraitM);

    S4 out_p("rolePhylo");
    out_p.slot("n") = Rcpp::clone(d.nTipsP);
    out_p.slot("e") = Rcpp::clone(Rcpp::wrap(d.edgesP)); //Rcpp::as<T>(obj)
    out_p.slot("l") = Rcpp::clone(Rcpp::wrap(d.lengthsP));
    out_p.slot("alive") = Rcpp::clone(d.aliveP);
    out_p.slot("tipNames") = Rcpp::clone(d.tipNamesP);
    out_p.slot("scale") = Rcpp::clone(d.scaleP);

    S4 out_d("roleData");
    out_d.slot("localComm") = out_l;
    out_d.slot("metaComm") = out_m; 
    out_d.slot("phylo") = out_p;
    
    return(out_d);
}
//' @title iterModelCpp
//' @name iterModelCpp
//' @param local local
//' @param meta meta
//' @param phylo phylo
//' @param params params
//' @param print print
//
//
// [[Rcpp::export]]
List iterModelCpp(RObject local, RObject meta, RObject phylo, List params, bool print) {

    
    if(print){ Rcout << "save niter and niterTimestep"<< "\n";}
    
    // save niter and niterTimestep
    int niter = as<IntegerVector>(params["niter"])[0];
    int niter_timestep = as<IntegerVector>(params["niterTimestep"])[0];
    
    if(print){ Rcout << "make cpp objects for data and params"<< "\n";}
    // make cpp objects for data and params
    roleDataCpp d(local,meta,phylo);
    roleParamsCpp p = roleParamsCpp(params, niter); // constructor samples/stretches
    
    // setup rng distr using params
    // JACOB  - allow these to time vary - make new distr every iteration when any of these params change
    //  change to unorm and multiple by stdev
    //d.tnorm = std::normal_distribution<double>(0, p.trait_sigma[0] / (p.speciation_meta[0] + p.extinction_meta[0]));
    //d.sptnorm = std::normal_distribution<double>(0, p.trait_sigma[0]);
    
    // save n_indv to use easily throughout
    int n_indv = p.individuals_local[0];
    
    // if not neutral, calculate initial probs of death from environmental filtering and set to data
    NumericVector env_filter_probs (n_indv,1.0);
    if(p.neut_delta[0] != 1)
    {
        env_filter_probs = 1 - exp(-1/p.env_sigma[0] * pow(d.indTraitL - 0, 2));
    }
    d.envFilterProbs = env_filter_probs;
    
    // save the env_sigma of the first iter as the trivial prev for the first titer
    d.prevEnvSigma = p.env_sigma[0];
    
    // create out array to hold timeseries data 
    List out((niter / niter_timestep) + 1);
    
    // save the initial state to index 0 
    out[0] = role_data_from_cpp(d);
    
    if(print){ Rcout << "loop from 0 to niter - 1" << "\n";}
    // loop from 0 to niter - 1 
    for(int i = 0; i < niter; i++) {

        if(print){ Rcout << "call death" << "\n";}
        // METHOD - call death
        int dead_index = call_death(i,d,p);
        
        // set dispersal var for use in speciation, which is slightly different depending
        bool dispersed_this_iter = false; 
        
        if(print){ Rcout << "call birth or dispersal" << "\n";}
        // check for birth (prob = 1 - dispersal_prob) 
        if(d.unif(d.rng) >= p.dispersal_prob[i]){ //R::runif(0,1)

            // sample for the parent
            int parent_indv = sample_zero_to_x(p.individuals_local[i]);
            
            // METHOD - call birth
            call_birth(i, dead_index, parent_indv, d, p, print);
        }
        else{ // ROCKS NOTE - add new param for birth_prob where b + d <= 1, then dispersal/birth are allowed to NOT happen, leaving a rock
            dispersed_this_iter = true;
            
            // sample a parent index using species abundances as probs
            int parent_index = sample_index_using_probs(d.spAbundM);
            
            // METHOD - call dispersal
            call_dispersal(i,dead_index, parent_index, d, p, print);
        }
        // ROCKS NOTE - else add a rock to the dead_index and prevent speciation from occurring
        
        if(print){ Rcout << "call speciation if it occurs" << "\n";}
        // randomly decide if speciation occurs
        if(d.unif(d.rng) < p.speciation_local[i])
        {
            // METHOD - call speciation and get the chosen species 
            call_speciation(i, dead_index, d, p, dispersed_this_iter, print);
        }
        
        if(print){ Rcout << "if non-neutral, update objects used for comp and filtering" << "\n";}
        // if non-neutral, update objects used for comp and filtering
        if(p.neut_delta[0] < 1){
            
            // JACOB NOTE - changed this as part of new death ifs 
            // if filtering only
            if(p.env_comp_delta[0] == 1){
                // METHOD - update the probs of environmental filtering ONLY if non-neutral and either for all indv or just the new one
                update_env_filter_probs(i, dead_index,d,p);
            }
            else if(p.env_comp_delta[0] == 0){ // if competition only 
                // METHOD - update the squared differences between traits at the dead_index
                update_trait_diffs_sq(dead_index,d,p);
            }
            else{
                // update both
                update_env_filter_probs(i, dead_index,d,p);
                update_trait_diffs_sq(dead_index,d,p);
            }
        }
        
        // METHOD - update the local species sum of reciprocals
        update_local_sp_sum_recipr(i, dead_index, d, p);
        
        // save if i is 0, 9, 19 ... 99 
        if((i + 1) % niter_timestep == 0)
        {
            // if iter is 0, add to 0th index
            if(i == 0){
                // METHOD - convert to R S4 object
                out[i] = role_data_from_cpp(d);
            }
            // otherwise add to the correct timestep index
            else{
                out[(i + 1) / niter_timestep] = role_data_from_cpp(d); 
            }
        }
        
        // JACOB NOTE - if a certain amount of equilibrium achieved, return(out)
    } 
    
    return(out);
};

// 
// // these funs, one per return data type, are R wrappers around multiple Cpp functions
// // these avoid having to wrap all ~20 functions individually 
// // ONLY used for testing and nothing else
// //' @name intFunCpp
// //' @title intFunCpp
// //' @param fun_Name fun_Name
// //' @param probs probs
// //' @param x x
// //
// // [[Rcpp::export]]
// int intFunCpp(Rcpp::StringVector fun_name,
//                 NumericVector probs=NULL, 
//                 int x=NULL) {
//     std::string fn = Rcpp::as<std::string>(fun_name(0));
//    
//     
//     // tried switch, didn't work but may revisit
//     if(fn == "sample_index_using_probs"){
//         return(sample_index_using_probs(probs));
//     }
//     if(fn == "sample_zero_to_x"){
//         return(sample_zero_to_x(x));
//     }
// }
// 
// 
// // need to update this so it can deal with params as a List
// 
// 
// //' @title dataFunCpp
// //' @name dataFunCpp
// //' @param fun_name fun_name
// //' @param local local
// //' @param meta meta
// //' @param phylo phylo
// //' @param params params
// //' @param niter niter
// //' @param i i
// //' @param dead_index dead_index
// //' @param parent_indv parent_indv
// //' @param dispersed_this_year dispersed_this_year
// //' @param speciation_sp speciation_sp
// //
// // [[Rcpp::export]]
// S4 dataFunCpp(Rcpp::StringVector fun_name,
//               RObject local=NULL, RObject meta=NULL,RObject phylo=NULL, //used universally
//               List params=NULL, int niter=NULL, int i=NULL, //used universally
//               int dead_index=NULL, // used universally
//               int parent_indv=NULL, // used by call_birth and call_dispersal
//               bool dispersed_this_iter=NULL, // used by call_speciation and update_speciation_local_meta
//               int speciation_sp=NULL) { // used in update_speciation_local_meta
// 
//     // create Cpp objects
//     roleDataCpp d(local,meta,phylo);
//     roleParamsCpp p(params,niter);
//     std::string fn = Rcpp::as<std::string>(fun_name(0));
// 
//     // tried switch, didn't work
//     if(fn == std::string("call_birth")){
//         call_birth(i,dead_index,parent_indv,d, p,true);
//         return(role_data_from_cpp(d));
//     }
//     if(fn == std::string("call_dispersal")){
//         call_dispersal(i,dead_index,parent_indv,d, p,true);
//         return(role_data_from_cpp(d));
//     }
//     if(fn == std::string("update_trait_diffs_sq")){
//         update_trait_diffs_sq(dead_index,d, p);
//         return(role_data_from_cpp(d));
//     }
//     if(fn == std::string("call_speciation")){
//         call_speciation(i, dead_index, d, p, dispersed_this_iter, true);
//         return(role_data_from_cpp(d));
//     }
//     if(fn == std::string("update_speciation_local_meta")){
//         update_speciation_local_meta(i,dead_index,d,p,dispersed_this_iter,speciation_sp,true);
//         return(role_data_from_cpp(d));
//     }
//     if(fn == std::string("update_speciation_phylo")){
//         update_speciation_phylo(i,d,p,speciation_sp);
//         return(role_data_from_cpp(d));
//     }
// }
// //' @title vectFunCpp
// //' @name vectFunCpp
// //' @param fun_name fun_name
// //' @param local local
// //' @param meta meta
// //' @param phylo phylo
// //' @param params params
// //' @param niter niter
// //' @param i i 
// // [[Rcpp::export]]
// NumericVector vectFunCpp(Rcpp::StringVector fun_name,
//                          RObject local=NULL, RObject meta=NULL,RObject phylo=NULL, // used universally 
//                          List params=NULL, int niter=NULL, int i = NULL){ // used universally 
//     // create Cpp objects
//     roleDataCpp d(local,meta,phylo);
//     roleParamsCpp p(params,niter);
//     std::string fn = Rcpp::as<std::string >(fun_name(0));
//     
//     // tried switch, didn't work
//     if(fn == "get_filtering_death_probs"){
//         return(get_filtering_death_probs(i,d,p));
//     }
//     if(fn == "get_comp_death_probs"){
//         return(get_comp_death_probs(i,d,p));
//     }
//     if(fn == "get_speciation_probs"){
//         return(get_speciation_probs(i,d,p));
//     }
// }
// 
// // export Rcpp modules for use in R
// // only 4 functions get exported - the core loop and the test wrappers
// RCPP_MODULE(iterModelCpp) {
//     function("iterModelCpp", &iterModelCpp);
//     function("intFunCpp", &intFunCpp);
//     function("vectFunCpp", &vectFunCpp);
//     // function("dataFunCpp", &dataFunCpp);
// }
