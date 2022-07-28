// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

#include "roleDataCpp.cpp"
#include "roleParamsCpp.cpp"

using namespace Rcpp;

// still missing trait_diffs for fast comp filtering computation 
// prob of selecting a parent for speciation is no longer metacomm abundance weighted by dispersal prob + local 
// death probs due to filtering/ comp
// replicates sample(1:x, 1)
int sample_zero_to_x(int x)
{
    return((int) (R::runif(0,1) * (double) x));
}

NumericVector get_zero_to_x_vector(int x)
{
    std::vector<int> v(x);
    std::iota(v.begin(), v.end(), 1);
    return(wrap(v));
}

int sample_index_using_probs(NumericVector probs){
    IntegerVector v = Rcpp::sample(probs.length(), 1, false, probs, false);
    return(v(0));
}

List iterModelCpp(RObject local, RObject meta, RObject phylo, RObject params, bool print) {
    if(print){Rcout << "iter loop started" << "\n";}
    
    // save niter and niterTimestep
    int niter = params.slot("niter");
    int niter_timestep = params.slot("niterTimestep");

    // make cpp objects for data and params
    roleDataCpp d(local,meta,phylo);
    roleParamsCpp p = roleParamsCpp(params,niter); // constructor samples/stretches
    
    int n_indv = p.individuals_local[0];
    
    // calculate initial probs of death from environmental filtering
    NumericVector env_filter_probs = 1 - exp(-1/p.env_sigma[0] * pow(d.indTraitL - 0, 2));
    double prev_env_sigma = p.env_sigma[0];
    
    // create out array to hold timeseries data 
    RObject out[(niter / niter_timestep) + 1]; 
    
    // loop from 0 to niter - 1 
    for(int i = 0; i < (int) params.slot("niter"); i++) {
        if(print){Rcout << "started iteration " << i << "\n";}
            
        // // sample for dead index 
        // calculate probs of death due to env filtering
        // picture curve with prob of death being 0 at trait optimum
        // curve equation is 1-exp(-1(xi - xopt)^2 / env_sigma
        // env_sigma determines if curve is wide or narrow
        // traits optimum is 0 
        
        // note that env_filter_probs is calculated before the loop and updated when indTraitL or env_sigma changes
        
        // sum of
        // sigma comp makes the curve wider, making more dissimilar individuals having greater impact on death
        // lower sigma comp means individual is only competing with itself and very close indv is trait space 
        // recall that we called the optimum trait_z for the use case where trait optima changes as the environment changes 
        //NumericVector c_probs = Rcpp::as<NumericVector>(Rcpp::wrap(1/p.individuals_local * sum(exp((-1/p.comp_sigma(i)) * d.traitDiffsSq,0)));
        //Rcout << "trait diffs sq " << d.traitDiffsSq << "\n";
        //Rcout << "arma exp " << arma::exp((-1/p.comp_sigma(i)) * d.traitDiffsSq) << "\n";
        //Rcout << "arma sum " << arma::sum(arma::exp((-1/p.comp_sigma(i)) * d.traitDiffsSq),1) << "\n";
        //Rcout << "type of arma sum" << typeid(arma::sum(arma::exp((-1/p.comp_sigma(i)) * d.traitDiffsSq),1)).name() << "\n";
        //Rcout << "arma sum / indv" << arma::sum(arma::exp((-1/p.comp_sigma(i)) * d.traitDiffsSq),1) / n_indv << "\n";
        
        arma::colvec comp_probs_a = arma::sum(arma::exp((-1/p.comp_sigma(i)) * d.traitDiffsSq),1) / n_indv;
        NumericVector comp_probs = Rcpp::wrap(comp_probs_a.begin(),comp_probs_a.end());
                                                               
        //Rcout << "comp probs " << comp_probs << "\n";
        //Rcout << "env filter probs " << "\n" << env_filter_probs << "\n";
        // 1/indv_local averages
        // save squares of traitdiffs instead of traitdiffs
        
        // calculate combined death probs
        NumericVector death_probs = env_filter_probs + comp_probs;
        //Rcout << "dead probs " << death_probs << "\n";
        
        // get dead index and dead species
        int dead_index = sample_index_using_probs(death_probs);
        if(print){Rcout << "dead_index: " << dead_index << "\n";}
        int dead_species = d.indSpeciesL(dead_index);
        if(print){Rcout << "dead_species: " << dead_species << "\n";}
        
        // subtract from abundance of dead species
        d.spAbundL(dead_species) = d.spAbundL(dead_species) - 1;
        
        // check for extinction
        // extinct if species abundance of dead individual is now <0 
        bool extinct_in_local = d.spAbundL(dead_species) <= 0;
        // not in meta if species of the dead individual is not in meta
        bool not_in_meta = dead_species > d.spAbundM.length();
        // check for extinction
        if(extinct_in_local & not_in_meta){
            if(print){Rcout << "extinction occured in local, species: " << dead_species << "\n";}
            d.aliveP(dead_species) = false;
        }
        
        // check for birth (prob = 1 - dispersal_prob) 
        if(R::runif(0,1) >= p.dispersal_prob(i)){
            if(print){Rcout << "birthed, dispersal prob: " << p.dispersal_prob(i) << "\n";}
            
            // sample for the parent
            int parent_indv = sample_zero_to_x(p.individuals_local(i));
            if(print){Rcout << "parent indv: " << parent_indv << "\n";}
            
            //int birthed_species = d.indSpTrtL(parent_indv,0); 
            int birthed_species = d.indSpeciesL(parent_indv); 
            
            // set the species of the new indv to that of the parent 
            //d.indSpTrtL(dead_index,0) = birthed_species;
            d.indSpeciesL(dead_index) = birthed_species;
            
            //if(print){Rcout << "set new indv to parent species: " << d.indSpTrtL(dead_index,0) << " from" << d.indSpTrtL(parent_indv,0)<< "\n";}
            
            // add to the abundance of the species matrix
            d.spAbundL(birthed_species) = d.spAbundL(birthed_species) + 1; 
            
            double var = p.trait_sigma(1) / (p.speciation_meta(1) + p.extinction_meta(1));
                                                 
            // calculate trait change from parent
            NumericVector trait_change = Rcpp::rnorm(1, 0, p.trait_sigma(1) / (p.speciation_meta(1) + p.extinction_meta(1)));
            if(print){Rcout << "trait change from parent: " << trait_change << "\n";}
            
            // add new trait value
            d.indTraitL(dead_index) = d.indTraitL(parent_indv) + trait_change(0);
            if(print){Rcout << "new trait value: " <<  d.indTraitL(dead_index) << "\n";}
            
            // update env_filter_probs
            
            // update traitDiffsSq by changing the row and column of the new individual
        }
        else{
            if(print){Rcout << "dispersed, dispersal prob: " <<  p.dispersal_prob(i) << "\n";}
    
            // sample a parent index using species abundances as probs
            int parent_index = sample_index_using_probs(d.spAbundM);
            if(print){Rcout << "sampled parent index:" <<  parent_index << "\n";}
            
            // set the species to the parent species from meta
            d.indSpeciesL(dead_index) = parent_index;
            
            // update trait    
            // calculate trait change from parent
            NumericVector trait_change = Rcpp::rnorm(1, 0, p.trait_sigma(1) / (p.speciation_meta(1) + p.extinction_meta(1)));
            if(print){Rcout << "trait change from parent: " << trait_change << "\n";}
            
            // add new trait value
            d.indTraitL(dead_index) = d.spTraitM(parent_index) + trait_change(0);
        }
        
        // update traitDiffsSq using two for loops 
        for(int r = 0; r < n_indv; r++)
        {
            d.traitDiffsSq(r,dead_index) = pow(d.indTraitL(r) - d.indTraitL(dead_index),2);
        }
        for(int c = 0; c < n_indv; c++)
        {
            d.traitDiffsSq(dead_index,c) = pow(d.indTraitL(dead_index) - d.indTraitL(c),2);
        }
        
        // if new env_sigma, must recalculate entirely
        if(prev_env_sigma != p.env_sigma[i]){
            env_filter_probs = 1 - exp((-1/p.env_sigma(i)) * pow(d.indTraitL - 0, 2));
        }
        //otherwise just update the probs of the new individual
        else{
            env_filter_probs(dead_index) = 1 - exp((-1/p.env_sigma(i)) * pow(d.indTraitL(dead_index), 2));
        }
        
        // if speciation occurs
        if(R::runif(0,1) < p.speciation_local(i))
        {
            if(print){Rcout << "speciated, prob:" <<  p.speciation_local(i) << "\n";}
            //d.indSpTrtL(dead_index,0) = d.nTipsP + 1;
            //d.nTipsP = d.nTipsP + 1;
        }
        
        // save if i is 0, 9, 19 ... 99 
        if(i == 0 || (i + 1) % niter_timestep == 0) //i == 0 || (i + 1) % niter_timestep == 0
        {
            
            if(print){Rcout << "saving, niter: " << i << "\n";}
            
            S4 out_l("localComm");
            //out_l.slot("indSppTrt") = Rcpp::clone(d.indSpTrtL);
            out_l.slot("indSpecies") = Rcpp::clone(d.indSpeciesL);
            out_l.slot("indTrait") = Rcpp::clone(d.indTraitL);
            
            if(print){Rcout << "created l" << i << "\n";}
            
            S4 out_m("metaComm");
            //out_m.slot("sppAbundTrt") = Rcpp::clone(d.spAbundTrtM);
            out_m.slot("spAbund") = Rcpp::clone(d.spAbundM);
            out_l.slot("spTrait") = Rcpp::clone(d.spTraitM);
            
            if(print){Rcout << "created m" << "\n";}
            
            S4 out_p("rolePhylo");
            out_p.slot("n") = Rcpp::clone(d.nTipsP);
            out_p.slot("e") = Rcpp::clone(d.edgesP);
            out_p.slot("l") = Rcpp::clone(d.lengthsP);
            out_p.slot("alive") = Rcpp::clone(d.aliveP);
            out_p.slot("tipNames") = Rcpp::clone(d.tipNamesP);
            out_p.slot("scale") = Rcpp::clone(d.scaleP);
            if(print){Rcout << "created p" << "\n";}
            
            S4 out_d("roleData");
            out_d.slot("localComm") = out_l;
            out_d.slot("metaComm") = out_m; 
            out_d.slot("phylo") = out_p;
            
            if(print){Rcout << "created new data" << "\n";}
            
            // if iter is 0, add to 0th index
            if(i == 0){
                out[i] = out_d; 
            }
            // otherwise add to the correct timestep index
            else{
                out[(i + 1) / niter_timestep] = out_d; 
            }
        
            if(print){Rcout << "added to out" << "\n";}
            
        }
    } 
    
    List out_list = List::create(out[0]);
    for(int j = 1; j < sizeof(out)/sizeof(out[0]); j++)
    {
        out_list.push_back(out[j]);
    }
    return(out_list);
};

RCPP_MODULE(iterModelCpp) {
    function("iterModelCpp", &iterModelCpp);
}