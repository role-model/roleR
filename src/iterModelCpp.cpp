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


// [[Rcpp::export]]

List iterModelCpp(RObject local, RObject meta, RObject phylo, RObject params, bool print) {
    if(print){Rcout << "iter loop started" << "\n";}
    
    // save niter and niterTimestep
    int niter = params.slot("niter");
    int niter_timestep = params.slot("niterTimestep");

    // make cpp objects for data and params
    roleDataCpp d(local,meta,phylo);
    roleParamsCpp p = roleParamsCpp(params,niter); // constructor samples/stretches
    
    // save n_indv to use easily throughout
    int n_indv = p.individuals_local[0];
    
    // if not neutral, calculate initial probs of death from environmental filtering
    NumericVector env_filter_probs (1.0,n_indv);
    if(p.neut_delta[0] != 1)
    {
        env_filter_probs = 1 - exp(-1/p.env_sigma[0] * pow(d.indTraitL - 0, 2));
        if(print){Rcout << "calculated initial probs of death from env filtering: " << env_filter_probs << "\n";}
    }
    double prev_env_sigma = p.env_sigma[0];
    
    // create out array to hold timeseries data 
    RObject out[(niter / niter_timestep) + 1]; 
    
    // save niter 0 
    if(print){Rcout << "saving, niter: 0" << "\n";}
    S4 out_l("localComm");
    //out_l.slot("indSppTrt") = Rcpp::clone(d.indSpTrtL);
    out_l.slot("indSpecies") = Rcpp::clone(d.indSpeciesL);
    out_l.slot("indTrait") = Rcpp::clone(d.indTraitL);
    out_l.slot("spAbund") = Rcpp::clone(d.spAbundL);
    S4 out_m("metaComm");
    //out_m.slot("sppAbundTrt") = Rcpp::clone(d.spAbundTrtM);
    out_m.slot("spAbund") = Rcpp::clone(d.spAbundM);
    out_l.slot("spTrait") = Rcpp::clone(d.spTraitM);
    S4 out_p("rolePhylo");
    out_p.slot("n") = Rcpp::clone(d.nTipsP);
    out_p.slot("e") = Rcpp::clone(d.edgesP);
    out_p.slot("l") = Rcpp::clone(d.lengthsP);
    out_p.slot("alive") = Rcpp::clone(d.aliveP);
    out_p.slot("tipNames") = Rcpp::clone(d.tipNamesP);
    out_p.slot("scale") = Rcpp::clone(d.scaleP);
    S4 out_d("roleData");
    out_d.slot("localComm") = out_l;
    out_d.slot("metaComm") = out_m; 
    out_d.slot("phylo") = out_p;
    out[0] = out_d; 
    
    // loop from 0 to niter - 1 
    for(int i = 0; i < (int) params.slot("niter"); i++) {
        if(print){Rcout << "started iteration " << i << "\n";}
            
        // NOTES
        // sample for dead index
        // picture curve with prob of death being 0 at trait optimum
        // curve equation is 1-exp(-1(xi - xopt)^2 / env_sigma
        // env_sigma determines if curve is wide or narrow
        // traits optimum is 0
        // env_filter_probs is calculated before the loop and updated when indTraitL or env_sigma changes
        // sigma comp makes the curve wider, making more dissimilar individuals having greater impact on death
        // lower sigma comp means individual is only competing with itself and very close indv is trait space
        // recall that we called the optimum trait_z for the use case where trait optima changes as the environment changes
        // either find param values that make things neutral, or togglable parameter that makes things neutral
        
        // neutral death probs is rep(1,n_indv)
        NumericVector death_probs(1.0,n_indv);
        
        // if not neutral, calculate initial probs of death from environmental filtering
        if(p.neut_delta[0] != 1)
        {
            NumericVector env_filter_probs = 1 - exp(-1/p.env_sigma[0] * pow(d.indTraitL - 0, 2));
            if(print){Rcout << "calculated initial probs of death from env filtering: " << env_filter_probs << "\n";}
            double prev_env_sigma = p.env_sigma[0];

            arma::colvec comp_probs_a = arma::sum(arma::exp((-1/p.comp_sigma(i)) * d.traitDiffsSq),1) / n_indv;
            NumericVector comp_probs = Rcpp::wrap(comp_probs_a.begin(),comp_probs_a.end());
            if(print){Rcout << "calculated comp probs: " << comp_probs << "\n";}
        
            // calculate combined death probs
            death_probs = p.neut_delta[0] + (1-p.neut_delta[0]) * (env_filter_probs + comp_probs);
            //death_probs = env_filter_probs + comp_probs;
            if(print){Rcout << "combined env and comp probs: " << death_probs << "\n";}
            //Rcout << "dead probs " << death_probs << "\n";
        }
        
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
            if(print){Rcout << "extinction occured in local only species, species: " << dead_species << "\n";}
            d.aliveP(dead_species) = false;
        }
        
        // set dispersal var for use in speciation, which is slightly different depending
        bool dispersed_this_iter = false; 
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
                                                 
            // calculate trait change from parent
            NumericVector trait_change = Rcpp::rnorm(1, 0, p.trait_sigma(1) / (p.speciation_meta(1) + p.extinction_meta(1)));
            if(print){Rcout << "trait change from parent: " << trait_change << "\n";}
            
            // add new trait value
            d.indTraitL(dead_index) = d.indTraitL(parent_indv) + trait_change(0);
            if(print){Rcout << "new trait value: " <<  d.indTraitL(dead_index) << "\n";}
        }
        else{
            if(print){Rcout << "dispersed, dispersal prob: " <<  p.dispersal_prob(i) << "\n";}
            dispersed_this_iter = true;
            
            // sample a parent index using species abundances as probs
            int parent_index = sample_index_using_probs(d.spAbundM);
            if(print){Rcout << "sampled parent index:" <<  parent_index << "\n";}
            
            // set the species to the parent species from meta
            d.indSpeciesL(dead_index) = parent_index;
            // add to the abundance of the species matrix
            d.spAbundL(d.indSpeciesL(dead_index)) = d.spAbundL(d.indSpeciesL(dead_index)) + 1;

            // calculate trait change from parent
            NumericVector trait_change = Rcpp::rnorm(1, 0, p.trait_sigma(1) / (p.speciation_meta(1) + p.extinction_meta(1)));
            if(print){Rcout << "trait change from parent: " << trait_change << "\n";}
            
            // add new trait value
            d.indTraitL(dead_index) = d.spTraitM(parent_index) + trait_change(0);
        }
        
        if(print){Rcout << "updating trait diffs sq" << "\n";}
        // update traitDiffsSq using two for loops 
        for(int r = 0; r < n_indv; r++)
        {
            d.traitDiffsSq(r,dead_index) = pow(d.indTraitL(r) - d.indTraitL(dead_index),2);
        }
        for(int c = 0; c < n_indv; c++)
        {
            d.traitDiffsSq(dead_index,c) = pow(d.indTraitL(dead_index) - d.indTraitL(c),2);
        }
        
        if(print){Rcout << "checking if new env sigma" << "\n";}
        
        if(p.neut_delta[0] != 1)
        {
            // if new env_sigma, must recalculate entirely
            if(prev_env_sigma != p.env_sigma[i]){
                if(print){Rcout << "recalculating full env sigma" << "\n";}
                env_filter_probs = 1 - exp((-1/p.env_sigma(i)) * pow(d.indTraitL - 0, 2));
            }
            //otherwise just update the probs of the new individual
            else{
                if(print){Rcout << "recalculating env sigma for new individual" << "\n";}
                env_filter_probs(dead_index) = 1 - exp((-1/p.env_sigma(i)) * pow(d.indTraitL(dead_index), 2));
            }
            prev_env_sigma = p.env_sigma[i];
        }
        
        // if speciation occurs
        if(R::runif(0,1) < p.speciation_local(i))
        {
            if(print){Rcout << "speciated, prob:" <<  p.speciation_local(i) << "\n";}
            
            float dp = p.dispersal_prob(i);
            if(print){Rcout << "dp" << dp << "\n";}
            if(print){Rcout << "spAbundM" << d.spAbundM << "\n";}
            if(print){Rcout << "spAbundL" << d.spAbundL << "\n";}
            if(print){Rcout << "computing probs of speciation as sum of weighted meta and local abundances" << "\n";}
            NumericVector probs = dp * (d.spAbundM / sum(d.spAbundM)) + 
                (1 - dp) * (d.spAbundL / sum(d.spAbundL));
            if(print){Rcout << "computed probs" << probs << "\n";}
            probs = (abs(probs)+probs)/2;
            IntegerVector s = sample(d.spAbundM.length(), 1, false, probs);
            // make i from 0 to phylo.n - 1 (previously 1 to phylo.n)
            int chosen_species = s(0) - 1;
            if(print){Rcout << "chosen speciation species: " << chosen_species << "\n";}
            
            //calculate deviation of trait from previous species trait
            float trait_dev = R::rnorm(0, p.trait_sigma(i));
            if(print){Rcout << "calculated trait dev: " << trait_dev << "\n";}
            float parent_trait_val = 0;
            if(dispersed_this_iter){
                if(print){Rcout << "getting parent trait val from meta" << "\n";}
                parent_trait_val = d.spTraitM(chosen_species);
            }
            else{
                if(print){Rcout << "getting parent trait val from meta" << "\n";}
                parent_trait_val = d.spTraitM(chosen_species);
                //parent_trait_val = d.indTraitL(chosen_species);
            }
            d.indTraitL(dead_index) = parent_trait_val + trait_dev;
            if(print){Rcout << "set new trait val";}
            
            // the indv that was birthed or immigrated becomes the new species
            d.indSpeciesL(dead_index) = d.nTipsP(0); 
            
            if(print){Rcout << "starting phylo element of speciation" << "\n";}
            bool print_matrices = false;
            
            // set easy named variables for phylogeny speciation
            NumericMatrix e = d.edgesP;
            //e = e(Rcpp::Range(0,110),_);
            
            NumericVector l = d.lengthsP;
            int n = d.nTipsP(0);
            if(print){Rcout << "n tips p: " << n << "\n";}
            int i = chosen_species;
            if(print){Rcout << "chosen species: " << i << "\n";}
            LogicalVector alive = d.aliveP;
            
            if(print_matrices){Rcout << "edge matrix at start: " << e << "\n";}
            
            // nrows of the edge matrix
            int eMax = e.nrow();
            if(print){Rcout << "got number of rows of edge matrix: " << eMax << "\n";}
            
            // find index of where unrealized edges in edge matrix start
            // eNew <- min(which(e[, 1] == -1))
            int eNew = -1;
            
            for (int k = 0; k < eMax; k++) {
                if (e(k, 0) == -2) {
                    eNew = k;
                    break;
                }
            }
            
            // add if statement to catch whether eNew >= eMax
            // if eNew >= eMax, augment e with addition eMax rows
            // else leave as is
            
            if(print){Rcout << "got index of where unrealized edges in edge matrix start: " << eNew << "\n";}
            
            // index of the edge matrix of where to add new edge
            // j <- which(e[, 2] == i)
            int j = -1;
            for (int k = 0; k < eMax; k++) {
                if (e(k, 1) == i) {
                    j = k;
                    break;
                }
            }
            
            if(print){Rcout << "got index of the edge matrix of where to add new edge: " << j << "\n";}
            // add one to internal node indices
            //e[e > n] <- e[e > n] + 1
            
            for (int r = 0; r < eNew; r++) {
                for (int c = 0; c < 2; c++){
                    if (e(r, c) >= n) {
                        e(r, c) ++;
                    }
                }
            }
            if(print_matrices){Rcout << "edge matrix after adding to internal node indices: " << e << "\n";}
            
            // add new internal node
            int newNode = 2 * n; // index of new node n+1
            if(print){Rcout << "got index of new internal node: " << newNode << "\n";}
            
            e(eNew, 0) = newNode;
            e(1 + eNew, 0) = newNode; // do this more elegantly
            if(print_matrices){Rcout << "edge matrix after adding new parent nodes: " << e << "\n";}
            
            // add tips
            e(eNew, 1) = e(j, 1); // add old tip
            e(eNew + 1, 1) = n; // add new tip
            if(print_matrices){Rcout << "edge matrix after adding new child nodes: " << e << "\n";}
            
            // update ancestry of internal nodes
            e(j, 1) = newNode;
            if(print_matrices){Rcout << "edge matrix after updating internal node ancestry: " << e << "\n";}
            
            // augment edge lengths
            l[eNew] = 0;
            l[1 + eNew] = 0;
            
            // e[Rcpp::range(0, eNew), ]
            // Rcout << "saving, niter: " << e(0:eNew, ) << "\n";
            
            // increase all tip edge lengths by 1 time step
            for (int r = 0; r <= eNew + 1; r++) {
                if (e(r, 1) <= n + 1) { //n+1
                    l(r) ++;
                }
            }
            
            // update alive vector
            alive(n) = TRUE; // double check this

            // increment nTipsP
            d.nTipsP(0) = d.nTipsP(0) + 1;
        }
        
        // save if i is 0, 9, 19 ... 99 
        if((i + 1) % niter_timestep == 0) //i == 0 || (i + 1) % niter_timestep == 0
        {
            
            if(print){Rcout << "saving, niter: " << i << "\n";}
            
            S4 out_l("localComm");
            //out_l.slot("indSppTrt") = Rcpp::clone(d.indSpTrtL);
            out_l.slot("indSpecies") = Rcpp::clone(d.indSpeciesL);
            out_l.slot("indTrait") = Rcpp::clone(d.indTraitL);
            out_l.slot("spAbund") = Rcpp::clone(d.spAbundL);
            if(print){Rcout << "created l" << "\n";}
            
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