rocks changes
split birth disp param into two
if neither birth or dispersal happens, replace the dead indv with a rock
rock properties
	species id is -1 or indv vector value set to 0 
	dont have traits

matrix of species effects on one another
desnsity dependent term for inter and intra 
	these must be independent 
	because mess is based on trait diffs these terms are always symmetric 

how do growth rates, carrying capacities, and  get picked for a new species from speciation? 
delta (prop death due to chance) and gamma (relative comp vs filtering) are already in place (???) 

basically, if lv params are specified change the death probs completely (?) 

iterModel(roleModel model, roleParams p)
	for i in 0:niter
		## calculate death probs
		# probs from env filtering
		# probs from competition
		
		## replace individual chosen for death using birth OR dispersal
		
		## speciate if it happens
		

iterModel(roleModel model, roleParams p) 
    
    // default is neutrality with equal probs of death due to env filtering
    env_filter_probs = rep(1,nindv)
	
	// if not neutral, set initial probs of death from environmental filtering
    if p.neut_delta[0] != 1
    {
        env_filter_probs = 1 - exp(-1/p.env_sigma[0] * pow(d.indTraitL - 0, 2));
    }
    prev_env_sigma = p.env_sigma[0];
    
    // loop from 0 to niter - 1 
    for i in 0:niter

        death_probs =  rep(1,nindv)
        
        // if not neutral, calculate initial probs of death from environmental filtering
        if p.neut_delta[0] != 1
            env_filter_probs = 1 - exp(-1/p.env_sigma[0] * pow(d.indTraitL - 0, 2));
         
            prev_env_sigma = p.env_sigma[0];

            comp_probs_a = arma::sum(arma::exp((-1/p.comp_sigma(i)) * d.traitDiffsSq),1) / n_indv;
            comp_probs = Rcpp::wrap(comp_probs_a.begin(),comp_probs_a.end());
        
            // calculate combined death probs
            death_probs = p.neut_delta[0] + (1-p.neut_delta[0]) * ((p.env_comp_delta[0] * env_filter_probs) + ((1-p.env_comp_delta[0]) * comp_probs));

        
        // get dead index and dead species
        int dead_index = sample_index_using_probs(death_probs);
        
        int dead_species = d.indSpeciesL(dead_index);
       
        // subtract from abundance of dead species
        d.spAbundL(dead_species) = d.spAbundL(dead_species) - 1;
        
        // check for extinction
        // extinct if species abundance of dead individual is now <0 
        extinct_in_local = d.spAbundL(dead_species) <= 0;
		
        // not in meta if species of the dead individual is not in meta
        not_in_meta = dead_species > d.spAbundM.length();
		
        // check for extinction
        if(extinct_in_local & not_in_meta){
            if(print){Rcout << "extinction occured in local only species, species: " << dead_species << "\n";}
            d.aliveP(dead_species) = false;
        }
        
        // set dispersal var for use in speciation, which is slightly different depending
        dispersed_this_iter = false; 
		
        // check for birth (prob = 1 - dispersal_prob) 
        if(R::runif(0,1) >= p.dispersal_prob(i)){
            
            // sample for the parent
            int parent_indv = sample_zero_to_x(p.individuals_local(i));
            
            //int birthed_species = d.indSpTrtL(parent_indv,0); 
            int birthed_species = d.indSpeciesL(parent_indv); 
            
            // set the species of the new indv to that of the parent 
            d.indSpeciesL(dead_index) = birthed_species;
                        
            // add to the abundance of the species matrix
            d.spAbundL(birthed_species) = d.spAbundL(birthed_species) + 1; 
                                                 
            // calculate trait change from parent
            NumericVector trait_change = Rcpp::rnorm(1, 0, p.trait_sigma(1) / (p.speciation_meta(1) + p.extinction_meta(1)));
            
            // add new trait value
            d.indTraitL(dead_index) = d.indTraitL(parent_indv) + trait_change(0);
        }
        else{
            dispersed_this_iter = true;
            
            // sample a parent index using species abundances as probs
            int parent_index = sample_index_using_probs(d.spAbundM);
            
            // set the species to the parent species from meta
            d.indSpeciesL(dead_index) = parent_index;
            // add to the abundance of the species matrix
            d.spAbundL(d.indSpeciesL(dead_index)) = d.spAbundL(d.indSpeciesL(dead_index)) + 1;

            // calculate trait change from parent
            NumericVector trait_change = Rcpp::rnorm(1, 0, p.trait_sigma(1) / (p.speciation_meta(1) + p.extinction_meta(1)));
            
            // add new trait value
            d.indTraitL(dead_index) = d.spTraitM(parent_index) + trait_change(0);
        }
        
        // update traitDiffsSq at the new individual only using two for loops 
        for(int r = 0; r < n_indv; r++)
        {
            if(print){Rcout << "updating trait diffs sq" << "\n";}
            d.traitDiffsSq(r,dead_index) = pow(d.indTraitL(r) - d.indTraitL(dead_index),2);
        }
        for(int c = 0; c < n_indv; c++)
        {
            d.traitDiffsSq(dead_index,c) = pow(d.indTraitL(dead_index) - d.indTraitL(c),2);
        }
        
		// if non neutral, recalculate env_filter_probs
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
            float dp = dispersal_prob(i);

            NumericVector probs = dp * (d.spAbundM / sum(d.spAbundM)) + 
                (1 - dp) * (d.spAbundL / sum(d.spAbundL));

            probs = (abs(probs)+probs)/2;
            IntegerVector s = sample(d.spAbundM.length(), 1, false, probs);
			
            // make i from 0 to phylo.n - 1 (previously 1 to phylo.n)
            int chosen_species = s(0) - 1;
            
            //calculate deviation of trait from previous species trait
            float trait_dev = R::rnorm(0, p.trait_sigma(i));
            float parent_trait_val = 0;
            if(dispersed_this_iter){
                parent_trait_val = d.spTraitM(chosen_species);
            }
            else{
                parent_trait_val = d.spTraitM(chosen_species);
            }
            d.indTraitL(dead_index) = parent_trait_val + trait_dev;
            if(print){Rcout << "set new trait val";}
            
            // the indv that was birthed or immigrated becomes the new species
            d.indSpeciesL(dead_index) = d.nTipsP(0); 
			
            bool print_matrices = false;
            
            // set easy named variables for phylogeny speciation
            NumericMatrix e = d.edgesP;
            
            NumericVector l = d.lengthsP;
            int n = d.nTipsP(0);
            int i = chosen_species;
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
