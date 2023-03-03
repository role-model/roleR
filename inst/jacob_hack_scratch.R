# number of individuals and species in the local and meta
individuals_local = 1 # 'function',
individuals_meta = 1 # 'integer',
species_meta = 1 # 'integer',

# probs of speciation, extinction, dispersal
speciation_local = 1# 'function',
speciation_meta = 1# 'double',
extinction_meta = 1# 'double',
dispersal_prob = 1# 'function',

trait_sigma = 1# 'double',
env_sigma = 1 # 'double',
comp_sigma = 1# 'double',
neut_delta = 1# 'double',
env_comp_delta = 1# 'double',

# gene evolution simulations
mutation_rate = 1# 'double',
equilib_escape = 1# 'double',
alpha = 1# 'function',
num_basepairs = 1# 'integer',

init_type = 'bridge_island' #'character', 

# iterations 
niter = 1# 'integer', 
niterTimestep = 1# 'integer'

# NEW
individuals_local = 100
individuals_meta = 1000
species_meta = 100
speciation_local = 0.1
speciation_meta = 0.1
extinction_meta = 0.05
dispersal_prob = 0.1
trait_sigma=1
env_sigma=1
comp_sigma = 0.5
neut_delta=1
env_comp_delta=1
mutation_rate=0
equilib_escape = 1
alpha=50
num_basepairs = 250
init_type = 'oceanic_island'
niter = 100
niterTimestep = 10

all_params <- list(individuals_local,
                   individuals_meta,
                   species_meta,
                   
                   speciation_local,
                   speciation_meta,
                   extinction_meta,
                   dispersal_prob,
                   
                   trait_sigma,
                   env_sigma,
                   comp_sigma,
                   neut_delta,
                   env_comp_delta,
                   
                   mutation_rate,
                   equilib_escape,
                   alpha,
                   num_basepairs,
                   
                   init_type, 
                   niter, 
                   niterTimestep)

p <- roleParams2(individuals_local = 100, individuals_meta = 1000, species_meta = 10, 
                speciation_local = 0.1, speciation_meta = 0.1, extinction_meta = 0.05, dispersal_prob = 0.1,
                
                trait_sigma=1, env_sigma=1, comp_sigma = 0.5, neut_delta=1, env_comp_delta=1,
                mutation_rate=0,equilib_escape = 1, alpha=0, num_basepairs = 250,
                init_type = 'oceanic_island', niter = 100, niterTimestep = 10)
