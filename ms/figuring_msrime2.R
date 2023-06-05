#remotes::install_github("rstudio/reticulate")
library(reticulate)

conda_list()
use_condaenv("r-reticulate")
py_list_packages()
msprime <- reticulate::import('msprime')
newick <- reticulate::import('python-newick')
collections <- reticulate::import('collections')

metaSAD <- list(A = 1000, B = 1000, C = 1000)
lambda <- function() mean(unlist(metaSAD))
initSize <- collections$defaultdict(lambda)
initSize$update(reticulate::dict(metaSAD))

x <- ape::read.tree(text='((A:2, B:2):1, C:3);')
plot(x)
ape::read
d <- msprime$Demography$from_species_tree('((A:2, B:2):1, C:3)', 
                                          time_units = 'myr', 
                                          initial_size = initSize, 
                                          generation_time = 1, growth_rate = 0)
# after this, contains all SPECIES including those originally in meta that moved to local, 
# and those that originated in local 
# and those that went extinct in local 

# what needs to be added are branches for populations OF THE SAME SPECIES that split off meta nad went to local
# ie dispersal

# add 
d$add_population(name = 'A_meta', initial_size = 10000)
d$add_population(name = 'A_local', initial_size = 1000)
d$add_population_split(time = 1, derived = c('A_meta', 'A_local'), 
                       ancestral = 'A')

d$sort_events()

ts <- msprime$sim_ancestry(reticulate::dict(list(A_meta = 3, A_local = 4, B = 2, C = 1)), demography = d, ploidy = 1)

cat(ts$draw_text())


d$add_population_split()


d$add_population_parameters_change

d$add_population(name = 'A_loc', initial_size = 1000)
d$add_population_split



# remember add_census

