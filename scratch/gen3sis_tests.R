
library(gen3sis)

sim <- run_simulation(config = "data/gen3sis/config/SouthAmerica/config_southamerica.R",
                      landscape = "data/gen3sis/landscape/SouthAmerica", output_directory = tempdir(),
                      call_observer = 1, verbose = 0)

# read rds problem

# 1. run gen3sis simulation for one step
# 2. create role metacommunity for each (square? species?) 
# 3. 

# simulate traits mostly? 

# would role just simulate the internal dynamics of each single-species population?
# then change the pops

# or simulate linked populations in tandem? like a kernel or something? 
# pretty easy 

