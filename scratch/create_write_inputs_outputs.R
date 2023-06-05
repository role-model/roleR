# create some example inputs and outputs
role <- dummyExperiment()
#role <- dummyExperiment(run=TRUE)

# writes test_params.roleparams AND metadata file test_params.txt 
writeRoleExperiment(role, dir="data/example_out_in", fileName="test_exp")

# writes test_params.roleparams AND metadata file test_params.txt 
writeRoleParams(role@params,dir="data/example_out_in",fileName="test_params")

# read in params and sim
experiment_in = readRDS("data/example_out_in/test_exp.roleexperiment")
params_in = readRDS("data/example_out_in/test_params.roleparams")

# access the first modelRun from the experiment
experiment_in@modelRuns[[1]]

# access the params from the experiment (the params used in this experiment)
experiment_in@params

# access params from the params file
params_in

# this does not work currently (working thru a bug) but this is how you would use a params object to 
#   start or replicate an experiment, creating a new roleExperiment object by running a simulation using the params 
test_role <- roleExperiment(params=params_in)