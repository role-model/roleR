# example timeseries for connor
library(roleR)
# read in serialized experiment
test_exp <- readRDS("data/example_out_in/test.roleexperiment")

# add hill summary stats at entropies 3,4,5
test_exp <- addHillStats(test_exp,entrScales=c(3,4,5))

# get a timeseries of all summary stats for the experiment
ts <- getTimeseries(test_exp,runNum=1,type="summary_stat")
# a vector of the hill number based on abundance with entropy 3 
# NOTE - timeseries copying is still not working, so these values are copies of the first time step and are thus identical, but will not be in the future
ts$abundance_3

# alternatively, can just do 
ts <- getTimeseries(test_exp,runNum=1,type="summary_stat",valueName="abundance_3")
# to get vector directly in one line 