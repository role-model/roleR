# create some example inputs and outputs
sim <- dummySim()
writeSim(sim,fileName="test_sim")

writeParams(sim@params)

# read in params

params_in = readRDS("test_params")
