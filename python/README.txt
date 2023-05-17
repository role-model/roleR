Python files and jupyter notebooks for implementing msprime simulations
within the RoLE model. There are 4 files in this directory which are explained below.

role_msprime.py - A python file containing the function `py_msprime_simulate`. This
takes in a bunch of different information necessary to generate the simulations. It
is expected that the bookkeeping for driving the call to this function will mostly
be handled by R (see the `r-role-msprime.ipynb` for the necessary R function).

r-role-msprime.ipynb - An R notebook demonstrating the use of the `sim_seqs` function
to run the python code. `sim_seqs` takes one argument which is a roleModel object which
has been run. For each of the `modelSteps` it strips out all the necessary information
and calls `py_msprime_simulate` and then populates the roleModel with the resulting
pi and haplotype data.

role-msprime-debug.ipynb - A simple-as-possible notebook for debugging the msprime
simulations. One cell invokes R and runs a role model and exports all the necessary
data into python. One cell (python) runs `py_msprime_simulate` and shows how to 
access a bunch of debug information, and also how to plot the resulting tree/mutations.

role-msprime-dev.ipynb - A notebook for developing the role_msprime.py code. This is
included for completeness so the development process is transparent, but if you are not
Isaac you probably should not look at this with an intention of trying to undestand it
(it's messy).

