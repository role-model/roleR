import msprime
import newick
import numpy as np
import pandas as pd

from collections import Counter, OrderedDict, defaultdict

def py_msprime_simulate(J_m,
                        J,
                        curtime,
                        metaTree,
                        metaAbund,
                        localAbund,
                        localTDiv,
                        alpha,
                        sequence_length,
                        mu,
                        seed=None,
                        return_debug=False,
                        verbose=False
                       ):

    if verbose:
        print(INIT_MSG.format(J_m=J_m, J=J, curtime=curtime, metaTree=metaTree,
                            metaAbund=metaAbund, localTDiv=localTDiv,
                            alpha=alpha, sequence_length=sequence_length, mu=mu))

    ## Sanity check
    if alpha <= 0:
        print("alpha must be >= 1, setting alpha == 1")
        alpha = 1

    ## Enumerate +1 for 1-indexed species names
    meta_sad = {x+1:y*J_m for x,y in enumerate(metaAbund)}

    ## Local community
    local_sad = {int(x):y for x,y in Counter(localAbund).items()}

    ## -1 so lidx indices from 1-idexed species names properly index indo localTDiv array
    lidx = np.array(list(local_sad.keys()), dtype=int)-1
    ## Get tdiv in generations before present
    ## Subtract localTDiv from curtime (current time in iterations) and divide by iterations per generation (J/2)
    localTDiv = np.array(localTDiv)[lidx]
    tdiv = {x+1:(curtime-y)/(J/2) for x,y in zip(lidx, localTDiv)}
    
    ## Make dataframe from list of dictionaries where all keys are species IDs
    full_df = pd.DataFrame([local_sad, meta_sad, tdiv], index=["local_abund", "meta_abund", "tdiv"])
    ## Sort columns ascending
    full_df = full_df[sorted(full_df.columns)]
    ## Update column names to match species tree tip labels
    full_df.columns = ["t{}".format(x) for x in full_df.columns]
    ## Create a new dict for just the species present in the local community
    local_df = full_df.dropna(subset=["local_abund"], axis=1)

    if verbose: print(local_df)
    ## format the metacommunity abundances as a dictionary to pass in to msprime
    meta_Nes = (full_df.loc["meta_abund"]*alpha).to_dict()
    ## Create a default dict so internal nodes have a default Ne. Arbitrarily set to 10,000.
    ## TODO: What is a reasonable default Ne. This is necessary.
    meta_Nes = defaultdict(lambda: 10000, meta_Nes)

    ## Create msprime demography from newick tree
    ## TODO: generation_time is fixed to 1 here. This should be more flexible.    
    demography = msprime.Demography.from_species_tree(_force_ultrametric(metaTree),
                                                  initial_size=meta_Nes,
                                                  time_units="myr",
                                                  generation_time=1)
    
    ## Update msprime demography with local populations
    ## The newick tree only knows about the metacommunity, so we need to update
    ## it to add the branches for the local community pops.
    for sp in local_df:
        meta_sp = f"{sp}_m"
        local_sp = f"{sp}_l"
        local_Ne = local_df[sp]["local_abund"]*alpha
        tdiv = local_df[sp]["tdiv"]
        demography.add_population(name=meta_sp, initial_size=meta_Nes[sp])
        demography.add_population(name=local_sp, initial_size=local_Ne)
        demography.add_population_split(time=tdiv+0.1, derived=[meta_sp, local_sp], ancestral=sp)
        ## Strong colonization bottleneck after colonization
        ##demography.add_instantaneous_bottleneck(time=tdiv, strength=2*local_Ne, population=local_sp)
        demography.add_simple_bottleneck(time=tdiv, population=local_sp, proportion=1)


    ## Make sure demographic events are in proper order ascending order (backward in time)
    demography.sort_events()

    ## Simulate the genealogy
    ## TODO: local sample size is hacked and fixed at 5
    ts = msprime.sim_ancestry(
        samples={f"{sp}_l":5 for sp in local_df},
                demography=demography,
                sequence_length=sequence_length,
                ploidy=1,
                random_seed=seed)
    ## Simulate mutations on the genealogy
    ts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)

    ## Simulation is done, now pull out the information we want. Most of the methods
    ## for extracting data from a treesequence require node IDs, so here we
    ## get a dictionary mapping population ids to the list of node IDs per population
    nodeIDs = defaultdict(lambda: [])
    for n in ts.nodes():
        if n.is_sample():
            nodeIDs[ts.population(n.population).metadata["name"]].append(n.id)

    ## Return all the simulated genotypes and some sumstats as a dataframe
    res = {}
    for popname, idxs in nodeIDs.items():
        ## Split off the _l so the pop name agrees with the names in roleModel localComm
        pname = popname.split("_")[0]
        res[pname] = []
        res[pname].append(ts.diversity(sample_sets=idxs))
        res[pname].append(ts.Tajimas_D(sample_sets=idxs))
        res[pname].append(list(ts.haplotypes(samples=idxs)))
    res_df = pd.DataFrame(res, index=["pi", "TajD", "gtypes"])

    if return_debug:
        return res_df, demography, ts
    else:
        return res_df


def _force_ultrametric(tree):
    """
    Force all leaf nodes to time 0
    Input a newick tree in string format
    Output is a modified newick tree with all leaves having equal total length
    """
    # Parse the newick tree string.
    parsed = newick.loads(tree)
    if len(parsed) == 0:
        raise ValueError(f"Not a valid newick tree: '{tree}'")
    root = parsed[0]

    # Set node depths (distances from root).
    stack = [(root, 0)]
    max_depth = 0
    while len(stack) > 0:
        node, depth = stack.pop()
        if depth > max_depth:
            max_depth = depth
        node.depth = depth
        for child in node.descendants:
            stack.append((child, depth + child.length))
    for node in root.walk():
        if node.is_leaf:
            ## Add offset to node.length to foce all nodes to fall at time 0
            ## The offset is the difference between the depth of this node
            ## and the max_depth of the deepest leaf node.
            node.length = node.length + (max_depth - node.depth)
    return root.newick


INIT_MSG = """
    J_m  - {J_m}
    J - {J}
    curtime - {curtime}
    metaTree - {metaTree}
    metaAbund - {metaAbund}
    localTDiv - {localTDiv}
    alpha - {alpha}
    sequence_length - {sequence_length}
    mu - {mu}
"""

