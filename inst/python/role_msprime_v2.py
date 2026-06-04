import msprime
import newick
import numpy as np
import pandas as pd
from collections import Counter, defaultdict
from typing import List, Dict, Tuple, Optional, Union

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
            ## Add offset to node.length to force all nodes to fall at time 0
            ## The offset is the difference between the depth of this node
            ## and the max_depth of the deepest leaf node.
            node.length = node.length + (max_depth - node.depth)
    return root.newick

def load_tree(newick_path: str) -> str:
    """
    Load a phylogenetic tree from a Newick format file.
    
    Args:
        newick_path: Path to the Newick file
        
    Returns:
        A newick tree string
        
    Raises:
        ValueError: If the tree cannot be loaded
    """
    try:
        with open(newick_path, 'r') as f:
            tree_str = f.read().strip()
        return tree_str
    except Exception as e:
        raise ValueError(f"Error loading tree: {e}")

def sim_role(J_m,
             J,
             curtime,
             metaTree,
             metaAbund,
             localAbund,
             spAbundHarmMean,
             localTDiv,
             alpha,
             sequence_length,
             mu,
             seed=None,
             return_debug=False,
             verbose=False):
    """
    Simulate population genetics based on the ROLE model framework.
    
    Args:
        J_m: Meta community size
        J: Local community size
        curtime: Current time in iterations
        metaTree: Newick tree string or path to Newick file
        metaAbund: Meta community abundances
        localAbund: Local community abundances 
        spAbundHarmMean: Harmonic mean of species abundances
        localTDiv: Colonization times for local populations 
        alpha: Scaling factor for converting abundance to Ne (>= 1)
        sequence_length: Length of the genomic region to simulate
        mu: Mutation rate per base
        seed: Random seed for reproducibility
        return_debug: Whether to return additional debug information
        verbose: Whether to print detailed information
        
    Returns:
        Pandas DataFrame containing diversity statistics and genotypes
    """
    if verbose:
        print(f"""
    J_m  - {J_m}
    J - {J}
    curtime - {curtime}
    metaTree - {metaTree}
    metaAbund - {metaAbund}
    localAbundHmean - {spAbundHarmMean}
    localTDiv - {localTDiv}
    alpha - {alpha}
    sequence_length - {sequence_length}
    mu - {mu}
    """)
    
    # Ensure tree is ultrametric
    metaTree = _force_ultrametric(metaTree)
    
    # Enumerate +1 for 1-indexed species names
    meta_sad = {x+1:y*J_m for x,y in enumerate(metaAbund)}

    # Local community. Raw abundances. In current form only used to get the
    # species ids for the local community (which will be used to grab harmonic means below)
    local_sad = {int(x):y for x,y in Counter(localAbund).items()}
    if verbose: 
        print("local_sad Raw - ", local_sad)

    # -1 so lidx indices from 1-indexed species names properly index into localTDiv array
    # and the spAbundHarmMean array (sorted for the sake of sanity)
    lidx = np.sort(np.array(list(local_sad.keys()), dtype=int)-1)

    # Maintain the dictionary format for the local community mapping species ids to 
    # harmonic mean abundances
    # If harmonic mean species abundance is all zeros then fall back to using the raw abundance
    # (primarily for time zero).
    if np.sum(spAbundHarmMean) > 0:
        hmeans = np.array(spAbundHarmMean)[lidx]
        local_sad = {x+1:y for x,y in zip(lidx, hmeans)}
        if verbose: 
            print("local_sad Hmean - ", local_sad)

    # Get tdiv in generations before present
    # Subtract localTDiv from curtime (current time in iterations) and divide by iterations per generation (J/2)
    localTDiv = np.array(localTDiv)[lidx]
    tdiv = {x+1:(curtime-y)/(J/2) for x,y in zip(lidx, localTDiv)}
    
    # Make dataframe from list of dictionaries where all keys are species IDs
    full_df = pd.DataFrame([local_sad, meta_sad, tdiv], index=["local_abund", "meta_abund", "tdiv"])
    # Sort columns ascending
    full_df = full_df[sorted(full_df.columns)]
    # Update column names to match species tree tip labels
    full_df.columns = ["t{}".format(x) for x in full_df.columns]
    # Create a new dict for just the species present in the local community
    local_df = full_df.dropna(subset=["local_abund"], axis=1)

    if verbose: 
        print(local_df)
        
    # Format the metacommunity abundances as a dictionary to pass in to msprime
    meta_Nes = (full_df.loc["meta_abund"]*alpha).to_dict()
    # Create a default dict so internal nodes have a default Ne. Arbitrarily set to 10,000.
    dNe = J / len(metaAbund) * alpha
    meta_Nes = defaultdict(lambda: dNe, meta_Nes)

    # Create msprime demography from newick tree
    demography = msprime.Demography.from_species_tree(
        metaTree,
        initial_size = meta_Nes,
        time_units = "myr",
        generation_time = 1
    )
    
    # Update msprime demography with local populations
    # The newick tree only knows about the metacommunity, so we need to update
    # it to add the branches for the local community pops.
    for sp in local_df:
        meta_sp = f"{sp}_m"
        local_sp = f"{sp}_l"
        local_Ne = local_df[sp]["local_abund"]*alpha
        split_time = local_df[sp]["tdiv"]
        
        # Add meta and local populations
        demography.add_population(name=meta_sp, initial_size=meta_Nes[sp])
        demography.add_population(name=local_sp, initial_size=local_Ne)
        
        # Add split between meta and local populations
        demography.add_population_split(
            time=split_time+0.1, 
            derived=[meta_sp, local_sp], 
            ancestral=sp
        )
        
        # Strong colonization bottleneck after colonization
        # Using ROLE implementation for bottleneck
        demography.add_simple_bottleneck(
            time=split_time, 
            population=local_sp, 
            proportion=1
        )

    # make sure demographic events are in proper order ascending order 
    # (backward in time)
    demography.sort_events()

    # get sample size from abund
    samples = {}
    for sp in local_df:
        local_sp = f"{sp}_l"
        
        # get the sample size from local_sad
        sample_size = int(local_df[sp]["local_abund"])
        
        # ensure at least 1 sample...need to?
        # sample_size = max(1, sample_size)
        
        samples[local_sp] = sample_size
        
    # simulate the genealogy
    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=sequence_length,
        ploidy=1,
        random_seed=seed
    )
    
    
    
    # Simulate mutations on the genealogy
    ts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)

    # Get node IDs for each population
    nodeIDs = defaultdict(list)
    for n in ts.nodes():
        if n.is_sample():
            nodeIDs[ts.population(n.population).metadata["name"]].append(n.id)

    # Return all the simulated genotypes and some sumstats as a dataframe
    # Following ROLE implementation
    res = {}
    for popname, idxs in nodeIDs.items():
        # Split off the _l so the pop name agrees with the names in roleModel localComm
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

# Example usage (for documentation purposes)
if __name__ == "__main__":
    # This code only runs if the script is executed directly (not imported)
    print("Example usage of sim_role function:")
    print("results = sim_role(")
    print("    J_m=1000,")
    print("    J=400,")
    print("    curtime=5000,")
    print("    metaTree='((t1:1,t2:1):1,t3:2);',")
    print("    metaAbund=[0.5, 0.3, 0.2],")
    print("    localAbund=[1, 1, 2, 2, 3],")
    print("    spAbundHarmMean=[50, 30, 20],")
    print("    localTDiv=[3000, 2000, 4000],")
    print("    alpha=2,")
    print("    sequence_length=10000,")
    print("    mu=1e-8,")
    print("    seed=42")
    print(")")
