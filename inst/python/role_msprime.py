import msprime
import tskit
import newick
import numpy as np
import pandas as pd
from collections import defaultdict


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



def sim_role(Jm, 
             curtime, 
             nwk, 
             meta_abund, 
             local_sad, 
             local_hm, 
             J_harm,
             gens, 
             tdiv, 
             alpha, 
             m, 
             sequence_length, 
             mu, 
             seed = None):
    """
    Simulate population genetics based on the ROLE model framework.
    
    Args:
        Jm: Meta community size (int)
        curtime: Current time in iterations (double)
        nwk: Newick tree string 
        meta_abund: Meta community abundances (list)
        local_sad: Local community abundances (numpy array)
        local_hm: Local community harmonic mean abund (numpy array)
        J_harm: Harmonic mean of total local community size (double)
        tdiv: Colonization times for local populations (numpy array)
        alpha: Scaling factor for converting abundance to Ne (>= 1)
        m: migration rate
        sequence_length: Length of the genomic region to simulate
        mu: Mutation rate per base
        seed: Random seed for reproducibility
        
    Returns:
        Pandas DataFrame containing diversity statistics and genotypes
    """
    
    # convert matrix-like objects with species abundances info to list of dicts
    local_sad_ts = local_sad.to_dict("records")
    local_hm_ts = local_hm.to_dict("records")
    
    # enumerate +1 because species are 1-indexed from R
    meta_Nes = {x+1:y*alpha for x,y in enumerate(meta_abund)}
    
    # create a default dict so internal nodes have a default Ne
    dNe = round(Jm / len(meta_Nes) * alpha)
    meta_Nes = defaultdict(lambda: dNe, meta_Nes)

    # proportional metacomm abundances (needed for calc imm rate)
    meta_p = {f"{x+1}":(y / Jm) for x,y in enumerate(meta_abund)}

    # NOTE: all time objects **must** be sorted from ancient to current
    tdiv = pd.DataFrame(tdiv, columns = local_sad.columns)
    local_tdiv_ts = tdiv.to_dict("records")

    # Ensure tree is ultrametric
    phylo_meta = _force_ultrametric(nwk)
    
    # create msprime demography from newick tree
    demography = msprime.Demography.from_species_tree(
        phylo_meta,
        initial_size = meta_Nes,
        time_units = "myr",
        generation_time = 1
    )
    
    # update msprime demography with local populations budding off from their
    # metacommunity parental populations and create a list of sample sets that
    # determine when observations are made

    sampset = [] # empty list to hold sample sets
    
    for sp in meta_p.keys():
        # we track the ancestral population with `this_anc` which will update
        # after each addition of a local population
        # we hard code `t` because we pass a phylo with all "t" prefixes
        this_anc = f"t{sp}"

        # object to track current local and meta pops (only initializing here)
        this_loc = f"l_0_{sp}"
        this_meta = f"m_0_{sp}"

        for t in range(len(local_sad_ts)):
            # `t` is the *index* for snapshots

            # extract abundance for this sp in this t
            this_abund = local_sad_ts[t][sp]

            # only bother with adding to demography if non-0 abundance
            if this_abund > 0:
                # tags to note _m_eta/_l_ocal and _t_ime
                meta_sp = f"m_{t}_{sp}"
                local_sp = f"l_{t}_{sp}"
                
                # since abund is > 0 we will be updating migration rate, so
                # calcualte it here for future use
                m_rate = m * meta_p[sp] * J_harm / (2 * local_hm_ts[t][sp])


                # set value of `add_yn` (determine if adding populations) based
                # on whether t == 0 (always add in initial time step)
                add_yn = True if t == 0 else False

                # if t!= 0, determine if we still add a new pop for this lineage
                if (not add_yn) and local_tdiv_ts[t][sp] != local_tdiv_ts[t - 1][sp]:
                    add_yn = True

                if add_yn:
                    # track the pops with `this_loc` and `this_meta` which will
                    # update at each addition of a local comm population; we
                    # have to track separately because we might need to take a
                    # sample of this same pop across multiple snapshots
                    this_loc = local_sp
                    this_meta = meta_sp

                    # needed info for adding pops
                    local_Ne = local_hm_ts[t][sp] * alpha # harmonic mean = Ne
                    
                    # gens accumulates C++ floating-point error and can
                    # slightly exceed curtime, making the difference negative
                    split_time = max(1e-09, curtime - local_tdiv_ts[t][sp])

                    demography.add_population(
                        name = this_meta,
                        initial_size = meta_Nes[sp]
                    )

                    demography.add_population(
                        name = this_loc,
                        initial_size = local_Ne
                    )

                    # when adding a pop need to make sure the ancestral pop is
                    # the previous metacomm pop, this may not be the same as
                    # `sp` so we use `this_anc` which initializes as `sp` but
                    # then updates
                    demography.add_population_split(
                        time = split_time,
                        derived = [this_meta, this_loc],
                        ancestral = this_anc # note we use `this_anc`
                    )

                    # strong colonization bottleneck cause only one individual
                    # at time of colonization
                    demography.add_simple_bottleneck(
                        time = split_time,
                        population = this_loc,
                        proportion = 1
                    )

                    # add migration
                    # source is local, dest is meta (cause backward time)
                    demography.set_migration_rate(
                        source = this_loc,
                        dest = this_meta,
                        rate = m_rate
                    )

                    # only now update `this_anc`
                    this_anc = this_meta

                else:
                    # if weʻre not adding a lineage, but abund is non-zero, we
                    # still need to update population size and migration rate

                    demography.add_population_parameters_change(
                        time = curtime - gens[t],
                        population = this_loc,
                        initial_size = local_hm_ts[t][sp] * alpha
                    )

                    demography.add_migration_rate_change(
                        time = curtime - gens[t],
                        source = this_loc,
                        dest = this_meta,
                        rate = m_rate
                    )

                # create sample set for this sp in this time step whether or not
                # we are adding a new lineage
                # if we add a new lineage, `this_loc` will reflect current timestep
                # if we donʻt add new, `this_loc` will reflect last timestep
                # where we did add a lineage

                samp_time = curtime - gens[t] # because backward time

                # samp_time must be strictly less than the population's split
                # time; when a species first appears gens[t] == origin gen so
                # samp_time == split_time -- clip by a tiny epsilon
                split_time_this_loc = max(1e-09, curtime - local_tdiv_ts[t][sp])
                if samp_time >= split_time_this_loc:
                    samp_time = split_time_this_loc - 1e-09

                this_samp = msprime.SampleSet(
                    this_abund,
                    population = this_loc,
                    time = samp_time
                )

                sampset.append(this_samp)


    # make sure demographic events are in proper ascending order
    # (backward in time)
    demography.sort_events()

    # simulate ancestry with specific sampling times
    ts = msprime.sim_ancestry(
        samples = sampset,
        demography = demography,
        ploidy = 1,
        sequence_length = sequence_length, 
        recombination_rate = 0
    )
    
    # simulate mutations
    mts = msprime.sim_mutations(ts, rate = mu)

    # node table with all info
    ntab = mts.tables.nodes

    # dataframe of needed info for sample nodes only
    df_ind = pd.DataFrame({
        "ids": ntab.individual,
        "pops": ntab.population,
        "time": ntab.time})
    df_ind = df_ind[df_ind["ids"] > -1] # NOTE!! this only works if `ploidy = 1`

    # flip time from backward direction to forward
    df_ind["time"] = curtime - df_ind["time"]

    # snap forward times to nearest gens value; the epsilon clipping above
    # produces times like gens[t] + 1e-9 which must map back to gens[t]
    gens_arr = np.array(gens)
    df_ind["time"] = df_ind["time"].apply(
        lambda t: float(gens_arr[np.argmin(np.abs(gens_arr - t))])
    )

    # get actual names of the populations
    pop_names = pd.Series([pop.metadata["name"] for pop in mts.tables.populations])
    df_ind["name"] = df_ind["pops"].map(pop_names)

    # get just species names from "name" column
    df_ind["sp_id"] = df_ind["name"].str.split('_').str[-1]

    # get simulated sequences (this hangs a suprising while)
    df_ind["seq_alignment"] = list(mts.alignments(
        samples = df_ind["ids"].tolist(),
        reference_sequence = tskit.random_nucleotides(ts.sequence_length)
    ))

    # aggregate df so rows are species and column of ids is a column of lists
    df_spp = df_ind.groupby(["time", "sp_id"])['ids'].apply(list).reset_index()

    # get diversities and Tajima for each pop in each time
    df_spp["pi"] = mts.diversity(df_spp["ids"].tolist())
    df_spp["tajimas_D"] = mts.Tajimas_D(df_spp["ids"].tolist())

    # all we need from `df_ind` is time, sp_id, and seqs
    df_ind = df_ind[["time", "sp_id", "seq_alignment"]]

    # all we need from `df_ind` is time, sp_id, pi, and taj
    df_spp = df_spp[["time", "sp_id", "pi", "tajimas_D"]]


    return df_ind, df_spp


def test_force_ultrametric(tree):
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
