import msprime
import numpy as np
import pandas as pd


def sim_single(
    N_island=1_000,
    T_found=1_000,
    n_samples=100,
    seq_length=500,
    mu=1e-5,
    seed=None,
):
    """
    Simulate a single-island model with msprime.

    Args:
        N_island: Haploid effective size of the island population.
        n_samples: Number of haploid genomes sampled from the island.
        seq_length: Sequence length in base pairs.
        mu: Mutation rate per bp per generation.
        seed: Random seed for reproducibility.

    Returns:
        DataFrame with 1 row and with columns:
            N_mainland      - number of sampled haploid genomes
            N_island     - sequence length (bp)
            m_main_to_island        - number of segregating sites
            T_found             - nucleotide diversity per site
            n_samples            - the msprime TreeSequence with mutations
            seq_length
            mu
            n_sites
            pi
    """
    # --- Demography -----------------------------------------------------------
    demography = msprime.Demography()
    demography.add_population(name="island", initial_size=N_island)
    
    demography.add_simple_bottleneck(time = T_found,
                                     population = "island",
                                     proportion = 1)
    
    # --- Ancestry -------------------------------------------------------------
    ts = msprime.sim_ancestry(
        samples={"island": n_samples},
        ploidy=1,
        demography=demography,
        sequence_length=seq_length,
        recombination_rate=0,
        random_seed=seed,
    )

    # --- Mutations ------------------------------------------------------------
    mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)

    out = pd.DataFrame([{
        "N_island": N_island,
        "n_samples": n_samples,
        "seq_length": int(mts.sequence_length),
        "mu": mu,
        "n_sites": mts.num_sites,
        "pi": float(mts.diversity())
    }])
    
    return out








def sim_island(
    N_mainland=100_000,
    N_island=1_000,
    m_main_to_island=1e-4,
    T_found=1_000,
    n_samples=100,
    seq_length=500,
    mu=1e-5,
    alpha=100,
    seed=None,
):
    """
    Simulate a simple mainland-island model with msprime.

    Args:
        N_mainland: Haploid effective size of the mainland source population.
        N_island: Haploid effective size of the island population.
        m_main_to_island: Forwards-time immigration rate mainland -> island.
        T_found: Generations ago the island was founded from the mainland.
        n_samples: Number of haploid genomes sampled from the island.
        seq_length: Sequence length in base pairs.
        mu: Mutation rate per bp per generation.
        seed: Random seed for reproducibility.

    Returns:
        DataFrame with 1 row and with columns:
            N_mainland      - number of sampled haploid genomes
            N_island     - sequence length (bp)
            m_main_to_island        - number of segregating sites
            T_found             - nucleotide diversity per site
            n_samples            - the msprime TreeSequence with mutations
            seq_length
            mu
            n_sites
            pi
    """
    # --- Demography -----------------------------------------------------------
    demography = msprime.Demography()
    demography.add_population(name="mainland", initial_size=N_mainland, 
                              initially_active = True)
    
    
    demography.add_population(name="island", initial_size=N_island)
    
    demography.add_simple_bottleneck(time = T_found,
                                     population = "island",
                                     proportion = 1)
                                     
    # Island founded from the mainland T_found generations ago
    demography.add_population_split(time=T_found, derived=["island"],
                                    ancestral="mainland")

    
                    
    # msprime migration rates are defined BACKWARDS in time: the rate from
    # `source` to `dest` is the per-generation rate at which an ancestral lineage
    # sitting in `source` jumps to `dest` as we look back. Forwards in time that
    # is movement of individuals from `dest` into `source`. We want forwards-time
    # immigration mainland -> island, so backwards in time island lineages must
    # move to the mainland:
    
    demography.set_migration_rate(source="island", dest="mainland",
                                  rate=m_main_to_island)
    


    # --- Ancestry -------------------------------------------------------------
    ts = msprime.sim_ancestry(
        samples={"island": n_samples},
        ploidy=1,
        demography=demography,
        sequence_length=seq_length,
        recombination_rate=0,
        random_seed=seed,
    )

    # --- Mutations ------------------------------------------------------------
    mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed)

    out = pd.DataFrame([{
        "N_mainland": N_mainland,
        "N_island": N_island,
        "m_main_to_island": m_main_to_island,
        "effective_m": m_main_to_island,
        "T_found": T_found,
        "n_samples": n_samples,
        "seq_length": int(mts.sequence_length),
        "mu": mu,
        "n_sites": mts.num_sites,
        "pi": float(mts.diversity()),
        "tajimas_D": float(mts.Tajimas_D())
    }])
    
    return out

