from multiprocessing.connection import wait
from dendropy.simulate import treesim
# import msprime
import numpy as np
from joblib import Parallel, delayed
import pandas as pd
import operator
from functools import *

t = treesim.birth_death_tree(birth_rate=1.0, death_rate=0.5, ntax=100)
nsims = 200
lam = 1
mus = np.arange(.05, 1, .05)
heights = np.arange(1, 5, .5)
njobs = 12


def sim_tree_covs(m: float, h: float):
    """returns list of records with child node height (age) and branch length.

    Args:
        m (mu): death rate
        h (): max tree height

    Returns:
        list: _description_
    """
    c = np.empty(nsims)
    for j in range(nsims):
        t = treesim.birth_death_tree(birth_rate=lam, death_rate=m, max_time=h)
        # output of this func is inconsistent with ordering of internal_nodes iterator
        taxa = len(t)
        t.calc_node_ages()
        root_age = t.nodes()[0].age
        x = np.empty([2, taxa-1])
        for i, n in enumerate(t.internal_nodes()):
            bl = n.edge.length
            h = root_age-n.age
            x[:, i] = (bl, h)
        c[j] = np.corrcoef(x)[1, 0]
        if np.isnan(c[j]):
            print(np.corrcoef(x))
            print(taxa, t.as_string('newick'))
    return c


def sim_trees(m: float, h: float):
    """returns list of records with child node height (age) and branch length.

    Args:
        m (mu): death rate
        h (): max tree height

    Returns:
        list: _description_
    """
    ages = []
    for _ in range(nsims):
        t = treesim.birth_death_tree(birth_rate=lam, death_rate=m, max_time=h)
        # output of this func is inconsistent with ordering of internal_nodes iterator
        taxa = len(t)
        t.calc_node_ages()
        root_age = t.nodes()[0].age
        for n in t.internal_nodes():
            record = {
                'mu': m,
                'height': h,
                'branch_length': n.edge.length,
                'age': n.age,
                'depth': root_age-n.age,
                'taxa': taxa}
            ages.append(record)
    return ages


sims = Parallel(njobs)(delayed(sim_trees)(m, h) for m in mus for h in heights)
sims = pd.DataFrame.from_records(reduce(operator.add, sims))
sims.to_pickle('bd_trees_200.pd.gz')

c = (sims
     .query('taxa>2')
     .groupby(['mu', 'height'])[['depth', 'branch_length']
                                ]
     .corr('spearman')
     .depth
     .xs('branch_length', level=-1)
     #  .unstack()
     )


initial_size = 1e6
t.scale_edges(initial_size)
demography = msprime.Demography.from_species_tree(
    t.as_string('newick'), initial_size)

# mutation model
pi = np.array([0.1, 0.2, 0.3, 0.4])
hky = msprime.HKY(kappa=0.75, equilibrium_frequencies=pi)
P = hky.transition_matrix
u = 2.5e-7
mu = u / (1 - np.sum(pi * np.diag(P)))

print(f"Mutation rate should be increased by a factor of {mu/u:.2f}.")

ts = msprime.sim_ancestry(1000,
                          demography=demography,
                          population_size=1000,
                          sequence_length=1e7,
                          recombination_rate=1e-8,
                          random_seed=5)
mts = msprime.sim_mutations(ts, rate=mu, model=hky, random_seed=27)
theta = mts.diversity()
print(f"Genetic diversity: {theta}.")


def cov(t1: float, t2: float, sigma2: float):
    """ t1: mrca of species A,B
        t2: root of tree
        sigma2: """
    ibl = t2-t1
    e1 = np.exp(t1)
    e2 = np.exp(t2)
    waiting_time = np.exp(-ibl)
    return sigma2*((1-waiting_time)*(e2*ibl/(e2-e1)) + waiting_time/3)
