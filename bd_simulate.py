
import re
import operator
from Bio.Phylo.Applications import RaxmlCommandline
import subprocess
from importlib import reload
from itertools import combinations
from scipy.stats import multivariate_normal
import dendropy
from collections import Counter
from pathlib import Path
from scipy import stats
from functools import *
import pandas as pd
from joblib import Parallel, delayed
import numpy as np
from dendropy.simulate import treesim
from multiprocessing.connection import wait


nsims = 250
lam = 1
mus = np.arange(.1, 1, .1)
heights = np.arange(1, 8, .5)
njobs = 12


def sim_tree_covs(mu: float, tree_height: float):
    """returns list of records with child node height (age) and branch length.

    Args:
        m (mu): death rate
        h (): max tree height

    Returns:
        list: _description_
    """
    c = np.empty(nsims)
    sizes = np.empty(nsims)
    for j in range(nsims):
        t = treesim.birth_death_tree(
            birth_rate=lam, death_rate=mu, max_time=tree_height)
        # output of this func is inconsistent with ordering of internal_nodes iterator
        taxa = len(t)
        t.calc_node_ages()
        root_age = next(t.ageorder_node_iter(descending=True)).age
        x = np.empty([2, taxa-1])
        for i, n in enumerate(t.internal_nodes()):
            bl = n.edge.length
            h = root_age-n.age
            x[:, i] = (bl, h)
        c[j] = stats.spearmanr(x, axis=1).correlation
        sizes[j] = taxa
    return c, sizes


def simulate_traits(n_genes: list, sigma: list, data_dir: Path, nsamps=100):
    for g, sigma in zip(n_genes, sigmas):
        fn = data_dir/f'samples_ngenes{g}.npz'
        if not fn.exists():
            samples = Parallel(12)(
                delayed(multivariate_normal.rvs)
                (np.zeros(len(s)), s, nsamps, i) for i, s in enumerate(sigma)
            )
            np.savez(fn, *samples)
    return samples


def simulate_c_star_matrices(parent_trees: list, data_dir: Path):
    sigmas = []
    try:
        for g in n_genes:
            s = np.load(data_dir/f'sigmas_ngenes{g}.npz')
            sigmas.append([s[f] for f in s.files])
        print('loaded sigmas')
    except FileNotFoundError:
        print('simulating')
        for parent_tree in parent_trees:
            containing_taxa = parent_tree.taxon_namespace
            gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
                containing_taxon_namespace=parent_tree.taxon_namespace,
                num_contained=1)

            gt_params = dict(containing_tree=parent_tree,
                             gene_to_containing_taxon_map=gene_to_species_map)
            gene_trees = Parallel(12)(delayed(treesim.contained_coalescent_tree)(
                **gt_params) for i in range(max_genes))

            covs = Parallel(12)(delayed(phylogenetic_covariance)(gt)
                                for gt in gene_trees)
            s = [np.sum(covs[:g], 0) for g in n_genes]
            sigmas.append(s)
        sigmas = list(zip(*sigmas))
        for g, s in zip(n_genes, sigmas):
            np.savez(data_dir/f'sigmas_ngenes{g}.npz', *s)
    return sigmas


def phylogenetic_covariance(gene_tree):
    """calculate phylogenetic covariance matrix from rooted tree"""
    gene_tree.calc_node_root_distances()
    ns = gene_tree.taxon_namespace
    ntaxa = len(ns)
    dmat = np.zeros((ntaxa, ntaxa))
    for t1 in range(ntaxa):
        for t2 in range(t1+1):
            cov = gene_tree.mrca(taxa=[ns[t1], ns[t2]]).root_distance
            dmat[t1, t2] = dmat[t2, t1] = cov
    return dmat


def sim_bd_trees(
    mu=.5,
    tree_height=6,
    n_parent_trees=100,
):
    """simulate n trees with constant lam.
    NOTE: joblib's loky backend is ok with rng."""

    ptree_params = dict(birth_rate=lam, death_rate=mu, max_time=tree_height)
    parent_trees = Parallel(12)(delayed(treesim.birth_death_tree)(
        **ptree_params) for _ in range(n_parent_trees))
    return parent_trees


parent_dir = Path('/N/project/phyloML/rate_timescaling/data/')

n_parent_trees = 100
n_genes = [50, 100, 250, 500, 1000]
max_genes = max(n_genes)


for mu in np.arange(.3, .9, .1):
    for tree_height in [7, 8, 9]:
        data_dir = parent_dir/f'm{mu:.1f}_h{tree_height}'
        print(data_dir)
        data_dir.mkdir(exist_ok=True)

        try:
            with open(data_dir/'parent_trees.phy', 'r') as f:
                parent_trees_full = [
                    dendropy.Tree.get_from_string(t, 'newick') for t in f]
            if len(parent_trees_full) < n_parent_trees:
                raise IOError('not enough trees, regenerating...')

        except Exception as e:
            print(e, 'simulating...')
            parent_trees_full = sim_bd_trees(mu, tree_height, n_parent_trees)
            with open(data_dir/'parent_trees.phy', 'w') as f:
                f.writelines([t.as_string('newick')
                             for t in parent_trees_full])

        ntaxa = pd.Series(list(map(len, parent_trees_full)))
        parent_trees = [t for t in parent_trees_full if len(t) > 3]
        with open(data_dir/'parent_trees_4plus.phy', 'w') as f:
            f.writelines([t.as_string('newick') for t in parent_trees])

        sigmas = simulate_c_star_matrices(parent_trees, data_dir)
        samples = simulate_traits(n_genes, sigmas, data_dir)

        # tall_trees = [tree.clone(1) for tree in parent_trees]
        # stretch = 5

        # tall_height = tree_height*stretch
        # tall_dir = data_dir.parent/(data_dir.name+f'_s{stretch}')
        # tall_dir.mkdir(exist_ok=True)

        # for tree in tall_trees:
        #     tree.scale_edges(stretch)
        # sigma_tall = simulate_c_star_matrices(tall_trees, data_dir=tall_dir)

        # tall_dir = data_dir.parent/(data_dir.name+f'_s{stretch}')

        # tall_samples = simulate_traits(n_genes, sigma_tall, tall_dir)
        # print(f'finished mu {mu} height {tree_height}')
