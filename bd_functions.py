import json
from importlib import reload
from dendropy import Tree
from dendropy.model import continuous, reconcile
import os
from functools import *
import random, string
from importlib import reload
from typing import Union
from Bio import Phylo
from io import StringIO
import operator
import re
from collections import Counter, defaultdict
from functools import *
from itertools import *
from multiprocessing.connection import wait
from pathlib import Path
from sys import argv

import dendropy
from dendropy.calculate import treecompare
import joblib  # import dump,load
import matplotlib.pyplot as plt
import numba as nb

# import msprime
import numpy as np
import pandas as pd
import seaborn as sns
from dendropy.simulate import treesim
from joblib import Parallel, delayed
from scipy import linalg, stats
from scipy.special import factorial, gammaln
from scipy.stats import multivariate_normal
from wrapt_timeout_decorator import timeout


import re

rx = re.compile("\/l([\d.]+)\/m([\d.]+)\/")

nsims = 250
lam = 1
mus = np.arange(0.1, 1, 0.1)
heights = np.arange(1, 8, 0.5)
njobs = 24

n_parent_trees = 100
n_genes = [5, 10, 25, 50]
nsims = 500
njobs = 24

TIMEOUT = 15 * 60


def get_branch_lengths_and_heights(t):
    """returns array of branch lengths and depths (distance from root) of the lower node."""
    t.calc_node_ages()
    taxa = len(t)
    root_age = next(t.ageorder_node_iter(descending=True)).age
    x = np.empty([2, taxa - 1])
    for i, n in enumerate(t.internal_nodes()):
        bl = n.edge.length
        h = root_age - n.age
        x[:, i] = (bl, h)
    return x

def symmetric_log_det_distance(X, Y, normalize=False):
    """https://lear.inrialpes.fr/people/cherian/papers/jbld-pami_final.pdf"""
    r = np.log(linalg.det((X + Y) / 2)) - np.log(linalg.det(X @ Y)) / 2
    if normalize:
        r /= len(X)
    return r


def frobenius_distance(X, Y, normalize=False):
    r = linalg.norm(X - Y, "fro")
    if normalize:
        r /= len(X)
    return r


def spectral_distance(X, Y, normalize=False):
    r = linalg.norm(X - Y, 2)
    if normalize:
        r /= len(X)
    return r


metrics = {
    # 'ld':symmetric_log_det_distance,
    "fro": frobenius_distance,
    "spec": spectral_distance,
}


@lru_cache(maxsize=300)
def rho(k, t):
    return rho_numba(k, t)


@nb.njit
def rho_numba(k, t, Ne=1):
    return np.exp((-(k**2) + k) * t / (2 * Ne))


@lru_cache(maxsize=300)
def upto(a, j):
    return upto_numba(a, j)


@nb.njit(parallel=True)
def upto_numba(a, j):
    if j < 1:
        return 1
    r = 1
    for k in nb.prange(a, a + j):
        r *= k
    return r


@lru_cache(maxsize=300)
def downto(a, j):
    return downto_numba(a, j)


@nb.njit(parallel=True)
def downto_numba(a, j):
    if j < 1:
        return 1
    r = 1
    for k in nb.prange(a - j + 1, a + 1):
        r *= k
    return r


@lru_cache(maxsize=300)
def lupto(a, j):
    return lupto_numba(a, j)


@nb.njit
def lupto_numba(a, j):
    if j < 1:
        return 0
    return np.sum(np.log(np.arange(a, a + j)))


@lru_cache(maxsize=300)
def ldownto(a, j):
    return ldownto_numba(a, j)


@nb.njit
def ldownto_numba(a, j):
    if j < 1:
        return 0
    return np.sum(np.log(np.arange(a, a - j, -1)))


# @njit(parallel=True)


@lru_cache(maxsize=400)
def gterm(i, j, k):
    if not j <= i:
        raise ValueError(f"i={i},j={j}")

    if j > 1:
        sgn = (k - j) % 2 and -1 or 1
    else:
        sgn = k % 2 and -1 or 1

    if i < 15:
        try:
            res = (2 * k - 1) * downto(i, k) / upto(i, k)
            if j > 1:
                res *= upto(j, k - 1) / (factorial(j) * factorial(k - j))
            # print(i,j,k,res)
            return res * sgn
        except:
            print(i, j, k)
            raise
    else:
        try:
            res = np.log(2 * k - 1) + ldownto(i, k) - lupto(i, k)
            if j > 1:
                res += lupto(j, k - 1) - gammaln(j + 1) - gammaln(k - j + 1)
            return np.exp(res) * sgn
        except:
            print(i, j, k)
            raise


# @njit(parallel=True)


@lru_cache(maxsize=400)
def g(i, j, t):
    # tavare p132
    if j == 1:
        try:
            return 1 - sum(gterm(i, j, k) * rho(k, t) for k in range(2, i + 1))
        except:
            print(i, j, t)
            raise
    # print( [rho(k,t) for k in range(j,i+1)])
    return sum(gterm(i, j, k) * rho(k, t) for k in range(j, i + 1))


def lineages_through_time_nodes(
    tree: dendropy.Tree, Ne: int = 1, delta: float = 1000, epsilon: float = 1e-6
):
    children = []  # list of #nlineages:proba dictionaries
    n_edges = 0
    for n in tree.postorder_internal_node_iter():
        edges = n.child_edges()
        for e in edges:
            if e.is_leaf() or e.length > delta:
                children.append({1: 1})
        e_a, e_b = edges

        a = children.pop()
        b = children.pop()
        probs = defaultdict(lambda: 0.0)
        for ai in a:
            for bi in b:
                for i_a in range(1, ai + 1):
                    for i_b in range(1, bi + 1):
                        p_a = g(ai, i_a, e_a.length) * a[ai]
                        p_b = g(bi, i_b, e_b.length) * b[bi]
                        probs[i_a + i_b] += p_a * p_b
                    # # print(ai+bi, j, e.length,'p',p)
                    # if epsilon>0 and p<epsilon: #and p0<epsilon :
                    #     break # don't know if p's are increasing or decreasing
                    # probs[j] += p*a[ai]*b[bi]
                    # print(probs[j])
                    # p0=p
        children.append(probs)
        n_edges += 1
        # if n_edges>400:break

    probs = children.pop()
    print(probs, sum(probs.values()))
    return sum(probs[k] * (1 - 1 / k) for k in probs) * 4 * Ne


def lineages_through_time(
    tree: dendropy.Tree, Ne: int = 1, delta: float = 2000, epsilon: float = 1e-8
):
    for n in tree.preorder_internal_node_iter(exclude_seed_node=True):
        if n.edge.length > delta / (2 * Ne):
            # since n affects \sum_k=1^n P(A_t=1|A_0=k), epsilon depends on number of descendents.
            # 1 - exp( -delta / (2Ne) ) >= Pr(A_t=2 | A_0=\infty); ET\approx 4N
            tree.prune_subtree(n)
    children = []  # list of #nlineages:proba dictionaries
    n_edges = 0
    # for leaf_edge_iter
    for e in tree.postorder_edge_iter():
        if e.is_leaf():  # or e.length > delta:
            children.append({1: 1})
        else:
            a = children.pop()
            b = children.pop()

            # probs[i] = Pr(i lineages at time t)
            probs = defaultdict(lambda: 0.0)
            for ai in a:
                for bi in b:
                    p_AB = a[ai] * b[bi]
                    if epsilon > 0 and p_AB < epsilon:
                        continue
                    p = 1
                    if (
                        e.tail_node is None
                    ):  # seed node, only care about lineages entering
                        probs[ai + bi] += p_AB
                    else:
                        for j in range(1, ai + bi + 1):
                            p = g(ai + bi, j, e.length)
                            # print(ai+bi, j, e.length,'p',p)
                            if epsilon > 0 and p < epsilon:
                                continue  # p's will increase and then decrease
                            probs[j] += p * p_AB
            # s = sum(probs.values())
            # if s<.9999999999:print(probs, s)
            children.append(probs)
        n_edges += 1
    probs = children.pop()
    # print(probs, sum(probs.values()))
    return 4 * Ne * (1 - sum(probs[k] / k for k in probs)), probs

    # return sum(probs[k]*(1-1/k) for k in probs)*4*Ne, probs # identical to the above only if sum(p_k)=1


def phylogenetic_covariance(gene_tree):
    """calculate phylogenetic covariance matrix from rooted tree.
    taxa pairs appear in same order as gene_tree's *lexicographically sorted* namespace."""
    gene_tree.calc_node_root_distances()
    ns = sorted(gene_tree.taxon_namespace)
    ntaxa = len(ns)
    dmat = np.zeros((ntaxa, ntaxa))
    for t1 in range(ntaxa):
        for t2 in range(t1 + 1):
            cov = gene_tree.mrca(taxa=[ns[t1], ns[t2]]).root_distance
            dmat[t1, t2] = dmat[t2, t1] = cov
    return dmat


def calculate_uncorrected_covariance(
    tree: dendropy.Tree,
    L: int,
    Ne: int = 1,
    sig2=1,
    delta=1000,
    epsilon=1e-5,
    as_df=False,
):
    tree.calc_node_root_distances()
    # H=tree.max_distance_from_root()
    L_sig = L * sig2
    taxa = [n.taxon for n in tree.leaf_nodes()]
    covs = {}
    for a, b in combinations_with_replacement(taxa, 2):
        covs[(a.label, b.label)] = L_sig * tree.mrca(taxa=(a, b)).root_distance
        covs[(b.label, a.label)] = covs[(a.label, b.label)]
    if as_df:
        return pd.Series(covs).unstack()
    return covs


def calculate_covariance(
    tree: dendropy.Tree,
    L: int,
    Ne: int = 1,
    sig2=1,
    delta=1000,
    epsilon: float = 1e-10,
    as_df=False,
):
    """use the recursive method to calculate an approximation to the theoretical C^* matrix"""
    #     TODO: memoize somehow, since this takes a long time
    E_k, _ = lineages_through_time(tree, Ne=Ne, delta=delta, epsilon=epsilon)
    # tree.calc_node_ages()
    tree.calc_node_root_distances()
    H = tree.max_distance_from_root()
    L_sig = L * sig2
    # print(L,sig2,H,E_k,Ne)
    # if H<=Ne:print('WARNING: covariances will be negative')
    summand = E_k - 2 * Ne
    taxa = [n.taxon for n in tree.leaf_nodes()]
    covs = {}
    for a in tree.leaf_node_iter():
        covs[(a.taxon.label, a.taxon.label)] = L_sig * (
            E_k + a.root_distance
        )  # identical for ultrametric trees
    for a, b in combinations(taxa, 2):
        covs[(a.label, b.label)] = L_sig * max(
            0, summand - tree.mrca(taxa=(a, b)).root_distance
        )
        covs[(b.label, a.label)] = covs[(a.label, b.label)]
    if as_df:
        return pd.Series(covs).unstack()
    return covs


def get_branch_lengths_and_heights_numb_branches(t):
    """returns array of <branch_length, height>, with height measured in
    number of branches"""
    t.calc_node_root_distances()
    taxa = len(t)
    x = np.empty([2, taxa - 1])
    for i, n in enumerate(t.preorder_internal_node_iter(exclude_seed_node=True)):
        bl = n.edge.length
        h = n.root_distance
        x[:, i] = (bl, h)
    return x


def get_branch_lengths_and_heights(t):
    """returns array of <branch_length, dist_from_root>"""
    t.resolve_node_depths()
    taxa = len(t)
    res = []
    for i, n in enumerate(t.internal_nodes(exclude_seed_node=True)):
        bl = n.edge.length  # edge subtending n
        h = n.depth
        res.append((bl, h))
    return np.stack(res, 1)


# @timeout(TIMEOUT)
def sim_tree_covs_branches(mu: float, tree_height: float, lam: float, nsims: int = 100):
    """returns correlation of branch lengths and heights for nsims simulated trees.
    If raw is True, computes rho across all simulated trees (this makes sense since many conditions produce very small trees).  Otherwise, returns average rho across all trees.

    Args:
        mu (float): _death rate_
        tree_height (float): _description_
        lam (float): _birthrate_

    Returns:
        list: _description_
    """
    out = []
    n = 0
    bl_min = []
    while n < 100 * nsims:
        t = treesim.birth_death_tree(
            birth_rate=lam, death_rate=mu, max_time=tree_height, num_extant_tips=5000
        )
        # output of this func is inconsistent with ordering of internal_nodes iterator
        taxa = len(t)
        if taxa < 3:  # skip trees with too few branches
            continue
        x = get_branch_lengths_and_heights(t)
        if np.any(x > 1e100):
            raise ValueError(f"branch lengths too large: {t.as_string('newick')}")
        bl_min.append(x[0].min())
        out.append(x)
        n += taxa - 1

    out = np.concat(out, 1)
    corr = stats.spearmanr(out, nan_policy="omit", axis=1)
    return corr.correlation, corr.pvalue, out[0].mean(), np.mean(bl_min)


# @timeout(TIMEOUT)
def sim_tree_covs(
    mu: float, tree_height: float, lam: float, raw=False, nsims: int = 100
):
    """returns correlation of branch lengths and heights for nsims simulated trees.
    If raw is True, computes rho across all simulated trees (this makes sense since many conditions produce very small trees).  Otherwise, returns average rho across all trees.

    Args:
        mu (float): _death rate_
        tree_height (float): _description_
        lam (float): _birthrate_

    Returns:
        list: _description_
    """
    if raw:
        out = []
    else:
        c = np.empty(nsims)
        pvals = np.empty(nsims)
    sizes = np.empty(nsims)
    for j in range(nsims):
        t = treesim.birth_death_tree(
            birth_rate=lam, death_rate=mu, max_time=tree_height, num_extant_tips=4000
        )
        # output of this func is inconsistent with ordering of internal_nodes iterator
        taxa = len(t)
        x = get_branch_lengths_and_heights(t)
        if raw:
            out.append(x)
        else:
            corr = stats.spearmanr(x, axis=1)
            c[j] = corr.correlation
            pvals[j] = corr.pvalue
        sizes[j] = taxa
    #         if np.isnan(c[j]):
    #             print(np.corrcoef(x))
    #         print(taxa)
    if raw:
        out = np.concat(out, 1)
        corr = stats.spearmanr(out, axis=1)
        return corr.correlation, corr.pvalue
    return c, sizes, pvals


def sim_bd_trees(
    mu=0.5,
    lam=1,
    tree_height=6,
    n_parent_trees=100,
):
    """simulate n trees with constant lam.
    NOTE: joblib's loky backend is ok with rng."""
    ptree_params = dict(
        birth_rate=lam,
        death_rate=mu,
        max_time=tree_height,
        num_extant_tips=2000,
        repeat_until_success=True,
    )
    parent_trees = Parallel(njobs)(
        delayed(treesim.birth_death_tree)(**ptree_params) for _ in range(n_parent_trees)
    )
    return parent_trees


@timeout(TIMEOUT)
def simulate_gene_trees(
    data_dir: Path, parent_trees: list = None, n_genes: list = [100], njobs=24
):
    """For a collection of species trees, simulate n_genes and
    calculate the phylogenetic covariance SUMMED over all genes"""
    max_genes = np.max(n_genes)
    sigmas = []
    if parent_trees is None:
        parent_trees = [
            dendropy.Tree.get_from_string(s, schema="newick", rooting="force-rooted")
            for s in (data_dir / "parent_trees_4plus.phy").open("r")
        ]  # TreeList.get(path=) will create a shared namespace

    with Parallel(njobs) as parallel:
        for i, parent_tree in enumerate(parent_trees):
            gene_to_species_map = (
                dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
                    containing_taxon_namespace=parent_tree.taxon_namespace,
                    num_contained=1,
                )
            )

            gt_params = dict(
                containing_tree=parent_tree,
                gene_to_containing_taxon_map=gene_to_species_map,
            )
            gene_trees = parallel(
                delayed(treesim.contained_coalescent_tree)(**gt_params)
                for i in range(max_genes)
            )
            covs = parallel(delayed(phylogenetic_covariance)(gt) for gt in gene_trees)

            s = [np.sum(covs[:g], 0) for g in n_genes]
            sigmas.append(s)

            if not i % 50:
                #                 swap 1st and 2nd dimensions: #parent_trees x #ngenes -> #ngenes x #parent_trees
                # break
                for g, s in zip(n_genes, list(zip(*sigmas))):
                    np.savez(data_dir / f"sigmas_ngenes{g}.npz", *s)

    sigmas = list(zip(*sigmas))
    for g, s in zip(n_genes, sigmas):
        np.savez(data_dir / f"sigmas_ngenes{g}.npz", *s)
    return sigmas


def simulate_c_star_matrices(
    data_dir: Path, n_genes: list, overwrite: bool = False, resume=False, njobs=24
):
    """load OR simulate gene trees and calculate SUMMED phylogenetic covariance.
    Output shape: #n_genes x ntaxa x ntaxa"""
    n_parent_trees = len(
        (data_dir / "parent_trees_4plus_pos_edges.phy").open("r").readlines()
    )
    try:
        if overwrite:
            raise FileNotFoundError
        sigmas = []
        for g in n_genes:
            s = np.load(data_dir / f"sigmas_ngenes{g}.npz")
            if resume and len(s) < n_parent_trees:
                print(f"not enough parent trees: {len(s)} in: {data_dir}")
                raise FileNotFoundError
            sigmas.append([s[f] for f in s.files])
        # print(f'loaded sigmas from {data_dir}')
    except FileNotFoundError:
        sigmas = simulate_gene_trees(data_dir=data_dir, n_genes=n_genes, njobs=njobs)
        print(f"simulated empirical sigma to {data_dir}")
    return sigmas

# ns = sorted(t.taxon_namespace)
# indices = {k:v for v,k in enumerate(ns)}


def subsample_tree_random(t: dendropy.Tree, eps: float = 0.05):
    """subsample tree by removing all but 1 randomly chosen child node from each internal edge with length < eps.
    Continues pruning until no edges are left with length < eps.
    """
    new_tree = t.clone()
    updated = True
    while updated:
        updated = False
        if len(new_tree)<=4: return new_tree
        for e in new_tree.postorder_internal_edge_iter(
            filter_fn=None, exclude_seed_edge=True
        ):
            if e.length < eps:
                head, tail = e.head_node, e.tail_node
                candidates = list(
                    tail.child_nodes()
                )  # [*tail.child_nodes(), *head.child_nodes()] # tail towards root, head towards leaves
                v = np.random.choice(candidates)
                new_tree.prune_subtree(v)
                updated = True
    return new_tree
