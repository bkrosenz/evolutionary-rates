import json
import operator
import re
import signal
from collections import Counter, defaultdict
from functools import *
from itertools import combinations, combinations_with_replacement
from multiprocessing.context import TimeoutError
from multiprocessing.connection import wait
from pathlib import Path
from sys import argv
from time import sleep

import dendropy
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

timeout_time = 15*60

n_parent_trees = 100
n_genes = [5, 10, 25, 50]
nsims = 250
njobs = 24


def symmetric_log_det_distance(X, Y, normalize=False):
    """https://lear.inrialpes.fr/people/cherian/papers/jbld-pami_final.pdf"""
    r = np.log(linalg.det((X+Y)/2)) - np.log(linalg.det(X@Y))/2
    if normalize:
        n = len(X)
        r /= (n**2-n)/2
    return r


def frobenius_distance(X, Y, normalize=False):
    r = linalg.norm(X-Y, 'fro')
    if normalize:
        n = len(X)
        r /= (n**2-n)/2
    return r


def spectral_distance(X, Y, normalize=False):
    r = linalg.norm(X-Y, 2)
    if normalize:
        n = len(X)
        r /= (n**2-n)/2
    return r


metrics = {
    # 'ld':symmetric_log_det_distance,
    'fro': frobenius_distance,
    'spec': spectral_distance
}


@lru_cache(maxsize=300)
def rho(k, t): return rho_numba(k, t)


@nb.njit
def rho_numba(k, t, Ne=1):
    return np.exp((-k**2+k)*t/(2*Ne))


@lru_cache(maxsize=300)
def upto(a, j): return upto_numba(a, j)


@nb.njit(parallel=True)
def upto_numba(a, j):
    if j < 1:
        return 1
    r = 1
    for k in nb.prange(a, a+j):
        r *= k
    return r


@lru_cache(maxsize=300)
def downto(a, j): return downto_numba(a, j)


@nb.njit(parallel=True)
def downto_numba(a, j):
    if j < 1:
        return 1
    r = 1
    for k in nb.prange(a-j+1, a+1):
        r *= k
    return r


@lru_cache(maxsize=300)
def lupto(a, j): return lupto_numba(a, j)


@nb.njit
def lupto_numba(a, j):
    if j < 1:
        return 0
    return np.sum(np.log(np.arange(a, a+j)))


@lru_cache(maxsize=300)
def ldownto(a, j): return ldownto_numba(a, j)


@nb.njit
def ldownto_numba(a, j):
    if j < 1:
        return 0
    return np.sum(np.log(np.arange(a, a-j, -1)))

# @njit(parallel=True)


@lru_cache(maxsize=400)
def gterm(i, j, k):
    if not j <= i:
        raise ValueError(f'i={i},j={j}')

    if j > 1:
        sgn = (k-j) % 2 and -1 or 1
    else:
        sgn = k % 2 and -1 or 1

    if i < 15:
        try:
            res = (2*k-1) * downto(i, k) / upto(i, k)
            if j > 1:
                res *= upto(j, k-1) / (factorial(j) * factorial(k-j))
            # print(i,j,k,res)
            return res * sgn
        except:
            print(i, j, k)
            raise
    else:
        try:
            res = np.log(2*k-1) + ldownto(i, k) - lupto(i, k)
            if j > 1:
                res += lupto(j, k-1) - gammaln(j+1) - gammaln(k-j+1)
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
            return 1 - sum(gterm(i, j, k)*rho(k, t) for k in range(2, i+1))
        except:
            print(i, j, t)
            raise
    # print( [rho(k,t) for k in range(j,i+1)])
    return sum(gterm(i, j, k)*rho(k, t) for k in range(j, i+1))


def lineages_through_time_nodes(tree: dendropy.Tree,
                                Ne: int = 1,
                                delta: float = 1000,
                                epsilon: float = 1e-6):
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
        probs = defaultdict(lambda: 0.)
        for ai in a:
            for bi in b:
                for i_a in range(1, ai+1):
                    for i_b in range(1, bi+1):
                        p_a = g(ai, i_a, e_a.length) * a[ai]
                        p_b = g(bi, i_b, e_b.length) * b[bi]
                        probs[i_a+i_b] += p_a * p_b
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
    return sum(probs[k]*(1-1/k) for k in probs)*4*Ne


def lineages_through_time(tree: dendropy.Tree,
                          Ne: int = 1,
                          delta: float = 2000,
                          epsilon: float = 1e-8):
                        #   Assumes population of size 2*Ne
    for n in tree.preorder_internal_node_iter(exclude_seed_node=True):
        
        if n.edge.length > delta/(2*Ne):
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
            probs = defaultdict(lambda: 0.)
            for ai in a:
                for bi in b:
                    p_AB = a[ai]*b[bi]
                    if epsilon > 0 and p_AB < epsilon:
                        continue
                    p = 1
                    if e.tail_node is None:  # seed node, only care about lineages entering
                        probs[ai+bi] += p_AB
                    else:
                        for j in range(1, ai+bi+1):
                            p = g(ai+bi, j, e.length)
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
    return 4*Ne*(1-sum(probs[k]/k for k in probs)), probs

# return sum(probs[k]*(1-1/k) for k in probs)*4*Ne


def phylogenetic_covariance(tree):
    """calculate phylogenetic covariance matrix from a rooted (species) tree.
    taxa pairs appear in same order as gene_tree's *lexicographically sorted* namespace."""
    if not tree.is_rooted:
        raise ValueError('tree must be rooted')

    tree.calc_node_root_distances()
    ns = sorted(tree.taxon_namespace)
    ntaxa = len(ns)
    dmat = np.zeros((ntaxa, ntaxa))
    for t1 in range(ntaxa):
        for t2 in range(t1+1):
            cov = tree.mrca(taxa=[ns[t1], ns[t2]]).root_distance
            dmat[t1, t2] = dmat[t2, t1] = cov
    return dmat


def calculate_uncorrected_covariance(tree: dendropy.Tree,
                                     L: int,
                                     Ne: int = 1,
                                     sig2=1,
                                     delta=1000,
                                     epsilon=1e-5,
                                     as_df=False):

    tree.calc_node_root_distances()
    # H=tree.max_distance_from_root()
    L_sig = L*sig2
    taxa = [n.taxon for n in tree.leaf_nodes()]
    covs = {}
    for a, b in combinations_with_replacement(taxa, 2):
        covs[(a.label, b.label)] = L_sig*tree.mrca(taxa=(a, b)).root_distance
        covs[(b.label, a.label)] = covs[(a.label, b.label)]
    if as_df:
        return pd.Series(covs).unstack()
    return covs


def calculate_covariance(tree: dendropy.Tree,
                         L: int,
                         Ne: int = 1,
                         sig2=1,
                         delta=1000,
                         epsilon=1e-10,
                         as_df=False):
    """calculate C*, the expected coalescent covariance matrix

    Args:
        tree (dendropy.Tree): species
        L (int): number of loci
        Ne (int, optional): pop size (2Ne will be used for all calcs). Defaults to 1.
        sig2 (int, optional): variance of each locus. Defaults to 1.
        delta (int, optional): branch length cutoff. Defaults to 1000.
        epsilon (_type_, optional): probability. Defaults to 1e-10.
        as_df (bool, optional): return df. Defaults to False.

    Returns:
        _type_: _description_
    """
    if not tree.is_rooted:
        raise ValueError('tree must be rooted')
    E_k, _ = lineages_through_time(tree, Ne=Ne, delta=delta, epsilon=epsilon)
    # tree.calc_node_ages()
    tree.calc_node_root_distances()
    H = tree.max_distance_from_root()
    L_sig = L*sig2
    # print(L,sig2,H,E_k,Ne)
    # if H<=Ne:print('WARNING: covariances will be negative')
    summand = E_k - 2*Ne
    taxa = [n.taxon for n in tree.leaf_nodes()]
    covs = {}
    for a in tree.leaf_node_iter():
        covs[(a.taxon.label, a.taxon.label)] = L_sig * \
            (E_k + a.root_distance)  # identical for ultrametric trees
    for a, b in combinations(taxa, 2):
        covs[(a.label, b.label)] = L_sig * \
            max(0, summand + tree.mrca(taxa=(a, b)).root_distance)
        covs[(b.label, a.label)] = covs[(a.label, b.label)]
    if as_df:
        return pd.Series(covs).unstack()
    return covs


def get_branch_lengths_and_heights(t):
    t.calc_node_ages()
    taxa = len(t)
    root_age = next(t.ageorder_node_iter(descending=True)).age
    x = np.empty([2, taxa-1])
    for i, n in enumerate(t.internal_nodes()):
        bl = n.edge.length
        h = root_age-n.age
        x[:, i] = (bl, h)
    return x


@timeout(timeout_time)
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
            birth_rate=lam,
            death_rate=mu,
            max_time=tree_height,
            ntax=4000)
        # output of this func is inconsistent with ordering of internal_nodes iterator
        taxa = len(t)
        x = get_branch_lengths_and_heights(t)
        c[j] = stats.spearmanr(x, axis=1).correlation
        sizes[j] = taxa
#         if np.isnan(c[j]):
#             print(np.corrcoef(x))
#         print(taxa)
    return c, sizes


def sim_bd_trees(
    mu=.5,
    lam=1,
    tree_height=6,
    n_parent_trees=100,
):
    """simulate n trees with constant lam.
    NOTE: joblib's loky backend is ok with rng."""
    ptree_params = dict(birth_rate=lam,
                        death_rate=mu,
                        max_time=tree_height,
                        num_extant_tips=2000,
                        repeat_until_success=True)
    parent_trees = Parallel(njobs)(
        delayed(treesim.birth_death_tree)(
            **ptree_params) for _ in range(n_parent_trees))
    return parent_trees


def simulate_gene_trees(data_dir: Path,
                        parent_trees: list = None,
                        n_genes: list = [100],
                        njobs=24):
    """For a collection of species trees, simulate n_genes and 
    calculate the phylogenetic covariance SUMMED over all genes"""
    max_genes = np.max(n_genes)
    sigmas = []
    if parent_trees is None:
        parent_trees = [
            dendropy.Tree.get_from_string(s, schema="newick", rooting="force-rooted") for s in (
                data_dir/'parent_trees_4plus_pos_edges.phy').open('r')]  # TreeList.get(path=) will create a shared namespace
    
    with Parallel(njobs) as parallel:
        for i, parent_tree in enumerate(parent_trees):
            if isinstance(parent_trees[0],str):
                parent_tree=dendropy.Tree.get_from_string(parent_tree, schema="newick", rooting="force-rooted")

            gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
                containing_taxon_namespace=parent_tree.taxon_namespace,
                num_contained=1)

            gt_params = dict(containing_tree=parent_tree,
                             gene_to_containing_taxon_map=gene_to_species_map)
            # units of 2Ne
            gene_trees = parallel(
                delayed(treesim.contained_coalescent_tree)(
                    **gt_params) for i in range(max_genes))
            covs = parallel(
                delayed(phylogenetic_covariance)(gt) for gt in gene_trees)

            s = [np.sum(covs[:g], 0) for g in n_genes]
            sigmas.append(s)

            if not i % 50:
                #                 swap 1st and 2nd dimensions: #parent_trees x #ngenes -> #ngenes x #parent_trees
                # break
                for g, s in zip(n_genes, list(zip(*sigmas))):
                    np.savez(data_dir/f'sigmas_ngenes{g}.npz', *s)

    sigmas = list(zip(*sigmas))
    for g, s in zip(n_genes, sigmas):
        np.savez(data_dir/f'sigmas_ngenes{g}.npz', *s)
    return sigmas


def simulate_c_star_matrices(data_dir: Path,
                             n_genes: list,
                             parent_tree_file:Path=None,
                             overwrite: bool = False,
                             resume=False,
                             njobs=24):
    """load OR simulate gene trees and calculate SUMMED phylogenetic covariance.
    Output shape: #n_genes x ntaxa x ntaxa """
    parent_trees=(parent_tree_file if parent_tree_file else data_dir/'parent_trees_4plus_pos_edges.phy').open('r').readlines()
    
    n_parent_trees = len(
        parent_trees)
    try:
        if overwrite:
            raise FileNotFoundError
        sigmas = []
        for g in n_genes:
            s = np.load(data_dir/f'sigmas_ngenes{g}.npz')
            if resume and len(s) < n_parent_trees:
                print(f'not enough parent trees: {len(s)} in: {data_dir}')
                raise FileNotFoundError
            sigmas.append([s[f] for f in s.files])
        # print(f'loaded sigmas from {data_dir}')
    except FileNotFoundError:
        sigmas = simulate_gene_trees(
            data_dir=data_dir, 
            parent_trees=parent_trees,
            n_genes=n_genes,
              njobs=njobs)
        print(f'simulated empirical sigma to {data_dir}')
    return sigmas



@timeout(timeout_time)
def simulate_traits(n_genes: list,
                    sigmas: list, # list of lists
                    data_dir: Path,
                    nsamps=100,
                    njobs:int=24,
                    overwrite:bool=False):
    """Input: sigmas: #ngenes x #parent_trees x ntaxa x ntaxa.
    NOTE: values in n_genes list should match ngenes parameter used to generate sigmas list"""

    for g, sig in zip(n_genes, sigmas):
        try:
            if overwrite:
                raise FileNotFoundError
            samples = list(np.load(data_dir/f'samples_ngenes{g}.npz').values())
        except:
            samples = Parallel(njobs)(
                delayed(multivariate_normal.rvs)
                (np.zeros(len(s)), s, nsamps, i) for i, s in enumerate(sig)
            )
            # samples={f'n{n}':s for n,s in zip(sig,samples)}
            np.savez(data_dir/f'samples_ngenes{g}.npz', *samples)
            print(f'simulated {g}-gene  traits to {data_dir}')
    return samples


if __name__ == '__main__':
    njobs = int(argv[1])
    if len(argv)>2:
        parent_dir = Path(argv[2])
    else:
        parent_dir = Path('/N/project/phyloML/rate_timescaling/data/tree_sims')
    # error :
#      ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
#   File "/N/project/phyloML/app/Quartz/miniconda3/envs/p310/lib/python3.11/site-packages/joblib/_parallel_backends.py", line 569, in wrap_future_result
#     raise TimeoutError from e
# multiprocessing.context.TimeoutError

    # for h in range(1, 9, 2):
    #     for d in parent_dir.glob(f'*/*/h{h}'):
    #         print(d)
    #         try:
    #             sigmas = simulate_c_star_matrices(
    #                 d, n_genes=n_genes,
    #                 overwrite=False,
    #                   njobs=njobs)
    #             samples = simulate_traits(
    #                 n_genes, sigmas,
    #                 data_dir=d,
    #                 njobs=njobs, overwrite=False)
    #         except TimeoutError:
    #             print('timeout, continuing...')
    #             continue

    dists_dict = {}  # defaultdict(dict)

    
#    for h in [5, 1, 3, 7]:
#        for i, d in enumerate(Path('/N/project/phyloML/rate_timescaling/data/tree_sims').rglob(f'*/*/h{h}')):
#            new_dir = d/'no_ils'
#            new_dir.mkdir(exist_ok=True)
#            if (new_dir / 'samples_ngenes1.npz').exists():continue
#            try:
#                with open(d/'parent_trees_4plus_pos_edges.phy','r') as f:
#                    parent_trees = [dendropy.Tree.get_from_string(t,'newick') for t in  f]
#                sigmas = *map(phylogenetic_covariance, parent_trees),
#                new_dir = d/'no_ils'
#                new_dir.mkdir(exist_ok=True)
#                print(d)
#                samples = simulate_traits([1], [sigmas], data_dir=new_dir, overwrite=False) #parent_trees x #nsamples x ntaxa
#            except FileNotFoundError:
#                continue
#    print('finished all control sims')

    for h in [5, 1, 3, 7]:
        for i, d in enumerate(parent_dir.rglob(f'*/*/h{h}')):
            dsamp = d/'samples_ngenes25.npz'
            dsig = d/'sigmas_ngenes25.npz'
            if dsamp.exists() and dsig.exists() and dsamp.stat().st_ctime < dsig.stat().st_ctime:
                dsamp.unlink()
            # if not any(d.glob('samp*npz')):  # don't simulate
            if not (d/'parent_trees_4plus_pos_edges.phy').exists():
                continue
            print('in', d)

            with open(d/'parent_trees_4plus_pos_edges.phy', 'r') as f:
                parent_trees = [
                    dendropy.Tree.get_from_string(
                        t, 'newick', rooting="force-rooted") for t in f]
            tree_sizes = list(map(len, parent_trees))
            if stats.mode(tree_sizes).mode[0] == 2000:
                continue
            try:
                dcov=d/'covs.gz'
                if dcov.stat().st_ctime < dsig.stat().st_ctime:
                    raise FileNotFoundError
                covs = joblib.load(dcov)
                heights = []
                for tree in parent_trees:
                    heights.append(tree.max_distance_from_root())
                # raise
            except FileNotFoundError:
                covs = []
                heights = []
                for tree in parent_trees:
                    c = calculate_covariance(
                        tree, L=1,
                        Ne=0.5,
                        epsilon=1e-8,
                        delta=200,
                        as_df=True)
                    # if c.max().max() == 0:
                    #     break
                    covs.append(c)
                    heights.append(tree.max_distance_from_root())
                joblib.dump(covs, d/'covs.gz')

            except KeyError as e:
                print(d)    
                raise e
                # print('computed covariances')
            # continue

            covs_u = []
            for tree in parent_trees:
                c = phylogenetic_covariance(
                    tree)
                if c.max().max() == 0:
                    break
                covs_u.append(c)
            joblib.dump(covs_u, d/'covs_u.gz')

            try:

                sigmas = simulate_c_star_matrices(
                    d, n_genes=n_genes, overwrite=False)

                # parent_trees x #nsamples x ntaxa
                samples = simulate_traits(
                    n_genes, sigmas, data_dir=d, overwrite=False)

                # print(len(sigmas),len(n_genes))
                for ngenes, sigs in zip(n_genes, sigmas):
                    for tree_ix, (h, c, c_u, s) in enumerate(zip(heights, covs, covs_u, sigs)):
                        n = len(c)
                        for fname, f in metrics.items():
                            dists_dict[('c', n, h, fname, ngenes, d.parent.name[1:], d.parent.parent.name[1:])] = f(
                                c_u, s/ngenes, normalize=True)
                            dists_dict[('c_star', n, h, fname, ngenes, d.parent.name[1:], d.parent.parent.name[1:])] = f(
                                c, s/ngenes, normalize=True)
                dists = pd.Series(dists_dict)  # .unstack(level=[-1,-2])
                dists.index.names = ['cov_type', 'n',
                                     'h', 'metric', 
                                     'ngenes', 'mu', 'lam']
                dists = dists.reset_index().rename(columns={0: 'value'})
                dists[['mu', 'lam']] = dists[['mu', 'lam']].applymap(float)
                dists['r'] = dists.lam-dists.mu
                dists['rho'] = dists.lam/dists.mu
                dists['val_h'] = dists.value / dists.h
                dists.to_pickle(parent_dir / f'c_pred_v_actual_{d.name}.pd.gz')
            except TimeoutError:
                continue

        sns.relplot(
            # .query('ngenes==25'), #.query('metric=="fro" & ngenes==25'),
            data=dists.query('cov_type=="c"'),
            col='metric',
            hue='lam',
            style='lam',
            row='ngenes',
            x='r',
            y='value',
            kind='line',
        )
        plt.ylim(0, .7)
        plt.tight_layout()
        plt.savefig(parent_dir / f'c_v_actual_{d.name}.png')
        plt.clf()

        sns.relplot(
            # .query('ngenes==25'), #.query('metric=="fro" & ngenes==25'),
            data=dists.query('cov_type=="c_star"'),
            col='metric',
            hue='lam',
            style='lam',
            row='ngenes',
            x='r',
            y='value',
            kind='line')
        plt.ylim(0, .7)

        plt.tight_layout()
        plt.savefig(parent_dir / f'c_star_v_actual_{d.name}.png')
        plt.clf()


        v=(dists
           .query('ngenes==10')
           .groupby(['h','metric', 'ngenes', 'mu', 'lam','r']))
        rd=(v
                 .apply(
                     lambda s:np.mean(s.query('cov_type=="c"').val_h.values-s.query('cov_type=="c_star"').val_h.values)
                     )
                     .reset_index()
                     )
        rd=rd.rename(columns={0:'rel_dif'})


        rel=sns.relplot(
            # .query('ngenes==25'), #.query('metric=="fro" & ngenes==25'),
            data=rd,
            col='metric',
            hue='lam',
            style='lam',
            row='ngenes',
            x='r',
            y='rel_dif',
            kind='line')
        # plt.ylim(0, .7)
        rel.fig.suptitle(r'$d(C,\Sigma)-d(C^*,\Sigma)$')

        plt.tight_layout()
        plt.savefig(parent_dir / f'c_m_c_star_{d.name}.png')
        plt.clf()


        rd=(v
                 .apply(
                     lambda s:np.mean(s.query('cov_type=="c_star"').val_h.values/s.query('cov_type=="c"').val_h.values)
                     )
                     .reset_index()
                     )
        rd=rd.rename(columns={0:'rel_dif'})

        
        rel=sns.relplot(
            # .query('ngenes==25'), #.query('metric=="fro" & ngenes==25'),
            data= rd,
            col='metric',
            hue='lam',
            style='lam',
            row='ngenes',
            x='r',
            y='rel_dif',
            kind='line')
        # plt.ylim(0, .7)
        rel.fig.suptitle(r'$d(C^*,\Sigma)/d(C,\Sigma)$')

        plt.tight_layout()
        plt.savefig(parent_dir / f'c_v_c_star_{d.name}.png')
        plt.clf()
