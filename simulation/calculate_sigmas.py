import json
from importlib import reload
import ete3
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
from time import sleep

import dendropy
import joblib  # import dump,load
# import numba as nb

# import msprime
import numpy as np
import pandas as pd
from dendropy.simulate import treesim
from joblib import Parallel, delayed
from scipy import linalg, stats
from scipy.special import factorial, gammaln
from scipy.stats import multivariate_normal, hmean


import re

rx = re.compile("\/l([\d.]+)\/m([\d.]+)\/")

nsims = 250
lam = 1
mus = np.arange(0.1, 1, 0.1)
heights = np.arange(1, 8, 0.5)
njobs = 24


parent_dir = Path("/N/project/phyloML/rate_timescaling/data/tree_sims")
if not parent_dir.exists():
    parent_dir = Path("/home/b/rate_timescaling/tree_sims")
figdir = parent_dir.parent / "figures"

# mu=.8
# tree_height=9

n_parent_trees = 100
n_parent_trees = 100
nsims = 250

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


def get_branch_lengths_and_heights(t):
    """returns array of <branch_length, height>"""
    t.calc_node_root_distances()
    taxa = len(t)
    x = np.empty([2, taxa - 1])
    for i, n in enumerate(t.preorder_internal_node_iter(exclude_seed_node=True)):
        bl = n.edge.length
        h = n.root_distance
        x[:, i] = (bl, h)
    return x


def get_branch_lengths_and_heights2(t):
    """returns array of <branch_length, height>"""
    t.calc_node_ages()
    taxa = len(t)
    root_age = next(t.ageorder_node_iter(descending=True)).age
    x = np.empty([2, taxa - 1])
    for i, n in enumerate(t.internal_nodes(exclude_seed_node=True)):
        bl = n.edge.length
        h = root_age - n.age
        x[:, i] = (bl, h)
    return x


def subsample_tree_random(t: dendropy.Tree, eps: float = 0.05):
    """subsample tree by removing all but 1 randomly chosen child node from each internal edge with length < eps.
    Continues pruning until no edges are left with length < eps.
    """
    new_tree = t.clone()
    updated = True
    while updated:
        updated = False
        if len(new_tree) <= 4:
            return new_tree
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


ngenes = 1


def estimate_sigmas(
    t,
    chars: Union[np.array, dendropy.ContinuousCharacterMatrix],
    pic_calculator=continuous.PhylogeneticIndependentContrasts,
    max_replicates: int = None,
):
    n_replicates, n_taxa = chars.shape
    if max_replicates:
        n_replicates = min(n_replicates, max_replicates)

    for e in t.edges():
        if e.length == 0:
            e.length = 1e-9
    if not isinstance(chars, dendropy.ContinuousCharacterMatrix):
        taxa = sorted(t.taxon_namespace)
        chars = dict(zip(taxa, map(list, chars.T)))
        # print(chars)
        chars = dendropy.ContinuousCharacterMatrix.from_dict(chars)

    pic = pic_calculator(tree=t, char_matrix=chars)
    sigmas = np.empty(n_replicates)
    for i in range(n_replicates):
        ctree = pic.contrasts_tree(
            character_index=i,
            annotate_pic_statistics=True,
            state_values_as_node_labels=False,
            corrected_edge_lengths=False,
        )
        sigmas[i] = sigma_pic(ctree)
    return sigmas


def subtree_sigma(t: Tree, chars: np.array, eps=0.01, n_samples=20):
    samps = set()
    sigs = []
    if not isinstance(chars, dendropy.ContinuousCharacterMatrix):
        taxa = sorted(t.taxon_namespace)
        chars = dict(zip(taxa, map(list, chars.T)))
        # print(chars)
        chars = dendropy.ContinuousCharacterMatrix.from_dict(chars)

    for _ in range(n_samples):
        subtree = subsample_tree_random(t, eps)
        sigs.append(estimate_sigmas(subtree, chars))
    return sigs


def subsample_tree(t, eps=0.05):
    # TODO: need to fix recursion
    subtrees = set()
    above_threshold = True
    for e in t.postorder_internal_edge_iter(filter_fn=None, exclude_seed_edge=True):
        if e.length < eps:
            head, tail = e.head_node, e.tail_node
            for v in tail.child_node_iter():  # tail towards root, head towards leaves
                if v != head:
                    new_tree = t.clone()
                    new_tree.prune_subtree(v)
                    subtrees.update(subsample_tree(new_tree, eps))
            for v in head.child_node_iter():  # tail towards root, head towards leaves
                new_tree = t.clone()
                new_tree.prune_subtree(v)
                subtrees.update(subsample_tree(new_tree, eps))
            above_threshold = False
    if above_threshold:
        subtrees.add(t)
    return subtrees


def get_cherries(t):
    # TODO: make this more efficient
    filter_fn = lambda n: all(child.is_leaf() for child in n.child_node_iter())
    for i in t.postorder_internal_node_iter(
        filter_fn=filter_fn, exclude_seed_node=True
    ):
        yield i


def sigma_pic(ctree):
    """
    Computes the average of the squared contrasts divided by the variance.
    Leaves do not have contrasts (requires at least a cherry), so iterate over internal nodes only.
    must call _get_contrasts on ctree first."""
    return np.mean(
        [
            (nd.pic[0]["pic_contrast_raw"] ** 2 / nd.pic[0]["pic_contrast_variance"])
            for nd in ctree.postorder_internal_node_iter()
        ]
    )
    # return np.mean([nd.pic_contrast_standardized**2 for nd in ctree.postorder_internal_node_iter()])


class BrownianBridgePIC(continuous.PhylogeneticIndependentContrasts):
    def _get_contrasts(self, character_index):
        """
        Main work-horse method. If needed, adds an entry to
        self._character_constrants, with key being the character index, and a
        value being another dictionary that contains the constrast information.
        This second dictionary has the node's id as a key and as a value the a
        dictionary with the following:

                - ``pic_state_value``
                - ``pic_state_variance``
                - ``pic_contrast_raw``
                - ``pic_contrast_variance``
                - ``pic_contrast_standardized``
                - ``pic_edge_length_error``
                - ``pic_corrected_edge_length``

        """
        #         TODO: rewrite the state, var funcs to handle polytomies
        if character_index in self._character_contrasts:
            return self._character_contrasts[character_index]
        all_results = {}
        for nd in self._tree.postorder_node_iter():
            nd_results = {}
            child_nodes = nd.child_nodes()
            if len(child_nodes) == 0:
                nd_results["pic_state_value"] = self._char_matrix[nd.taxon][
                    character_index
                ]
                nd_results["pic_state_variance"] = None
                nd_results["pic_contrast_raw"] = None
                nd_results["pic_contrast_variance"] = None
                nd_results["pic_contrast_standardized"] = None
                nd_results["pic_edge_length_error"] = 0.0
                nd_results["pic_corrected_edge_length"] = nd.edge.length
            elif len(child_nodes) == 1:
                # root node?
                nd_results["pic_state_value"] = None
                nd_results["pic_state_variance"] = None
                nd_results["pic_contrast_raw"] = None
                nd_results["pic_contrast_variance"] = None
                nd_results["pic_contrast_standardized"] = None
                nd_results["pic_edge_length_error"] = None
                nd_results["pic_corrected_edge_length"] = None
            else:
                state_vals = []
                corrected_edge_lens = []
                actual_edge_lens = []
                for cnd in child_nodes:
                    state_vals.append(all_results[cnd._track_id]["pic_state_value"])
                    actual_edge_lens.append(cnd.edge.length)
                    if (
                        all_results[cnd._track_id]["pic_corrected_edge_length"]
                        is not None
                    ):
                        corrected_edge_lens.append(
                            all_results[cnd._track_id]["pic_corrected_edge_length"]
                        )
                    else:
                        corrected_edge_lens.append(cnd.edge.length)
                n = len(state_vals)
                # numerator_fn = lambda i : (1.0/corrected_edge_lens[i]) * state_vals[i]
                # denominator_fn = lambda i  : 1.0/corrected_edge_lens[i]
                # nd_results['pic_state_value'] = \
                #         sum(numerator_fn(i) for i in range(n)) \
                #         / sum(denominator_fn(i) for i in range(n))
                sum_of_child_edges = sum(corrected_edge_lens)
                prod_of_child_edges = reduce(operator.mul, corrected_edge_lens)

                numerator = (
                    corrected_edge_lens[0] ** 2 * state_vals[1]
                    + corrected_edge_lens[1] ** 2 * state_vals[0]
                )
                denominator = sum_of_child_edges**2 - prod_of_child_edges

                nd_results["pic_state_value"] = numerator / denominator
                e0, e1 = corrected_edge_lens
                # where did this come from?
                nd_results["pic_edge_length_error"] = (
                    prod_of_child_edges * (e0**3 + e1**3)
                ) / (e0**2 + prod_of_child_edges + e1**2) ** 2
                if nd.edge.length is not None:
                    nd_results["pic_corrected_edge_length"] = (
                        nd.edge.length + nd_results["pic_edge_length_error"]
                    )
                else:
                    nd_results["pic_corrected_edge_length"] = None
                nd_results["pic_state_variance"] = nd_results[
                    "pic_corrected_edge_length"
                ]

                if len(child_nodes) != 2:
                    if self._polytomy_strategy == "ignore":
                        nd_results["pic_contrast_raw"] = None
                        nd_results["pic_contrast_standardized"] = None
                        nd_results["pic_contrast_variance"] = sum_of_child_edges
                    else:
                        raise ValueError("Tree is not fully-bifurcating")
                else:
                    nd_results["pic_contrast_raw"] = state_vals[0] - state_vals[1]
                    nd_results["pic_contrast_standardized"] = nd_results[
                        "pic_contrast_raw"
                    ] / (sum_of_child_edges**0.5)
                    nd_results["pic_contrast_variance"] = sum_of_child_edges

            nd._track_id = id(nd)  # will get cloned
            all_results[nd._track_id] = nd_results
            try:
                nd.pic[character_index] = dict(nd_results)
            except AttributeError:
                nd.pic = {character_index: dict(nd_results)}
        self._character_contrasts[character_index] = dict(all_results)
        return self._character_contrasts[character_index]


class JamesSteinPIC(continuous.PhylogeneticIndependentContrasts):
    def _get_contrasts(self, character_index):
        """
        Main work-horse method. If needed, adds an entry to
        self._character_constrants, with key being the character index, and a
        value being another dictionary that contains the constrast information.
        This second dictionary has the node's id as a key and as a value the a
        dictionary with the following:

                - ``pic_state_value``
                - ``pic_state_variance``
                - ``pic_contrast_raw``
                - ``pic_contrast_variance``
                - ``pic_contrast_standardized``
                - ``pic_edge_length_error``
                - ``pic_corrected_edge_length``
        """
        # TODO: fix numerator, denominator -
        # need to compute y_mean, sigma_hat from all cherries BEFORE we iterate postorder
        if character_index in self._character_contrasts:
            return self._character_contrasts[character_index]
        all_results = {}
        for (
            nd
        ) in (
            self._tree.postorder_node_iter()
        ):  # this height property is max num branches to a leaf
            if nd.is_leaf():
                nd.height = 0
            else:
                nd.height = max(x.height for x in nd.child_nodes()) + 1

        dists = self._tree.calc_node_root_distances(
            return_leaf_distances_only=True
        )  # use to weight leaves if tree isnt ultrametric
        root_state_est = np.average(
            [
                self._char_matrix[nd.taxon][character_index]
                for nd in self._tree.leaf_node_iter()
            ],
            weights=dists,
        )  # x_bar
        i = 0

        node_list = sorted(
            self._tree.postorder_node_iter(), key=lambda nd: nd.height
        )  # need to preserve order within level
        sig_sq = sig_sq_num = nodes_seen_so_far = 0
        states = []
        h = []
        edge_sums = []
        level_nodes = []
        prev_height = 0
        for nd_num, nd in enumerate(node_list):  # node_list gets shorter
            if nd.height != prev_height and prev_height > 0 and len(level_nodes) > 3:
                # now we go through the previous level and recompute all the state values:
                states = np.array(
                    states
                )  # TODO: store everything in a size(tree) array and index appropriately
                if not np.allclose(
                    states, root_state_est
                ):  # avoid the weird situation where everything is close to the root
                    h = np.array(h)
                    # sig_sq /= n_nodes_this_level  # nodes in PREVIOUS level
                    # y = states / h / 2 # old
                    y = states  # new
                    denom = y - root_state_est
                    n_level_nodes = len(states)  # nodes in THIS level
                    if (
                        prev_height == 1
                    ):  # could update sig_sq for deeper nodes, but for now just use cherries
                        sig_sq = sig_sq_num / nodes_seen_so_far

                    sig_sq = sig_sq_num / nodes_seen_so_far
                    num = sig_sq * (n_level_nodes - 2)
                    # states = h * 2 * (
                    #     (1 - num / np.linalg.norm(denom))  * denom + root_state_est
                    # ) # old

                    edge_sums = np.array(edge_sums)
                    # positive part JS estimator
                    states = (
                        np.maximum(0, 1 - num / np.sum(denom**2 / h)) * denom
                        + root_state_est
                    )

                    for i, n in enumerate(
                        level_nodes
                    ):  # iterate again and reassign the PIC values
                        all_results[n._track_id]["pic_state_value"] = states[i]
                        # nd_results['pic_edge_length_error'] = # TODO calculate this for JS estimator
                # RESET EVERYTHING FOR NEXT LEVEL:

                states = []
                edge_sums = []
                h = []
                level_nodes = []
                prev_height = nd.height
            # continue on with the new level

            nd_results = {}
            child_nodes = nd.child_nodes()
            if len(child_nodes) == 0:
                nd_results["pic_state_value"] = self._char_matrix[nd.taxon][
                    character_index
                ]
                nd_results["pic_state_variance"] = None
                nd_results["pic_contrast_raw"] = None
                nd_results["pic_contrast_variance"] = None
                nd_results["pic_contrast_standardized"] = None
                nd_results["pic_edge_length_error"] = 0.0
                nd_results["pic_corrected_edge_length"] = nd.edge.length
                # states.append(nd_results["pic_state_value"])
            elif len(child_nodes) == 1:
                # root node?
                nd_results["pic_state_value"] = None
                nd_results["pic_state_variance"] = None
                nd_results["pic_contrast_raw"] = None
                nd_results["pic_contrast_variance"] = None
                nd_results["pic_contrast_standardized"] = None
                nd_results["pic_edge_length_error"] = None
                nd_results["pic_corrected_edge_length"] = None
            elif len(child_nodes) > 2:
                if self._polytomy_strategy == "ignore":
                    nd_results["pic_contrast_raw"] = None
                    nd_results["pic_contrast_standardized"] = None
                    nd_results["pic_contrast_variance"] = sum_of_child_edges
            else:  # ignore polytomies
                corrected_edge_lens = []
                actual_edge_lens = []
                state_vals = []
                for cnd in child_nodes:
                    state_vals.append(all_results[cnd._track_id]["pic_state_value"])
                    actual_edge_lens.append(cnd.edge.length)
                    if (
                        all_results[cnd._track_id]["pic_corrected_edge_length"]
                        is not None
                    ):
                        corrected_edge_lens.append(
                            all_results[cnd._track_id]["pic_corrected_edge_length"]
                        )
                    else:
                        corrected_edge_lens.append(cnd.edge.length)
                state_vals = np.array(state_vals)
                edges = np.array(corrected_edge_lens)
                num = (state_vals[1] - state_vals[0]) ** 2
                sum_of_child_edges = edges.sum()
                prod_of_child_edges = edges.prod()
                # use this as denominator of running est of sig_sq
                nodes_seen_so_far += 1
                sig_sq_num += (
                    num / sum_of_child_edges
                )  # we've got a new contrast.  TODO: use new raw contrasts after JS update
                nd_results["pic_edge_length_error"] = prod_of_child_edges / (
                    sum_of_child_edges
                )
                if nd.edge.length is not None:
                    # TODO: figure out what the variance of the JS est is, use that for edge length error.
                    nd_results["pic_corrected_edge_length"] = (
                        nd.edge.length + nd_results["pic_edge_length_error"]
                    )
                nd_results["pic_contrast_raw"] = state_vals[0] - state_vals[1]
                nd_results["pic_contrast_standardized"] = nd_results[
                    "pic_contrast_raw"
                ] / (sum_of_child_edges**0.5)
                nd_results["pic_contrast_variance"] = sum_of_child_edges
                nd_results["pic_state_value"] = np.sum(state_vals / edges) / np.sum(
                    1.0 / edges
                )
                states.append(nd_results["pic_state_value"])
                edge_sums.append(sum_of_child_edges)
                h.append(prod_of_child_edges / sum_of_child_edges)
                level_nodes.append(
                    nd
                )  # only nodes with 2 children are included in JS update

            nd._track_id = id(nd)  # will get cloned

            all_results[nd._track_id] = nd_results
            try:
                nd.pic[character_index] = dict(nd_results)
            except AttributeError:
                nd.pic = {character_index: dict(nd_results)}
            prev_height = nd.height
        self.sigma_sq = np.mean(
            [
                nd["pic_contrast_standardized"] ** 2
                for nd in all_results.values()
                if nd["pic_contrast_standardized"] is not None
            ]
        )
        self._character_contrasts[character_index] = dict(all_results)
        return self._character_contrasts[character_index]


class PIC(object):
    def __init__(
        self,
        t,
        chars: np.array,
        eps: float = None,
        coalescent_correction: bool = False,
        calculate_pic=continuous.PhylogeneticIndependentContrasts,
        max_replicates=20,
    ):
        """chars: n_replicates x n_taxa.
        If eps>0, will prune cherries / paths / branches with separation < eps.
        If coalescent_correction s True, adds the *expected* time to coalescence,
        assuming branch lengths are in coalescent units.  No effect for 'full' algorithm
        """
        for e in t.edges():
            if e.length == 0:
                e.length = 1e-9
        self.t = t
        self.eps = eps
        self.coalescent_correction = coalescent_correction

        self.taxa = sorted(t.taxon_namespace)
        self.chars = dict(zip(self.taxa, chars[:max_replicates, :].T))
        self.nchars = len(chars)
        self.max_replicates = max_replicates
        self.calculate_pic = calculate_pic

    def make_pic(self, tree=None):
        if tree is None:
            tree = self.t
        # print(len(self.taxa), self.chars)
        chars = {k: list(v) for k, v in self.chars.items()}
        chars = dendropy.ContinuousCharacterMatrix.from_dict(chars)
        self._pic = self.calculate_pic(tree=tree, char_matrix=chars)

    def estimate_sigmas(self, tree=None):
        """use full tree. If self.eps is provided, will delete paths with pendant edges < eps."""
        n_replicates = self.max_replicates
        if self.eps:
            tree = subsample_tree_random(self.t, self.eps)
        if tree is not None or not hasattr(self, "pic"):
            self.make_pic(tree)

        sigmas = np.empty(n_replicates)
        for i in range(n_replicates):
            ctree = self._pic.contrasts_tree(
                character_index=i,
                annotate_pic_statistics=True,
                state_values_as_node_labels=False,
                corrected_edge_lengths=False,
            )
            sigmas[i] = sigma_pic(ctree)
        return sigmas

    # def js_sigma(self):
    #     """deprecated. Pass JS calculator in constructor."""
    #     dists = self.t.calc_node_root_distances(
    #         return_leaf_distances_only=True
    #     )  # use to weight leaves if tree isnt ultrametric
    #     root_state_est = np.average(
    #         [self.chars[nd.taxon][character_index] for nd in self.t.leaf_node_iter()],
    #         weights=dists,
    #     )

    #     n_leaves = len(self.t)
    #     max_level = next(self.t.leaf_node_iter()).level()
    #     n_cherries = n_leaves - 1
    #     sig_sq = count = 0.0
    #     states = np.empty(n_cherries)
    #     h = np.empty(n_cherries)
    #     i = 0
    #     # last level only
    #     state_ests = {}
    #     c_sq = []
    #     for nd in self._tsqree.postorder_node_iter():
    #         child_nodes = nd.child_nodes()
    #         # if len(child_nodes) == 2 and all(child.is_leaf() for child in child_nodes):
    #         if nd.level() == level:
    #             c_sq.append(self.contrast_sq(nd))
    #             c1, c2 = child_nodes
    #             state_vals = np.array(
    #                 (
    #                     self._char_matrix[c1.taxon][character_index],
    #                     self._char_matrix[c2.taxon][character_index],
    #                 )
    #             )
    #             edges = np.array((c1.edge.length, c2.edge.length))
    #             num = (state_vals[1] - state_vals[0]) ** 2
    #             sum_of_child_edges = np.sum(edges)
    #             prod_of_child_edges = np.prod(edges)
    #             sig_sq += num / sum_of_child_edges

    #             states[i] = np.sum(state_vals / edges) / np.sum(1.0 / edges)
    #             h[i] = 2 * prod_of_child_edges / sum_of_child_edges
    #             i += 1
    #             self.t.prune(child_nodes)
    #         # elif nd.level() == level-1: break
    #     sig_sq /= count
    #     y = states / h
    #     denom = y - root_state_est
    #     num = sig_sq * (n_children - 3)
    #     states = h * ((1 - num / np.linalg.norm(denom)) * denom + y)

    def cherries_sigma(self, normalize=True):
        """use cherries only."""
        n_replicates = self.max_replicates
        t = self.t.clone()

        if self.eps:
            contrasts_sq = self.pruned_cherries(t)
        else:
            contrasts_sq = []
            for i in get_cherries(t):
                contrasts_sq.append(self.contrast_sq(i))

        return np.array(contrasts_sq).mean(axis=0)

    def contrast_sq(self, node):
        n1, n2 = node.child_nodes()
        tax1, tax2 = n1.taxon, n2.taxon
        v = n1.edge.length + n2.edge.length + 2 * self.coalescent_correction
        c = self.chars[tax1] - self.chars[tax2]
        return c**2 / v

    def pruned_cherries(self, t):
        t.calc_node_root_distances()
        cherries = list(get_cherries(t))
        taxa = {c: [n.taxon for n in c.leaf_nodes()] for c in cherries}
        found = True
        removed = 0
        while found:
            found = False
            for c1, c2 in combinations(cherries, 2):
                mrca = t.mrca(
                    taxa=(*taxa[c1], *taxa[c2])
                )  # [n.taxon for n in (*c1.leaf_nodes(), *c2.leaf_nodes())])
                if (c1.root_distance - mrca.root_distance < self.eps) or (
                    c2.root_distance - mrca.root_distance < self.eps
                ):
                    cherries.remove(c1)
                    found = True
                    removed += 1
                    # print(c1.root_distance - mrca.root_distance, c2.root_distance - mrca.root_distance)
                    break
        # print('removed', removed)

        return [self.contrast_sq(i) for i in cherries]

    def paired_lineages_sigma(self, normalize=True, max_replicates: int = None):
        """paired lineages.
        if eps is provided, will delete paths with pendant edges < eps.
        This is extreme - all we really need is > eps separation between any 2 paths."""
        n_replicates = self.nchars
        if max_replicates:
            n_replicates = min(n_replicates, max_replicates)

        t = self.t.clone()

        contrasts_sq = []
        filter_fn = lambda n: all(child.is_leaf() for child in n.child_node_iter())
        while len(t) > 1:
            for i in t.postorder_internal_node_iter(
                filter_fn=filter_fn, exclude_seed_node=False
            ):
                n1, n2 = i.child_nodes()

                if n1.is_leaf() and n2.is_leaf():
                    tax1, tax2 = n1.taxon, n2.taxon
                    if (not self.eps) or (i.edge.length > self.eps):
                        v = (
                            n1.edge.length
                            + n2.edge.length
                            + 2 * self.coalescent_correction
                        )
                        c = self.chars[tax1] - self.chars[tax2]
                        contrasts_sq.append(c**2 / v)
                    try:
                        t.prune_taxa([tax1, tax2])
                    except AttributeError:
                        pass
                    break
        #                 TODO: we miss the last contrast - attribute error bc dendropy wont allow zero length trees
        # print(len(contrasts_sq))
        return np.array(contrasts_sq).mean(axis=0)


if __name__ == "__main__":
    
    max_rep = 20
    ngenes = 1

    directories = np.random.permutation(
        [d for h in (7, 5, 3, 1)
            for d in parent_dir.glob(f"*/*/h{h}")]
    )
    
    #########################################################

    def process_dir_object(
        d: Path,
        overwrite: bool = True,
        ils=False,
        eps=None,
        ngenes=1,
        calculate_pic=continuous.PhylogeneticIndependentContrasts,
        max_replicates=None,
    ):
        tree_file = d / "parent_trees_4plus_pos_edges.phy"
        if ils:
            samp_file = d / f"samples_ngenes{ngenes}.npz"
            outdir = d
        else:
            samp_file = d / "no_ils/samples_ngenes1.npz"
            if ngenes != 1:
                raise ValueError("ngenes must be 1 unless ILS is True")
            outdir = d / "no_ils"
        outfile = Path(
            f"pic_js_sigmas{max_rep}_ngenes{ngenes}_full"
            + (f"_{eps}" if eps else "")
            + ".npy"
        )

        if (
            (not overwrite and outfile.exists())
            or not samp_file.exists()
            or not tree_file.exists()
        ):
            return 1

        with open(tree_file, "r") as f:
            parent_trees = [
                dendropy.Tree.get_from_string(t, "newick", rooting="force-rooted")
                for t in f
            ]
        npzfile = np.load(samp_file)
        n = len(npzfile.files)
        samples = [v / ngenes for v in npzfile.values()]
        sigmas = np.empty((3, n, max_rep))
        for i, (t, s) in enumerate(zip(parent_trees, samples)):
            try:
                pic = PIC(
                    t=t,
                    chars=s,
                    calculate_pic=calculate_pic,
                    max_replicates=max_replicates,
                    eps=eps,
                )
                res = np.vstack(
                    (
                        pic.estimate_sigmas(),
                        pic.cherries_sigma(),
                        pic.paired_lineages_sigma(),
                    )
                )
                sigmas[:, i, :] = res
            except ZeroDivisionError:
                print("error", i, samp_file)
        for suffix, sig in zip(("full", "cherries", "paired"), sigmas):
            outfile = (
                f"pic_js_sigmas{max_rep}_ngenes{ngenes}_{suffix}"
                + (f"_{eps}" if eps else "")
                + ".npy"
            )
            np.save(outdir / outfile, sig)
        return 0


    def process(dirname):
        return process_dir_object(
            dirname,
            overwrite=False,
            ngenes=ngenes,
            max_replicates=max_rep,
            ils=False,
            calculate_pic=JamesSteinPIC,
        )

    with Parallel(n_jobs=njobs) as parallel:
        skipped = parallel(
            delayed(process)(d) for d in directories
        )
    print("JS skipped:", sum(skipped))

    ############################################

    def process_dir_object(
        d: Path,
        overwrite: bool = False,
        ils=False,
        eps=None,
        ngenes=1,
        calculate_pic=continuous.PhylogeneticIndependentContrasts,
        max_replicates=None,
    ):
        tree_file = d / "parent_trees_4plus_pos_edges.phy"
        if ils:
            samp_file = d / f"samples_ngenes{ngenes}.npz"
            outdir = d
        else:
            samp_file = d / "no_ils/samples_ngenes1.npz"
            if ngenes != 1:
                raise ValueError("ngenes must be 1 unless ILS is True")
            outdir = d / "no_ils"
        outfile = outdir / (
            f"pic_bbridge_sigmas{max_rep}_ngenes{ngenes}_full{'_{eps}' if eps else ''}.npy"
        )
        # print(outfile)

        if (
            (not overwrite and outfile.exists())
            or not samp_file.exists()
            or not tree_file.exists()
        ):
            return 1
        # print(samp_file)

        with open(tree_file, "r") as f:
            parent_trees = [
                dendropy.Tree.get_from_string(t, "newick", rooting="force-rooted")
                for t in f
            ]
        npzfile = np.load(samp_file)
        n = len(npzfile.files)
        samples = [v / ngenes for v in npzfile.values()]
        sigmas = np.empty((3, n, max_rep))
        for i, (t, s) in enumerate(zip(parent_trees, samples)):
            try:
                pic = PIC(
                    t=t,
                    chars=s,
                    calculate_pic=calculate_pic,
                    max_replicates=max_replicates,
                    eps=eps,
                )
                res = np.vstack(
                    (
                        pic.estimate_sigmas(),
                        pic.cherries_sigma(),
                        pic.paired_lineages_sigma(),
                    )
                )
                sigmas[:, i, :] = res
            except ZeroDivisionError:
                print(f"error on tree {i} in file {treefile}")
        for suffix, sig in zip(("full", "cherries", "paired"), sigmas):
            outfile = f"pic_bbridge_sigmas{max_rep}_ngenes{ngenes}_{suffix}{'_{eps}' if eps else ''}.npy"

            np.save(outdir / outfile, sig)
        return 0

    print(njobs, os.cpu_count(), parent_dir)

    def process(dirname):
        return process_dir_object(
            dirname,
            overwrite=False,
            ngenes=ngenes,
            max_replicates=max_rep,
            ils=False,
            calculate_pic=BrownianBridgePIC,
        )

    with Parallel(n_jobs=njobs) as parallel:
        skipped = parallel(
            delayed(process)(d) for d in np.random.permutation(directories)            
            # delayed(process)(d) for d in parent_dir.glob("tree_sims/*/*/h*")
        )
    print("skipped:", sum(skipped))
