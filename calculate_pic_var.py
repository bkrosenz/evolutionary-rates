# Import the necessary modules from SymPy
from sympy import symbols, simplify, factor, expand, lambdify, Matrix, prod
import numpy as np
import math
from joblib import Parallel, delayed
import itertools
import pickle, gzip

from itertools import combinations, permutations, combinations_with_replacement
from operator import mul
from functools import *

n_jobs=10

# Define the symbols used in the expression
# 'a', 'b', 'c', 'd' are general variables.
# 'e', 'eps' (epsilon), 'delta' are also variables, often representing small quantities or parameters.
(
    a,
    b,
    c,
    d,
    e,
    f,
    t0,
    t1,
    t2,
    t3,
    t4,
    sigma,
) = symbols("a b c d e f t0 t1 t2 t3 t4 sigma")

variables = [
    a,
    b,
    c,
    d,
    e,
    f,
]
branching_times = [t0, t1, t2, t3, t4]
root_to_tip = sum(branching_times)

# Define the given mathematical expression
# The expression is broken down into three main terms for clarity,
# then summed together.

contrast = lambda a, b, e: (a - b) ** 2 / e


def f_full3(a, b, c):
    """Tree: ((a:t0, b:t0):t1, c:(t0+t1));"""
    c1 = (a - b) ** 2 / (2 * t0)
    c2 = (c - (a + b) / 2) ** 2 / (2 * t1 + 3 * t0 / 2)
    return (c1 + c2) / 2


def trifurcating_full3(a, b, c):
    """Tree: ((a,b,c);"""
    return (a - b) ** 2 / (2 * t0) + (c - (a + b) / 2) ** 2 / (3 * t0 / 2)


def f_full(a, b, c, d, e, f):
    c1 = contrast(a, b, 2 * t0)
    c2 = contrast(c, d, 2 * (t0 + t1))
    c3 = contrast(e, f, 2 * (t0 + t1 + t2))
    x1 = (a + b) / 2
    x2 = (c + d) / 2
    x3 = (e + f) / 2
    e1 = t0 / 2 + t1 + t2 + t3
    e2 = (t0 + t1) / 2 + t2 + t3
    e3 = (t0 + t1 + t2) / 2 + t3 + t4
    denom = e1 + e2
    c4 = contrast(x1, x2, denom)
    x4 = (e1 * x2 + e2 * x1) / (e1 + e2)
    e4 = (e1 * e2) / (e1 + e2)
    c5 = contrast(x3, x4, e3 + e4 + t4)
    return (c1 + c2 + c3 + c4 + c5) / 5


def f_full_caterpillar(a, b, c, d, e, f):
    """Tree shape: (((a,b),c),d);"""
    #TODO: finish this
    c1 = contrast(a, b, 2 * t0)
    x1 = (a + b) / 2
    e1 = t0 / 2 + sum(branching_times[1:])
    e2 = (t0 + t1) / 2 + sum(branching_times[2:])
    e3 = (t0 + t1 + t2) / 2 + sum(branching_times[3:])
    denom = e1 + e2  # [:,np.newaxis]
    # print(x1.shape,c.shape,denom.shape, ep.shape)
    c4 = contrast(x1, c, denom)
    x4 = (e1 * c + e2 * x1) / (e1 + e2)
    e4 = (e1 * e2) / (e1 + e2)
    c5 = contrast(d, x4, e3 + e4 + 2 * t4 + t3)
    return (c1 + c2 + c3 + c4 + c5) / 5


def enumerate_k_partitions(n, k):
    """
    Enumerates all unique ways to partition n distinct elements into sets of size k.

    This function works by recursively finding sets of k elements. It fixes the
    first available element and forms a set of k elements by combining it with
    k-1 other available elements, then recursively finds partitions for the
    remaining elements.

    Args:
        n (int): The total number of elements to partition. Elements are represented by
                 integers from 0 to n-1.
        k (int): The desired size of each set in the partition.

    Returns:
        list: A list of lists, where each inner list represents a unique partition.
              Each partition is a list of tuples, where each tuple is a set of k elements
              (e.g., (0, 1, 2)).
              Returns an empty list if k does not divide n, or if n, k are invalid.

    Raises:
        ValueError: If n is negative or k is not a positive integer.
    """
    if n < 0:
        raise ValueError("Number of elements (n) cannot be negative.")
    if k <= 0:
        raise ValueError("Set size (k) must be a positive integer.")
    if n % k != 0:
        # It's impossible to partition n elements into sets of size k if k does not divide n.
        # print(f"Warning: Cannot partition {n} elements into sets of size {k} because {k} does not divide {n}.")
        return []
    if n == 0:
        # Base case: There's one way to partition an empty set (an empty partition).
        return [[]]

    # Initialize the list of elements as integers from 0 to n-1
    elements = list(range(n))
    all_partitions = []

    # The core recursive function
    def _find_partitions(current_elements):
        # Base case: If no elements left, we found a valid partition.
        if not current_elements:
            return [[]]

        # Pick the first available element to start a new set
        first_element = current_elements[0]
        # Elements available for forming the rest of the current set
        remaining_to_choose_from = current_elements[1:]

        partitions_for_current_level = []

        # We need to choose k-1 additional elements to form a set of size k
        # from the remaining_to_choose_from.
        # This is where itertools.combinations becomes useful.
        for combo_of_k_minus_1 in itertools.combinations(
            remaining_to_choose_from, k - 1
        ):
            # Form the current set, ensuring elements are sorted for uniqueness
            current_set = tuple(sorted((first_element,) + combo_of_k_minus_1))

            # Create the list of elements remaining after forming this set
            # This involves removing 'first_element' and all elements in 'combo_of_k_minus_1'
            elements_in_current_set = set(current_set)  # Use a set for efficient lookup
            elements_after_set = [
                e for e in current_elements if e not in elements_in_current_set
            ]

            # Recursively find partitions for the remaining elements
            sub_partitions = _find_partitions(elements_after_set)

            # Combine the current set with all sub-partitions found
            for sub_p in sub_partitions:
                partitions_for_current_level.append([current_set] + sub_p)

        return partitions_for_current_level

    all_partitions = _find_partitions(elements)
    return all_partitions


def generate_moments(variables, cov, max_moment=None):
    n = len(variables)
    r = range(n)
    moments = []
    for k in range(1, max_moment or n + 1):
        m = {}
        for x in combinations_with_replacement(r, k):
            indices = enumerate_k_partitions(k, 2)
            m[reduce(mul, (variables[i] for i in x))] = sum(
                reduce(mul, (cov[x[i], x[j]] for i, j in ix)) for ix in indices
            )
        moments.append(m)
    return moments

if __name__ == "__main__":
    # balanced 6 taxa

    cov = (
        Matrix(
            [
                [root_to_tip, root_to_tip - t0, t4, t4, 0, 0],
                [root_to_tip - t0, root_to_tip, t4, t4, 0, 0],
                [t4, t4, root_to_tip, t2 + t3 + t4, 0, 0],
                [t4, t4, t2 + t3 + t4, root_to_tip, 0, 0],
                [0, 0, 0, 0, root_to_tip, t3 + t4],
                [0, 0, 0, 0, t3 + t4, root_to_tip],
            ]
        )
        * sigma**2
    )
    moments = list(reversed(generate_moments(variables, cov, cov.shape[1] + 1)))

    pic = expand(f_full(*variables))
    pic_v = expand(pic**2 - sigma**4)
    fn = "pic_v_6balanced_terms.txt.gz"
    #    with open(fn, "wb") as f:
    #        for term in pic_v.args:
    #            f.write(str(term) + "\n")
    #        f.write('---------'
    #        )
    #    with gzip.open('pic_balanced6_terms.pkl.gz', 'wb') as f:
    #        pickle.dump(args, f)

    print(f"wrote terms to {fn}, reduction is {pic_v.func}")

    def substitute_moments(expr):
        for m in moments:
            expr = expr.subs(m, simultaneous=True)
        return expr

    f = pic_v.func

    with Parallel(n_jobs = n_jobs) as parallel:
        args = parallel(delayed(substitute_moments)(arg) for arg in pic_v.args)

    with gzip.open("pic_v.pkl.gz", "wb") as outf:
        pickle.dump(args, outf)

    print('wrote pic var obj')

    # with open("pic_v.pkl") as inf:
    #     y = pickle.load(inf)

    pic_v = simplify(f(*args))
    print(pic_v)

    with gzip.open("pic_simple.pkl.gz", "wb") as outf:
        pickle.dump(pic_v, outf)
