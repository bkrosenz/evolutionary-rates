import sympy
import math
from typing import List, Tuple, Optional, Dict, Union
from cov import *
from calculate_pic_var import *

# --- Data Structures and Parser (from previous version) ---


def _sympify_branch_length(length_str: str) -> sympy.Expr:
    """
    Converts a branch length string to a sympy object.

    If the length can be interpreted as a number, it becomes a sympy Number.
    Otherwise, it becomes a sympy Symbol. None or empty string is treated as zero.

    Args:
        length_str: The branch length value as a string (e.g., '0.1', 'a').

    Returns:
        A sympy object representing the branch length.
    """
    if not length_str:
        return sympy.sympify(0)

    try:
        # Try to convert to a number first
        return sympy.sympify(float(length_str))
    except (ValueError, TypeError):
        # If it fails, treat it as a symbol
        return sympy.Symbol(length_str)


class Node:
    """
    A class to represent a node in a phylogenetic tree.
    It includes attributes needed for PIC calculation.
    """

    def __init__(self, name: Optional[str] = None, parent: Optional["Node"] = None):
        self.name = name
        self.branch_length: sympy.Expr = sympy.sympify(0)
        self.parent = parent
        self.children: List["Node"] = []
        # Attribute to hold trait data for PIC
        self.trait_value: Optional[float] = None

    def __repr__(self) -> str:
        return f"Node(name='{self.name}', length={self.branch_length})"


def parse_newick(newick_string: str) -> Optional[Node]:
    """
    Parses a Newick string into a tree of Node objects without external libraries.

    Args:
        newick_string: The tree in Newick format.

    Returns:
        The root Node of the parsed tree, or None if parsing fails.
    """
    try:
        tokens = (
            newick_string.strip()
            .replace("(", " ( ")
            .replace(")", " ) ")
            .replace(",", " , ")
            .replace(":", " : ")
            .replace(";", " ; ")
            .split()
        )

        if not tokens or tokens[-1] != ";":
            print("Error: Newick string must end with a semicolon ';'.")
            return None

        root = Node(name="root_sentinel")
        current_node = root

        for i, token in enumerate(tokens):
            if token == "(":
                new_node = Node(parent=current_node)
                current_node.children.append(new_node)
                current_node = new_node
            elif token == ",":
                if current_node.parent:
                    current_node = current_node.parent
                    new_node = Node(parent=current_node)
                    current_node.children.append(new_node)
                    current_node = new_node
            elif token == ")":
                if current_node.parent:
                    current_node = current_node.parent
            elif token == ":":
                pass  # The next token is the branch length.
            elif token == ";":
                break
            else:  # It's a name or a branch length
                if i > 0 and tokens[i - 1] == ":":
                    current_node.branch_length = _sympify_branch_length(token)
                else:
                    current_node.name = token
        return root
    except Exception as e:
        print(f"An error occurred during parsing: {e}")
        return None


# --- Felsenstein's Independent Contrasts Method ---


def calculate_pic_rate(
    newick_string: str, tip_data: Dict[str, Union[int, float, sympy.Symbol]]
) -> Tuple[float, List[float]]:
    """
    Calculates the mutation rate (sigma^2) using Felsenstein's PIC method.

    This function implements a recursive post-order traversal of the tree to
    compute n-1 independent contrasts, where n is the number of taxa. The
    maximum likelihood estimate of the evolutionary rate under a Brownian
    motion model is the mean of the squares of these contrasts.

    Args:
        newick_string: The phylogenetic tree in Newick format with numeric
                       branch lengths.
        tip_data: A dictionary mapping taxon names (strings) at the tips of
                  the tree to their continuous trait values (floats or ints).

    Returns:
        A tuple containing:
        - The estimated mutation rate (sigma^2) as a float.
        - A list of the calculated independent contrasts.

    Raises:
        ValueError: If the tree is not bifurcating, if tip data is missing,
                    or if the number of tips does not match the data.
        TypeError: If branch lengths or trait values are not numeric.
    """
    root = parse_newick(newick_string)
    if not root:
        raise ValueError("Failed to parse the Newick tree string.")

    # Attach tip data to leaf nodes
    leaves = (
        [node for node in root.children if not node.children]
        if not root.children
        else get_leaves(root)
    )

    if len(leaves) != len(tip_data):
        raise ValueError(
            "The number of tips in the tree does not match the number of entries in tip_data."
        )

    for leaf in leaves:
        if leaf.name in tip_data:
            leaf.trait_value = tip_data[leaf.name]
            if isinstance(leaf.trait_value, str):
                leaf.trait_value = float(leaf.trait_value)
        else:
            raise ValueError(f"No trait data found for tip: {leaf.name}")

    contrasts: List[float] = []

    def _pic_recursion(node: Node) -> Tuple[float, float]:
        """
        Inner recursive function to compute contrasts.
        Returns a tuple of (estimated_ancestral_state, variance_at_node).
        """
        if not node.children:  # It's a leaf
            if node.trait_value is None:
                # This should be caught earlier, but as a safeguard:
                raise ValueError(f"Leaf node {node.name} is missing trait data.")
            return node.trait_value, 0.0

        if len(node.children) != 2:
            raise ValueError(
                "Phylogenetic Independent Contrasts requires a fully bifurcating (resolved) tree."
            )

        # Recurse on children
        child1, child2 = node.children
        x1, v1 = _pic_recursion(child1)
        x2, v2 = _pic_recursion(child2)

        # Check for symbolic branch lengths
        # if not child1.branch_length.is_Number or not child2.branch_length.is_Number:
        #     raise TypeError("Branch lengths must be numeric for PIC calculation.")

        # Calculate total variance from children (node variance + branch length)
        V1 = v1 + (child1.branch_length)
        V2 = v2 + (child2.branch_length)

        # Calculate and store the standardized contrast
        contrast = (x1 - x2) / sympy.sqrt(V1 + V2)
        contrasts.append(contrast)

        # Estimate the ancestral state at the current node
        ancestral_state = (x1 / V1 + x2 / V2) / (1 / V1 + 1 / V2)

        # Calculate the expected variance of the estimate at this node
        node_variance = (V1 * V2) / (V1 + V2)

        return ancestral_state, node_variance

    # Start the recursion from the root
    _pic_recursion(root)

    # The rate sigma^2 is the sum of squared contrasts divided by n-1
    # which is equivalent to the number of contrasts.
    if not contrasts:
        return 0.0, []

    sigma_squared = sum(c**2 for c in contrasts) / len(contrasts)

    return sigma_squared, contrasts


# --- NEW: Cherry-Based Contrasts Method ---


def find_cherries(node: Node) -> List[Node]:
    """
    Recursively finds all 'cherry' nodes in the tree.
    A cherry is an internal node whose children are both leaves (tips).
    """
    found_cherries = []

    # A leaf node cannot contain any cherries.
    if not node.children:
        return []

    # Check if the current node is a cherry.
    is_a_cherry = (
        len(node.children) == 2
        and not node.children[0].children
        and not node.children[1].children
    )
    if is_a_cherry:
        found_cherries.append(node)

    # Recursively search for cherries in the children's subtrees.
    for child in node.children:
        found_cherries.extend(find_cherries(child))

    return found_cherries


def calculate_cherry_pic_rate(
    newick_string: str, tip_data: Dict[str, Union[int, float]]
) -> Tuple[float, List[float]]:
    """
    Calculates mutation rate (sigma^2) using contrasts from only the
    'cherries' (sister taxa pairs) of the tree.

    Args:
        newick_string: The phylogenetic tree in Newick format.
        tip_data: A dictionary mapping taxon names to their trait values.

    Returns:
        A tuple containing:
        - The estimated mutation rate (sigma^2) as a float.
        - A list of the calculated cherry-based contrasts.

    Raises:
        ValueError: If parsing fails or tip data is missing for a cherry.
        TypeError: If branch lengths for a cherry are not numeric.
    """
    root = parse_newick(newick_string)
    if not root:
        raise ValueError("Failed to parse the Newick tree string.")

    # Attach tip data to all leaf nodes first
    all_leaves = get_leaves(root)
    for leaf in all_leaves:
        if leaf.name in tip_data:
            leaf.trait_value = float(tip_data[leaf.name])
        else:
            # We don't raise an error here yet, only if data for a cherry is missing
            pass

    # Find all cherries in the tree
    cherries = find_cherries(root)
    if not cherries:
        print("Warning: No cherries found in the tree.")
        return 0.0, []

    contrasts = []
    for cherry_node in cherries:
        child1, child2 = cherry_node.children

        x1 = child1.trait_value
        x2 = child2.trait_value

        if x1 is None or x2 is None:
            raise ValueError(
                f"Missing trait data for cherry with taxa '{child1.name}' and '{child2.name}'."
            )

        if not child1.branch_length.is_Number or not child2.branch_length.is_Number:
            raise TypeError(
                f"Branch lengths for cherry '{child1.name}, {child2.name}' must be numeric."
            )

        len1 = child1.branch_length
        len2 = child2.branch_length

        # Calculate the standardized contrast for this cherry
        # The variance for a tip is just its branch length.
        contrast = (x1 - x2) / sympy.sqrt(len1 + len2)
        contrasts.append(contrast)

    if not contrasts:
        return 0.0, []

    # Rate is the mean of the squared contrasts
    sigma_squared = sum(c**2 for c in contrasts) / len(contrasts)

    return sigma_squared, contrasts


def get_leaves(node: Node) -> List[Node]:
    """Helper function to recursively find all leaf nodes."""
    if not node.children:
        return [node]
    leaves = []
    for child in node.children:
        leaves.extend(get_leaves(child))
    return leaves


if __name__ == "__main__":
    # --- Example Usage ---
    # This example uses the tree and data from a well-known PIC tutorial
    # to verify the implementation.

    # Tree: ((A:0.1, B:0.1):0.2, C:0.3);
    # Trait values: A=1, B=2, C=4

    example_tree = "((A:0.1, B:0.1):0.2, C:0.3);"
    example_data = dict(zip("ABC", (sympy.Symbol(f"x{i}") for i in range(3))))
    example_data_values = {"A": 1.0, "B": 2.0, "C": 4.0}

    print("--- Felsenstein's Independent Contrasts Calculator ---")
    print(f"\nInput Newick Tree:\n{example_tree}")
    print(f"\nInput Tip Data:\n{example_data}\n")

    # Calculate the rate
    rate, calculated_contrasts = calculate_pic_rate(example_tree, example_data)

    print(f"Calculated Contrasts: {calculated_contrasts}")
    print(f"Number of Contrasts (n-1): {len(calculated_contrasts)}")
    print(f"\nEstimated Evolutionary Rate (sigma^2): {rate}\n")
    print(
        f"\nEstimated Evolutionary Rate substituted (sigma^2): {rate.subs({example_data[k]:example_data_values[k] for k in example_data})}\n"
    )

    # --- Felsenstein (1985) Example ---
    # A more complex example based on data from the original paper.
    felsenstein_tree = "(((A:1, B:1):1, C:2):1, D:3);"
    felsenstein_data = {"A": -0.447, "B": -0.566, "C": 0.019, "D": 0.948}

    print("=" * 50)
    print("\nExample from Felsenstein (1985):")
    print(f"Tree: {felsenstein_tree}")
    print(f"Data: {felsenstein_data}\n")

    rate_f, contrasts_f = calculate_pic_rate(felsenstein_tree, felsenstein_data)
    print(f"Calculated Contrasts: {contrasts_f}")
    print(f"Estimated Evolutionary Rate (sigma^2): {rate_f:.4f}")

    n_jobs = 16 

    with open("tree_topologies.nw", "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            s = line.strip()
            print(f"\nProcessing tree: {s}")

            cov, labs = newick_to_phylo_cov(s)
            cov = cov * sigma**2
            traits = {k: sympy.Symbol(f"x{i+1}") for i, k in enumerate(labs)}
            rate_f, contrasts_f = calculate_pic_rate(s, traits)
            moments = list(
                reversed(generate_moments(list(traits.values()), cov, cov.shape[1] + 1))
            )

            def substitute_moments(expr):
                for m in moments:
                    expr = expr.subs(m, simultaneous=True)
                return expr

            pic_e = sympy.expand(rate_f)
            func = pic_e.func
            with Parallel(n_jobs=n_jobs) as parallel:
                args = parallel(delayed(substitute_moments)(arg) for arg in pic_e.args)
            pic_e = simplify(func(*args))
            print(r"Expected value of PIC estimator ($\hat{\sigma}^2$):", f"{pic_e}")

            pic_v = sympy.expand(rate_f**2 - pic_e**2)
            func = pic_v.func
            with Parallel(n_jobs=n_jobs) as parallel:
                args = parallel(delayed(substitute_moments)(arg) for arg in pic_v.args)
            pic_v = simplify(func(*args))
            print(r"Variance of PIC estimator ($E[\hat{\sigma}^2-\sigma^2]^2$):", f"{pic_v}\n----\n")
