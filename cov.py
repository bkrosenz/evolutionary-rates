import sympy
import numpy as np
from typing import List, Tuple, Optional

def _sympify_branch_length(length_str: str):
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
    """A simple class to represent a node in a phylogenetic tree."""
    def __init__(self, name: Optional[str] = None, parent: Optional['Node'] = None):
        self.name = name
        self.branch_length = sympy.sympify(0)
        self.parent = parent
        self.children: List['Node'] = []

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
        # Tokenize the string for easier parsing
        tokens = newick_string.strip().replace('(', ' ( ').replace(')', ' ) ') \
            .replace(',', ' , ').replace(':', ' : ').replace(';', ' ; ').split()

        root = Node(name="root_sentinel")
        current_node = root
        
        for i, token in enumerate(tokens):
            if token == '(':
                new_node = Node(parent=current_node)
                current_node.children.append(new_node)
                current_node = new_node
            elif token == ',':
                if current_node.parent:
                    current_node = current_node.parent
                    new_node = Node(parent=current_node)
                    current_node.children.append(new_node)
                    current_node = new_node
            elif token == ')':
                if current_node.parent:
                    current_node = current_node.parent
            elif token == ':':
                # The next token is the branch length. It will be handled below.
                pass
            elif token == ';':
                break
            else: # It's a name or a branch length
                if i > 0 and tokens[i-1] == ':':
                    current_node.branch_length = _sympify_branch_length(token)
                else:
                    current_node.name = token

        # The actual root is the single child of our sentinel root
        return root
        # if len(root.children) == 1:
        #      final_root = root.children[0]
        #      final_root.parent = None
        #      return final_root
        # else:
        #     # This can happen for malformed strings, e.g., "(A,B)(C,D);"
        #     print("Error: Malformed Newick string detected.")
        #     return None
    except Exception as e:
        print(f"An error occurred during parsing: {e}")
        return None

def get_leaves(node: Node) -> List[Node]:
    """Recursively finds all leaf nodes (tips) in the subtree of a given node."""
    if not node.children:
        return [node]
    leaves = []
    for child in node.children:
        leaves.extend(get_leaves(child))
    return leaves

def newick_to_phylo_cov(newick_string: str) -> Tuple[Optional[np.ndarray], List[str]]:
    """
    Converts a Newick-formatted phylogenetic tree to a phylogenetic covariance matrix.

    This version uses a custom parser and NumPy instead of DendroPy and Pandas.

    The covariance between two taxa is the sum of branch lengths from the root of
    the tree to their most recent common ancestor (MRCA). The variance of a

    taxon (the diagonal elements) is the sum of branch lengths from the root to
    that taxon's tip.

    Args:
        newick_string: A string containing the tree in Newick format. Branch
                       lengths can be numeric or symbolic variables (e.g., "a", "b").

    Returns:
        A tuple containing:
        - A NumPy array representing the phylogenetic covariance matrix.
        - A sorted list of taxa labels corresponding to the matrix rows/columns.
        Returns (None, []) on failure.
    """
    root = parse_newick(newick_string)
    if not root:
        return None, []

    # 1. Get a sorted list of leaf nodes (taxa) for consistent matrix ordering
    leaf_nodes = sorted(get_leaves(root), key=lambda n: n.name or "")
    taxa_labels = [node.name for node in leaf_nodes]
    n_taxa = len(taxa_labels)

    # 2. Pre-calculate paths from each leaf to the root for efficiency
    paths_to_root = {}
    for leaf in leaf_nodes:
        path = []
        curr = leaf
        while curr:
            path.append(curr)
            curr = curr.parent
        paths_to_root[leaf.name] = path
    
    # 3. Initialize the covariance matrix with sympy objects
    cov_matrix = np.empty((n_taxa, n_taxa), dtype=object)

    # 4. Calculate pairwise covariances
    for i in range(n_taxa):
        for j in range(i, n_taxa):
            label1 = taxa_labels[i]
            label2 = taxa_labels[j]
            
            # Find the MRCA by finding the first common node in their paths to the root
            path1 = paths_to_root[label1]
            path2_nodes = set(paths_to_root[label2])
            mrca_node = None
            for node in path1:
                if node in path2_nodes:
                    mrca_node = node
                    break
            
            if mrca_node is None:
                # This should not happen in a valid tree structure
                return None, []

            # The covariance is the path length from the root to the MRCA
            path_length_to_root = sympy.sympify(0)
            curr_node = mrca_node
            while curr_node:
                path_length_to_root += curr_node.branch_length
                curr_node = curr_node.parent
            
            # The matrix is symmetric
            cov_matrix[i, j] = path_length_to_root
            cov_matrix[j, i] = path_length_to_root
            
    return cov_matrix, taxa_labels

if __name__ == '__main__':
    # --- Example Usage ---
    newick_tree = "(((Human:a, Chimp:0.01):b, Gorilla:c):d, Orangutan:e);"

    print("--- Phylogenetic Covariance Matrix Generator (NumPy/SymPy version) ---")
    print(f"\nInput Newick Tree:\n{newick_tree}\n")
    
    covariance_matrix, labels = newick_to_phylo_cov(newick_tree)

    if covariance_matrix is not None:
        print("Taxa Labels:", labels)
        print("\nResulting Covariance Matrix (NumPy array of SymPy expressions):")
        print(covariance_matrix)

    # --- Another Example: A star tree with a root branch length ---
    star_tree = "(T1:x, T2:y, T3:z):L;"
    print("\n" + "="*50)
    print(f"\nInput Newick Tree (Star Phylogeny):\n{star_tree}\n")
    star_cov_matrix, star_labels = newick_to_phylo_cov(star_tree)
    if star_cov_matrix is not None:
        print("Taxa Labels:", star_labels)
        print("\nResulting Covariance Matrix:")
        print(star_cov_matrix)

