"""
Clustering utilities for scikit-bio TreeNode objects.
"""

from __future__ import annotations

import numpy as np
from numba import njit
from skbio import TreeNode


# --- Helper: build distance matrix from skbio tree ---
def tree_to_distance_matrix(tree: TreeNode):
    """
    Convert a skbio TreeNode into a symmetric distance matrix of tip-to-tip patristic distances.

    Parameters
    ----------
    tree : TreeNode
        A rooted phylogenetic tree with branch lengths.

    Returns
    -------
    tuple[list[str], numpy.ndarray]
        Tip names in order and the corresponding symmetric distance matrix.
    """
    tips = list(tree.tips())
    names = [t.name for t in tips]
    n = len(tips)
    D = np.zeros((n, n), dtype=np.float64)
    for i, ti in enumerate(tips):
        for j in range(i + 1, n):
            d = ti.distance(tips[j])  # patristic distance
            D[i, j] = d
            D[j, i] = d
    return names, D


# --- Numba union–find clustering ---
@njit
def find(parent, x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]  # path compression
        x = parent[x]
    return x


@njit
def union(parent, rank, a, b):
    ra = find(parent, a)
    rb = find(parent, b)
    if ra == rb:
        return
    if rank[ra] < rank[rb]:
        parent[ra] = rb
    elif rank[ra] > rank[rb]:
        parent[rb] = ra
    else:
        parent[rb] = ra
        rank[ra] += 1


@njit()
def cluster(D, threshold):
    """
    Cluster indices of a distance matrix using union–find under a threshold.

    Parameters
    ----------
    D : ndarray (n x n)
        Symmetric distance matrix.
    threshold : float
        Maximum distance cutoff for clustering.

    Returns
    -------
    numpy.ndarray
        Cluster labels (starting at 1).
    """
    n = D.shape[0]
    parent = np.arange(n)
    rank = np.zeros(n, dtype=np.int64)

    # Union all pairs with distance <= threshold
    for i in range(n):
        for j in range(i + 1, n):
            if D[i, j] <= threshold:
                union(parent, rank, i, j)

    # Find canonical representatives
    reps = np.empty(n, dtype=np.int64)
    for i in range(n):
        reps[i] = find(parent, i)

    # Map representatives to consecutive cluster IDs
    rep_to_id = -np.ones(n, dtype=np.int64)
    cur_id = 1
    labels = np.empty(n, dtype=np.int64)
    for i in range(n):
        r = reps[i]
        if rep_to_id[r] == -1:
            rep_to_id[r] = cur_id
            cur_id += 1
        labels[i] = rep_to_id[r]

    return labels


# --- Wrapper: cluster tips directly from tree ---
def cluster_tips(tree: TreeNode, threshold: float):
    """
    Cluster tips of a skbio TreeNode by patristic distance threshold.

    Parameters
    ----------
    tree : TreeNode
        A rooted phylogenetic tree with branch lengths.
    threshold : float
        Maximum distance cutoff for clustering.

    Returns
    -------
    dict
        Mapping cluster_id -> list of tip names.
    """
    tip_names, D = tree_to_distance_matrix(tree)
    labels = cluster(D, float(threshold))
    cluster_dict = {}
    for name, cid in zip(tip_names, labels):
        cluster_dict.setdefault(int(cid), []).append(name)
    return cluster_dict


# --- Find lowest common taxa for clusters ---
def lowest_common_taxa(cluster_dict, tax_table):
    """Find the lowest common taxa name for each cluster."""
    lca_dict = {}
    lca_levels = []
    for cluster_id, tips in cluster_dict.items():
        taxa_list = [tax_table.loc[tip].values for tip in tips]
        taxa_transposed = list(zip(*taxa_list))
        lca = None
        level_names = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
        level_index = 0
        for taxa_level in taxa_transposed:
            if len(set(taxa_level)) == 1:
                lca = taxa_level[0]
            else:
                break
            level_index += 1
        level_index = min(level_index, len(level_names) - 1)
        lca_levels.append(level_names[level_index])
        lca_dict[cluster_id] = lca if lca is not None else "No common taxa"

    return lca_dict, np.array(lca_levels)


# --- Warm-up: trigger JIT compilation at import ---
def _warmup():
    D = np.zeros((2, 2), dtype=np.float64)
    labels = cluster(D, 0.0)
    return labels


_warmup()


if __name__ == "__main__":  # pragma: no cover - usage example
    newick = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"
    tree = TreeNode.read([newick])
    threshold = 0.5
    clusters = cluster_tips(tree, threshold)
    for cid, tips in clusters.items():
        print(f"Cluster {cid}: {tips}")
