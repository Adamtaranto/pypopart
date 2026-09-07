"""
Parsimony trees and ancestral-state sampling for the Parsimony Network.

Ports the relevant behaviour of PopART's tree/ParsimonyTree.cpp with
uniform substitution costs: a Fitch up-pass computes the parsimony
score and per-site state sets, and a randomised down-pass samples one
concrete ancestral sequence per internal node (PopART's
computeAncestors resolves ties randomly on every call, which is what
makes repeated edge sampling explore the ancestor space).

PopART receives its parsimony trees from a Nexus TREES block; PyPopART
has no tree input path yet, so trees are generated here by random-order
stepwise addition, keeping each taxon's best-scoring insertion point.
"""

import random
from typing import Dict, List, Sequence as SequenceType, Tuple

import networkx as nx

#: Wildcard state set used for gaps and ambiguity codes (PopART assigns
#: zero cost to every nucleotide for '-').
_ALL_STATES = frozenset('ACGT')


def _site_states(char: str) -> frozenset:
    """
    Map a sequence character to its Fitch state set.

    Parameters
    ----------
    char : str
        Upper-case sequence character.

    Returns
    -------
    frozenset
        Allowed nucleotide states at the site.
    """
    if char in ('U', 'u'):
        return frozenset('T')
    if char in 'ACGT':
        return frozenset(char)
    # Gaps and ambiguity codes cost nothing anywhere, as in SankoffUp
    return _ALL_STATES


class ParsimonyTree:
    """
    An unrooted tree over haplotype indices with Fitch machinery.

    Parameters
    ----------
    topology : nx.Graph
        Tree with leaf nodes 0..n-1 (haplotype indices) and internal
        nodes labelled with negative integers.
    sequences : list of str
        Condensed sequence per leaf index.
    weights : list of int
        Site weights for the condensed columns.
    """

    def __init__(
        self,
        topology: nx.Graph,
        sequences: SequenceType[str],
        weights: SequenceType[int],
    ):
        """
        Initialize the tree.

        Parameters
        ----------
        topology : nx.Graph
            Tree topology (leaves are haplotype indices >= 0).
        sequences : list of str
            Condensed leaf sequences.
        weights : list of int
            Site weights.
        """
        self.topology = topology
        self.sequences = list(sequences)
        self.weights = list(weights)
        self._n_sites = len(self.sequences[0]) if self.sequences else 0
        self._state_sets: Dict[int, List[frozenset]] = {}
        self._score: int = 0
        self._up_pass_done = False

    def compute_score(self) -> int:
        """
        Run the Fitch up-pass and return the weighted parsimony score.

        Returns
        -------
        int
            Weighted number of state changes implied by the tree.
        """
        root = self._pick_root()
        order = list(nx.dfs_postorder_nodes(self.topology, source=root))
        parent = nx.dfs_predecessors(self.topology, source=root)

        self._state_sets = {}
        score = 0
        for node in order:
            if node >= 0:  # leaf
                self._state_sets[node] = [_site_states(c) for c in self.sequences[node]]
                continue
            children = [
                neighbour
                for neighbour in self.topology.neighbors(node)
                if parent.get(neighbour) == node
            ]
            sets: List[frozenset] = []
            for site in range(self._n_sites):
                child_sets = [self._state_sets[c][site] for c in children]
                intersection = frozenset.intersection(*child_sets)
                if intersection:
                    sets.append(intersection)
                else:
                    # Pairwise Fitch union accumulation for >2 children
                    merged = child_sets[0]
                    for child_set in child_sets[1:]:
                        overlap = merged & child_set
                        if overlap:
                            merged = overlap
                        else:
                            merged = merged | child_set
                            score += self.weights[site]
                    sets.append(merged)
            self._state_sets[node] = sets

        self._root = root
        self._parent = parent
        self._score = score
        self._up_pass_done = True
        return score

    def sample_ancestors(self, rng: random.Random) -> Dict[int, str]:
        """
        Sample one concrete sequence per internal node (random Fitch).

        Ties are resolved randomly on every call, as in PopART's
        computeAncestors, so repeated calls explore alternative equally
        parsimonious ancestral reconstructions.

        Parameters
        ----------
        rng : random.Random
            Random source.

        Returns
        -------
        dict
            Mapping internal node id to its sampled sequence.
        """
        if not self._up_pass_done:
            self.compute_score()

        ancestors: Dict[int, str] = {}
        order = list(nx.dfs_preorder_nodes(self.topology, source=self._root))
        assigned: Dict[int, List[str]] = {}
        for node in order:
            if node >= 0:
                continue
            parent = self._parent.get(node)
            chars: List[str] = []
            for site in range(self._n_sites):
                states = self._state_sets[node][site]
                if parent is not None and parent < 0:
                    parent_char = assigned[parent][site]
                    if parent_char in states:
                        chars.append(parent_char)
                        continue
                chars.append(rng.choice(sorted(states)))
            assigned[node] = chars
            ancestors[node] = ''.join(chars)
        return ancestors

    def edge_sequences(self, ancestors: Dict[int, str]) -> List[Tuple[str, str]]:
        """
        List tree edges as (sequence, sequence) pairs.

        Parameters
        ----------
        ancestors : dict
            Sampled internal-node sequences from sample_ancestors.

        Returns
        -------
        list of tuple
            One (seq_from, seq_to) pair per tree edge.
        """

        def seq_of(node: int) -> str:
            return self.sequences[node] if node >= 0 else ancestors[node]

        return [(seq_of(u), seq_of(v)) for u, v in self.topology.edges()]

    def _pick_root(self) -> int:
        """
        Choose a root node for traversals (an internal node if any).

        Returns
        -------
        int
            Node id to root traversals at.
        """
        for node in self.topology.nodes:
            if node < 0:
                return node
        return next(iter(self.topology.nodes))


def stepwise_addition_tree(
    sequences: SequenceType[str],
    weights: SequenceType[int],
    rng: random.Random,
) -> ParsimonyTree:
    """
    Build a parsimony tree by random-order stepwise addition.

    Taxa are inserted in random order; each is grafted onto the edge
    that minimises the Fitch score, breaking ties randomly. This is the
    classic stepwise-addition heuristic and serves as PyPopART's tree
    source in place of PopART's user-supplied Nexus trees.

    Parameters
    ----------
    sequences : list of str
        Condensed sequence per haplotype index.
    weights : list of int
        Site weights.
    rng : random.Random
        Random source.

    Returns
    -------
    ParsimonyTree
        The constructed tree.
    """
    n = len(sequences)
    order = list(range(n))
    rng.shuffle(order)

    topology = nx.Graph()
    internal_counter = [0]

    def new_internal() -> int:
        internal_counter[0] -= 1
        return internal_counter[0]

    if n == 1:
        topology.add_node(order[0])
        return ParsimonyTree(topology, sequences, weights)
    if n == 2:
        topology.add_edge(order[0], order[1])
        return ParsimonyTree(topology, sequences, weights)

    # Seed with the first three taxa around one internal node
    hub = new_internal()
    for leaf in order[:3]:
        topology.add_edge(hub, leaf)

    for leaf in order[3:]:
        best_score = None
        best_edges: List[Tuple[int, int]] = []
        for u, v in list(topology.edges()):
            attach = new_internal()
            topology.remove_edge(u, v)
            topology.add_edge(u, attach)
            topology.add_edge(attach, v)
            topology.add_edge(attach, leaf)

            score = ParsimonyTree(topology, sequences, weights).compute_score()
            if best_score is None or score < best_score:
                best_score = score
                best_edges = [(u, v)]
            elif score == best_score:
                best_edges.append((u, v))

            topology.remove_node(attach)
            topology.add_edge(u, v)

        u, v = rng.choice(best_edges)
        attach = new_internal()
        topology.remove_edge(u, v)
        topology.add_edge(u, attach)
        topology.add_edge(attach, v)
        topology.add_edge(attach, leaf)

    return ParsimonyTree(topology, sequences, weights)


def sample_parsimony_trees(
    sequences: SequenceType[str],
    weights: SequenceType[int],
    n_trees: int,
    rng: random.Random,
) -> List[ParsimonyTree]:
    """
    Generate a set of stepwise-addition parsimony trees.

    Parameters
    ----------
    sequences : list of str
        Condensed sequence per haplotype index.
    weights : list of int
        Site weights.
    n_trees : int
        Number of trees to build.
    rng : random.Random
        Random source.

    Returns
    -------
    list of ParsimonyTree
        The sampled trees (scores already computed).
    """
    trees = []
    for _ in range(n_trees):
        tree = stepwise_addition_tree(sequences, weights, rng)
        tree.compute_score()
        trees.append(tree)
    return trees


__all__ = ['ParsimonyTree', 'sample_parsimony_trees', 'stepwise_addition_tree']
