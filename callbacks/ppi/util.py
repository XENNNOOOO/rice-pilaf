import os
import pickle

import networkx as nx

from ..constants import Constants

# ======================
# PPI Networks
# ======================

PPI_NETWORKS_VALUE_LABEL = [
    {
        "value": "STRING-Physical",
        "label": "STRING Network",
        "label_id": "string",
    },
]


def get_user_facing_network(network):
    """
    Returns the user-facing label of a PPI network given its internal value

    Parameters:
    - network: Internal value of the PPI network (e.g., "STRING")

    Returns:
    - User-facing label of the PPI network (e.g., "STRING Network")
    """
    for entry in PPI_NETWORKS_VALUE_LABEL:
        if entry["value"] == network:
            return entry["label"]

    return network


# =========================================================
# Utility functions for PPI network construction/statistics
# =========================================================


def get_network_file(network):
    return f"{Constants.NETWORKS}/{network}.txt"


def load_ppi_network(network):
    """
    Loads the full PPI network (e.g., STRING) as a NetworkX graph

    Parameters:
    - network: PPI network (e.g., "STRING")

    Returns:
    - NetworkX graph representation of the network, or None if the
      network file does not exist or could not be parsed
    """
    network_file = get_network_file(network)

    try:
        return nx.read_weighted_edgelist(network_file)
    except FileNotFoundError:
        return None
    except Exception:
        # Fall back to an unweighted edge list in case the network file
        # does not carry edge weights
        try:
            return nx.read_edgelist(network_file)
        except Exception:
            return None


_UNIPROT_TO_MSU_CACHE = {}


def get_uniprot_to_msu_mapping():
    """
    Loads (and caches in-process) the UniProt accession -> [MSU gene ID, ...]
    mapping built for the PPI module.
    """
    if "mapping" in _UNIPROT_TO_MSU_CACHE:
        return _UNIPROT_TO_MSU_CACHE["mapping"]

    path = f"{Constants.MSU_MAPPING}/uniprot_to_msu.pickle"
    if not os.path.exists(path):
        _UNIPROT_TO_MSU_CACHE["mapping"] = {}
        return {}

    with open(path, "rb") as f:
        mapping = pickle.load(f)

    _UNIPROT_TO_MSU_CACHE["mapping"] = mapping
    return mapping


def get_msu_to_uniprot_mapping():
    """
    Builds (and caches) the reverse index: MSU gene ID -> [UniProt accession, ...]
    The raw STRING network is labeled by UniProt accession, while query genes
    coming from the rest of RicePilaf are MSU IDs -- this reverse index is what
    lets us check whether a query gene is actually present in the network.
    """
    if "reverse" in _UNIPROT_TO_MSU_CACHE:
        return _UNIPROT_TO_MSU_CACHE["reverse"]

    uniprot_to_msu = get_uniprot_to_msu_mapping()
    msu_to_uniprot = {}
    for uniprot_acc, msu_genes in uniprot_to_msu.items():
        for gene in msu_genes:
            msu_to_uniprot.setdefault(gene, []).append(uniprot_acc)

    _UNIPROT_TO_MSU_CACHE["reverse"] = msu_to_uniprot
    return msu_to_uniprot


def induce_subnetwork(G, gene_ids):
    """
    Returns the subgraph of G induced by the given set of MSU genes,
    restricted to genes that are actually present in the network.

    The network itself is labeled by UniProt accession, so query genes
    (MSU IDs) are first translated to their corresponding UniProt
    accession(s) before checking membership.

    Parameters:
    - G: NetworkX graph representation of the full PPI network (UniProt-labeled)
    - gene_ids: MSU accessions of the query genes

    Returns:
    - Induced subgraph (still UniProt-labeled)
    - Set of query MSU genes that were not found in the network (either no
      UniProt mapping exists, or none of their mapped accessions are in G)
    """
    if G is None:
        return None, set(gene_ids)

    msu_to_uniprot = get_msu_to_uniprot_mapping()

    present_accessions = []
    unrecognized = set()
    for gene in gene_ids:
        accessions = [
            acc for acc in msu_to_uniprot.get(gene, []) if acc in G
        ]
        if accessions:
            present_accessions.extend(accessions)
        else:
            unrecognized.add(gene)

    return G.subgraph(present_accessions).copy(), unrecognized


def compute_network_stats(G):
    """
    Computes topological summary statistics for a PPI (sub)network

    Parameters:
    - G: NetworkX graph

    Returns:
    - Dictionary containing the node count, edge count, network density,
      and average node degree
    """
    if G is None or G.number_of_nodes() == 0:
        return {
            "num_nodes": 0,
            "num_edges": 0,
            "density": 0.0,
            "avg_degree": 0.0,
        }

    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()
    density = nx.density(G)
    avg_degree = (
        sum(dict(G.degree()).values()) / num_nodes if num_nodes > 0 else 0.0
    )

    return {
        "num_nodes": num_nodes,
        "num_edges": num_edges,
        "density": density,
        "avg_degree": avg_degree,
    }


def get_hub_genes(G, top_n=10):
    """
    Ranks genes in the (sub)network by degree centrality to identify hub genes.
    The network is UniProt-labeled internally; each hub is translated back to
    its MSU gene ID(s) for display, since that's what the rest of RicePilaf
    (and the person using it) works with.

    Parameters:
    - G: NetworkX graph (UniProt-labeled)
    - top_n: Number of top hub genes to return

    Returns:
    - List of (display_label, degree) tuples, sorted by degree in descending
      order. display_label is one or more MSU gene IDs (slash-separated if a
      single UniProt accession maps to more than one), falling back to the
      raw UniProt accession if no MSU mapping is available.
    """
    if G is None or G.number_of_nodes() == 0:
        return []

    uniprot_to_msu = get_uniprot_to_msu_mapping()
    degrees = sorted(dict(G.degree()).items(), key=lambda x: x[1], reverse=True)

    hub_genes = []
    for uniprot_acc, degree in degrees[:top_n]:
        msu_genes = uniprot_to_msu.get(uniprot_acc)
        label = "/".join(msu_genes) if msu_genes else uniprot_acc
        hub_genes.append((label, degree))

    return hub_genes


def get_query_network_summary(network, gene_ids):
    """
    Builds the subnetwork induced by the query genes and computes its
    topological statistics and hub genes in one pass. Used to populate
    the module's statistics panel and to defensively handle empty/invalid
    gene lists.

    Parameters:
    - network: PPI network (e.g., "STRING")
    - gene_ids: Accessions of the query genes

    Returns:
    - Dictionary with keys: stats, hub_genes, unrecognized_genes, is_empty
    """
    if not gene_ids:
        return {
            "stats": compute_network_stats(None),
            "hub_genes": [],
            "unrecognized_genes": set(),
            "is_empty": True,
        }

    G = load_ppi_network(network)
    subG, unrecognized = induce_subnetwork(G, gene_ids)

    stats = compute_network_stats(subG)
    hub_genes = get_hub_genes(subG)

    is_empty = stats["num_nodes"] == 0

    return {
        "stats": stats,
        "hub_genes": hub_genes,
        "unrecognized_genes": unrecognized,
        "is_empty": is_empty,
    }