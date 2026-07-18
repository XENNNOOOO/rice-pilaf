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


def induce_subnetwork(G, gene_ids):
    """
    Returns the subgraph of G induced by the given set of genes,
    restricted to genes that are actually present in the network

    Parameters:
    - G: NetworkX graph representation of the full PPI network
    - gene_ids: Accessions of the query genes

    Returns:
    - Induced subgraph
    - Set of query genes that were not found in the network
    """
    if G is None:
        return None, set(gene_ids)

    present = [gene for gene in gene_ids if gene in G]
    unrecognized = set(gene_ids) - set(present)

    return G.subgraph(present).copy(), unrecognized


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
    Ranks genes in the (sub)network by degree centrality to identify hub genes

    Parameters:
    - G: NetworkX graph
    - top_n: Number of top hub genes to return

    Returns:
    - List of (gene, degree) tuples, sorted by degree in descending order
    """
    if G is None or G.number_of_nodes() == 0:
        return []

    degrees = sorted(dict(G.degree()).items(), key=lambda x: x[1], reverse=True)

    return degrees[:top_n]


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