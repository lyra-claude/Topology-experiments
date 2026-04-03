#!/usr/bin/env python3
"""
Figure 1 for ECTA 2026 paper: the 8 directed graph families used in the
directed cycle count experiment. All graphs have n=8 nodes and m=16 arcs.

Edge lists are extracted directly from the Haskell source (src/Main.hs).
The Haskell Topology format is: V.fromList [predecessors_of_0, ..., predecessors_of_7],
meaning if j appears in list i, the edge is j -> i.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import networkx as nx

# ---- Graph definitions (from src/Main.hs) ----
# Format: dict mapping node -> list of predecessors (nodes that send TO this node)

GRAPHS = {
    "DAG-Layer": {
        "predecessors": {
            0: [], 1: [0], 2: [0, 1], 3: [0, 1], 4: [0, 1],
            5: [2, 3, 4], 6: [2, 3, 4], 7: [2, 3, 4],
        },
        "kappa": 0,
    },
    "DAG-Wide": {
        "predecessors": {
            0: [], 1: [0], 2: [0], 3: [0], 4: [0],
            5: [0, 1, 2, 3], 6: [0, 1, 2, 3, 4], 7: [0, 4, 5],
        },
        "kappa": 0,
    },
    "LowCyc-1": {
        "predecessors": {
            0: [5], 1: [], 2: [0, 1], 3: [0, 1], 4: [0, 1],
            5: [2, 3, 4], 6: [2, 3, 4], 7: [2, 3, 4],
        },
        "kappa": 3,
    },
    "Bidir-Ring": {
        "predecessors": {
            0: [1, 7], 1: [0, 2], 2: [1, 3], 3: [2, 4],
            4: [3, 5], 5: [4, 6], 6: [5, 7], 7: [0, 6],
        },
        "kappa": 10,
    },
    "Two-Cliques": {
        "predecessors": {
            0: [2, 3], 1: [0, 3], 2: [0, 1], 3: [1, 2],
            4: [6, 7], 5: [4, 7], 6: [4, 5], 7: [5, 6],
        },
        "kappa": 14,
    },
    "Mesh-Cyclic": {
        "predecessors": {
            0: [3, 4], 1: [0, 5], 2: [1, 6], 3: [2, 7],
            4: [0, 7], 5: [1, 4], 6: [2, 5], 7: [3, 6],
        },
        "kappa": 20,
    },
    "Dense-Tri": {
        "predecessors": {
            0: [2, 7], 1: [0, 7], 2: [1, 4], 3: [1, 2],
            4: [3, 6], 5: [3, 4], 6: [0, 5], 7: [5, 6],
        },
        "kappa": 29,
    },
    "Ring+Skip2": {
        "predecessors": {
            0: [6, 7], 1: [0, 7], 2: [0, 1], 3: [1, 2],
            4: [2, 3], 5: [3, 4], 6: [4, 5], 7: [5, 6],
        },
        "kappa": 47,
    },
}


def build_digraph(pred_dict):
    """Build a networkx DiGraph from a predecessor dictionary."""
    G = nx.DiGraph()
    G.add_nodes_from(range(8))
    for target, sources in pred_dict.items():
        for source in sources:
            G.add_edge(source, target)
    return G


def get_layout(name, G):
    """Return node positions. Use structure-aware layouts for DAGs,
    circular layout for cyclic graphs."""
    # Circular positions (default for ring-based topologies)
    r_circ = 1.1
    angles = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, 8, endpoint=False)
    circle = {i: (r_circ * np.cos(a), r_circ * np.sin(a)) for i, a in enumerate(angles)}

    if name == "DAG-Layer":
        # 3-layer layout: {0,1} top, {2,3,4} middle, {5,6,7} bottom
        # Wider spacing to reduce edge crossings
        return {
            0: (-0.4, 1.1), 1: (0.4, 1.1),
            2: (-1.1, 0.0), 3: (0.0, 0.0), 4: (1.1, 0.0),
            5: (-1.1, -1.1), 6: (0.0, -1.1), 7: (1.1, -1.1),
        }
    elif name == "DAG-Wide":
        # Fan-out from node 0, more spread out
        return {
            0: (0.0, 1.2),
            1: (-1.3, 0.3), 2: (-0.5, 0.3), 3: (0.5, 0.3), 4: (1.3, 0.3),
            5: (-0.9, -0.9), 6: (0.0, -0.9), 7: (0.9, -0.9),
        }
    elif name == "LowCyc-1":
        # Same as DAG-Layer (it's a modified DAG-Layer with 5->0 feedback)
        return {
            0: (-0.4, 1.1), 1: (0.4, 1.1),
            2: (-1.1, 0.0), 3: (0.0, 0.0), 4: (1.1, 0.0),
            5: (-1.1, -1.1), 6: (0.0, -1.1), 7: (1.1, -1.1),
        }
    elif name == "Two-Cliques":
        # Two 4-cliques clearly separated
        off = 0.65
        r = 0.45
        angs = np.linspace(np.pi / 2, np.pi / 2 - 2 * np.pi, 4, endpoint=False)
        left = {i: (-off + r * np.cos(a), r * np.sin(a)) for i, a in enumerate(angs)}
        right = {4 + i: (off + r * np.cos(a), r * np.sin(a)) for i, a in enumerate(angs)}
        return {**left, **right}
    elif name == "Mesh-Cyclic":
        # 2x4 mesh with more vertical spacing
        return {
            0: (-1.2, 0.6), 1: (-0.4, 0.6), 2: (0.4, 0.6), 3: (1.2, 0.6),
            4: (-1.2, -0.6), 5: (-0.4, -0.6), 6: (0.4, -0.6), 7: (1.2, -0.6),
        }
    else:
        # Circular for ring-based topologies
        return circle


def has_reverse_edge(G, u, v):
    """Check if both u->v and v->u exist."""
    return G.has_edge(v, u)


def draw_graph(ax, name, G, pos, kappa):
    """Draw a single directed graph on the given axes."""
    # Track which reciprocal pairs we've already drawn
    drawn_pairs = set()

    for u, v in G.edges():
        x0, y0 = pos[u]
        x1, y1 = pos[v]

        dx = x1 - x0
        dy = y1 - y0
        dist = np.sqrt(dx**2 + dy**2)
        if dist < 1e-6:
            continue

        # Node radius for offset
        r = 0.13
        sx = x0 + r * dx / dist
        sy = y0 + r * dy / dist
        ex = x1 - r * dx / dist
        ey = y1 - r * dy / dist

        # Use curvature to separate reciprocal edges
        pair = (min(u, v), max(u, v))
        if has_reverse_edge(G, u, v):
            rad = 0.2
        else:
            rad = 0.1

        ax.annotate(
            "",
            xy=(ex, ey),
            xytext=(sx, sy),
            arrowprops=dict(
                arrowstyle="-|>",
                color="#555555",
                lw=0.7,
                shrinkA=0,
                shrinkB=0,
                connectionstyle=f"arc3,rad={rad}",
                mutation_scale=8,
            ),
        )

    # Draw nodes
    for node, (x, y) in pos.items():
        circle = plt.Circle(
            (x, y), 0.13, facecolor="white", edgecolor="black",
            linewidth=1.2, zorder=5,
        )
        ax.add_patch(circle)
        ax.text(
            x, y, str(node), ha="center", va="center",
            fontsize=6, fontweight="bold", zorder=6,
        )

    # Label below the graph
    ax.text(
        0.5, -0.02, f"{name}  ($\\kappa = {kappa}$)",
        transform=ax.transAxes, ha="center", va="top",
        fontsize=7.5, fontweight="bold",
    )
    ax.set_aspect("equal")
    ax.set_xlim(-1.7, 1.7)
    ax.set_ylim(-1.6, 1.5)
    ax.axis("off")


def main():
    fig, axes = plt.subplots(2, 4, figsize=(7.0, 4.2))
    fig.subplots_adjust(wspace=0.12, hspace=0.45, left=0.02, right=0.98,
                        top=0.90, bottom=0.04)

    names = list(GRAPHS.keys())
    for idx, name in enumerate(names):
        row, col = divmod(idx, 4)
        ax = axes[row, col]

        info = GRAPHS[name]
        G = build_digraph(info["predecessors"])
        pos = get_layout(name, G)
        kappa = info["kappa"]

        # Verify edge count
        assert G.number_of_edges() == 16, (
            f"{name}: expected 16 edges, got {G.number_of_edges()}"
        )
        assert G.number_of_nodes() == 8, (
            f"{name}: expected 8 nodes, got {G.number_of_nodes()}"
        )

        draw_graph(ax, name, G, pos, kappa)

    fig.suptitle(
        "Eight directed graph families ($n=8$, $m=16$), ordered by cycle count $\\kappa(D)$",
        fontsize=9, fontweight="bold", y=0.97,
    )

    # Save
    out_dir = "/home/lyra/projects/Topology-experiments/results/directed_full/figures"
    fig.savefig(f"{out_dir}/fig0_graph_families.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(f"{out_dir}/fig0_graph_families.png", dpi=300, bbox_inches="tight")
    print(f"Saved to {out_dir}/fig0_graph_families.pdf and .png")
    plt.close()


if __name__ == "__main__":
    main()
