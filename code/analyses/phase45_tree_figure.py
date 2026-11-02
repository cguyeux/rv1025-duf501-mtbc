#!/usr/bin/env python3
"""Phase 45 (piste AA.5) - figure : arbre ML de DUF501, etat du site metal aux feuilles.

Objet       : cladogramme circulaire de la famille DUF501 avec deux anneaux, le domaine du
              vivant et l'etat du site metal a trois etats au niveau de la PROTEINE (complet, perdu sur une
              proteine quasi complete, manquant pour les fragments),
              les trois pertes confirmees et l'ancre Rv1025 reperees. Cladogramme et non
              phylogramme : le domaine ne fournit que 129 colonnes, les longueurs de branche
              profondes ne sont pas interpretables (cf. report de la phase 43).
Entrees     : arbre IQ-TREE (.treefile) ; résultats/phase44_tree_states/<prefixe>_states.tsv
Sorties     : article/figures/duf501_tree_triad.pdf et .png (600 dpi)
Usage       : phase45_tree_figure.py <arbre.treefile> <states.tsv> [sortie_sans_extension]
Reutilisable: oui pour tout arbre de famille a etat discret aux feuilles.
Projet      : mtbc/clos_soumis/Rv1025, piste AA.5
Date        : 2026-11-02
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import dendropy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Wedge

PROJ = Path(__file__).resolve().parent.parent
# palette Okabe-Ito (lisible en vision des couleurs atypique)
COL_DOM = {
    "Bacteria": "#0072B2",
    "Eukaryota": "#D55E00",
    "supprimé d'UniProt": "#BBBBBB",
    "non renseigne": "#EEEEEE",
    "non retourné par UniProt": "#EEEEEE",
}
COL_STATE = {"complet": "#009E73", "perdu": "#CC79A7", "manquant": "#CCCCCC"}
LOSSES = {"A0A1D3CV76", "U6GA88", "A0A2C6L2K9"}


def require(p: Path, label: str) -> Path:
    if not p.exists():
        raise SystemExit(f"entree introuvable ({label}) : {p}")
    return p


def main() -> int:
    if len(sys.argv) < 3:
        raise SystemExit("usage : phase45_tree_figure.py <arbre.treefile> <states.tsv> [sortie]")
    treefile = require(Path(sys.argv[1]), "arbre")
    states_tsv = require(Path(sys.argv[2]), "etats")
    out = Path(sys.argv[3]) if len(sys.argv) > 3 else PROJ / "article" / "figures" / "duf501_tree_triad"
    out.parent.mkdir(parents=True, exist_ok=True)

    st = {r["leaf"]: r for r in csv.DictReader(states_tsv.open(), delimiter="\t")}

    tree = dendropy.Tree.get(path=str(treefile), schema="newick", preserve_underscores=True)
    try:
        tree.reroot_at_midpoint(update_bipartitions=False)
    except Exception:
        pass  # arbre non enracinable au point median : on garde la lecture telle quelle

    leaves = [n for n in tree.leaf_node_iter()]
    n = len(leaves)
    if n == 0:
        raise SystemExit("arbre sans feuille")

    # angle par feuille, dans l'ordre du parcours
    ang = {}
    for i, lf in enumerate(leaves):
        ang[id(lf)] = 2 * math.pi * i / n
    for node in tree.postorder_node_iter():
        if not node.is_leaf():
            ks = [ang[id(c)] for c in node.child_node_iter()]
            ang[id(node)] = sum(ks) / len(ks)

    # rayon : profondeur en nombre de noeuds, normalisee (cladogramme)
    depth = {}
    for node in tree.preorder_node_iter():
        depth[id(node)] = 0 if node.parent_node is None else depth[id(node.parent_node)] + 1
    maxd = max(depth.values()) or 1
    R = {k: 0.18 + 0.62 * (v / maxd) for k, v in depth.items()}

    fig = plt.figure(figsize=(180 / 25.4, 180 / 25.4))
    ax = fig.add_axes((0.02, 0.02, 0.96, 0.96))
    ax.set_aspect("equal")
    ax.axis("off")

    # branches
    for node in tree.preorder_node_iter():
        if node.parent_node is None:
            continue
        p = node.parent_node
        a0, a1 = ang[id(p)], ang[id(node)]
        r0, r1 = R[id(p)], R[id(node)]
        # arc a rayon du parent
        aa = [a0 + (a1 - a0) * t / 24 for t in range(25)]
        ax.plot([r0 * math.cos(t) for t in aa], [r0 * math.sin(t) for t in aa],
                color="#555555", lw=0.25, solid_capstyle="butt", zorder=1)
        # segment radial
        ax.plot([r0 * math.cos(a1), r1 * math.cos(a1)], [r0 * math.sin(a1), r1 * math.sin(a1)],
                color="#555555", lw=0.25, solid_capstyle="butt", zorder=1)

    # anneaux
    dtheta = 360.0 / n
    for lf in leaves:
        lab = lf.taxon.label.replace(" ", "_")
        r = st.get(lab)
        if r is None:
            continue
        th = math.degrees(ang[id(lf)])
        ax.add_patch(Wedge((0, 0), 0.88, th - dtheta / 2, th + dtheta / 2, width=0.045,
                           facecolor=COL_DOM.get(r["domaine"], "#EEEEEE"), edgecolor="none", zorder=2))
        ax.add_patch(Wedge((0, 0), 0.94, th - dtheta / 2, th + dtheta / 2, width=0.045,
                           facecolor=COL_STATE[r["etat_site"]], edgecolor="none", zorder=2))
        if r["accession"] in LOSSES:
            a = ang[id(lf)]
            ax.plot([1.00 * math.cos(a)], [1.00 * math.sin(a)], marker="v", ms=4.5,
                    color="#CC79A7", mec="black", mew=0.4, zorder=4)
            nom = (r["organism_lineage"].split(" | ")[0] if r.get("organism_lineage") else r["accession"])
            nom = nom.split(" (")[0]                      # retire le nom vernaculaire entre parentheses
            k = sum(1 for x in ax.texts if getattr(x, "_perte", False))
            t = ax.annotate(nom, xy=(1.0 * math.cos(a), 1.0 * math.sin(a)),
                            xytext=(-0.62, -1.08 - 0.07 * k), fontsize=6, style="italic",
                            ha="right", va="center",
                            arrowprops=dict(arrowstyle="-", lw=0.4, color="#CC79A7",
                                            relpos=(1.0, 0.5), shrinkA=1, shrinkB=2))
            t._perte = True
        if r["pfam_name"].startswith("P96375_MYCTU"):
            a = ang[id(lf)]
            ax.plot([1.00 * math.cos(a)], [1.00 * math.sin(a)], marker="*", ms=9,
                    color="#F0E442", mec="black", mew=0.5, zorder=5)
            ax.text(1.07 * math.cos(a), 1.07 * math.sin(a), "Rv1025 (M. tuberculosis)",
                    fontsize=6.5, fontweight="bold",
                    ha="left" if math.cos(a) >= 0 else "right", va="center",
                    rotation=math.degrees(a) if math.cos(a) >= 0 else math.degrees(a) + 180,
                    rotation_mode="anchor")

    handles = [
        Line2D([], [], marker="s", ls="", ms=7, color=COL_DOM["Bacteria"], label="Bacteria"),
        Line2D([], [], marker="s", ls="", ms=7, color=COL_DOM["Eukaryota"], label="Eukaryota"),
        Line2D([], [], marker="s", ls="", ms=7, color=COL_DOM["supprimé d'UniProt"],
               label="lineage unassigned"),
        Line2D([], [], marker="s", ls="", ms=7, color=COL_STATE["complet"],
               label="Cys-His-Glu site complete"),
        Line2D([], [], marker="v", ls="", ms=6, color=COL_STATE["perdu"], mec="black",
               label="site lost (full-length check)"),
        Line2D([], [], marker="s", ls="", ms=7, color=COL_STATE["manquant"],
               label="site unknown (domain fragment)"),
        Line2D([], [], marker="*", ls="", ms=10, color="#F0E442", mec="black", label="Rv1025"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.02, 1.0), frameon=False,
              fontsize=6.5, handletextpad=0.5, labelspacing=0.5)
    ax.set_xlim(-1.32, 1.55)
    ax.set_ylim(-1.26, 1.08)

    fig.savefig(str(out) + ".pdf")
    fig.savefig(str(out) + ".png", dpi=600)
    plt.close(fig)
    print(f"ecrit : {out}.pdf et {out}.png ({n} feuilles)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
