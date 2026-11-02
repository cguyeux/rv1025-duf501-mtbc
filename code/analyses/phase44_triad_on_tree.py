#!/usr/bin/env python3
"""Phase 44 (piste AA.5) - etat de la triade cartographie sur l'arbre, et comptage des pertes.

Objet       : poser sur l'arbre ML de DUF501 l'etat du site metal, puis compter par
              parsimonie (Fitch) le nombre MINIMAL d'evenements de perte. Le caractere est
              l'etat du site de la PROTEINE (union des segments d'une meme accession), pas
              celui du segment aligne : un domaine detecte en deux hits par le HMM donne deux
              segments a site partiel alors que la proteine porte les trois residus
              (23 accessions sur 24, phase 42). Compter par segment fabriquerait des pertes.
              Le nombre de changements de Fitch ne depend pas de l'enracinement ; la POLARITE
              (perte plutot que gain) vient de l'etat bacterien, complet chez 1487/1490.
Entrees     : résultats/phase43_reduced/align_reduced.fasta (jeu de l'arbre)
              arbre IQ-TREE (.treefile) passe en argument
              résultats/phase40_tree/leaf_metadata.tsv
Sorties     : résultats/phase44_tree_states/<prefixe>_states.tsv
              résultats/phase44_tree_states/<prefixe>_report.txt
Usage       : phase44_triad_on_tree.py <arbre.treefile> [prefixe]
Reutilisable: oui pour tout caractere discret a cartographier sur un arbre de famille.
Projet      : mtbc/clos_soumis/Rv1025, piste AA.5
Date        : 2026-11-02
"""
from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import dendropy

PROJ = Path(__file__).resolve().parent.parent
META = PROJ / "résultats" / "phase40_tree" / "leaf_metadata.tsv"
OUTDIR = PROJ / "résultats" / "phase44_tree_states"


def require(p: Path, label: str) -> Path:
    if not p.exists():
        raise SystemExit(f"entree introuvable ({label}) : {p}")
    return p


def protein_motif(segs: list[dict]) -> str:
    out = []
    for i in range(3):
        aas = {s["motif"][i] for s in segs if s["motif"][i] != "-"}
        out.append(sorted(aas)[0] if aas else "-")
    return "".join(out)


def fitch(tree: dendropy.Tree, states: dict[str, str]) -> tuple[int, dict]:
    """Nombre minimal de changements (Fitch, arbre non racine via un passage ascendant
    sur l'arbre tel qu'il est lu ; le score ne depend pas du choix de la racine)."""
    changes = 0
    sets: dict[int, set[str]] = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            lab = node.taxon.label.replace(" ", "_")
            sets[id(node)] = {states[lab]}
        else:
            kids = [sets[id(c)] for c in node.child_node_iter()]
            inter = set.intersection(*kids)
            if inter:
                sets[id(node)] = inter
            else:
                sets[id(node)] = set.union(*kids)
                changes += 1
    return changes, sets


def main() -> int:
    if len(sys.argv) < 2:
        raise SystemExit("usage : phase44_triad_on_tree.py <arbre.treefile> [prefixe]")
    treefile = require(Path(sys.argv[1]), "arbre")
    prefix = sys.argv[2] if len(sys.argv) > 2 else treefile.stem
    OUTDIR.mkdir(parents=True, exist_ok=True)

    rows = list(csv.DictReader(require(META, "table des feuilles").open(), delimiter="\t"))
    byname = {r["newick_name"]: r for r in rows}
    byacc: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        byacc[r["accession"]].append(r)

    tree = dendropy.Tree.get(path=str(treefile), schema="newick", preserve_underscores=True)
    leaves = [n.taxon.label.replace(" ", "_") for n in tree.leaf_node_iter()]
    missing = [l for l in leaves if l not in byname]
    if missing:
        raise SystemExit(f"{len(missing)} feuilles absentes de la table : {missing[:5]}")

    seg_state, prot_state = {}, {}
    for l in leaves:
        r = byname[l]
        seg_state[l] = "complet" if r["motif"] == "ECH" else "incomplet"
        pm = protein_motif(byacc[r["accession"]])
        prot_state[l] = "complet" if pm == "ECH" else "incomplet"

    n_seg, _ = fitch(tree, seg_state)
    n_prot, sets = fitch(tree, prot_state)

    report = [
        f"Arbre : {treefile}",
        f"Feuilles : {len(leaves)}",
        "",
        "=== Caractere : etat du site ===",
        f"  par SEGMENT aligne   : {Counter(seg_state.values())}",
        f"  par PROTEINE (union) : {Counter(prot_state.values())}",
        "",
        "=== Nombre MINIMAL de changements (parsimonie de Fitch) ===",
        f"  caractere par segment  : {n_seg} changements",
        f"  caractere par proteine : {n_prot} changements  <-- le seul defendable",
        "",
        "Lecture : l'ecart entre les deux comptages est le nombre d'evenements que la",
        "segmentation du HMM fabrique si l'on cartographie l'etat du segment au lieu de",
        "l'etat de la proteine.",
        "",
        "=== Feuilles a site incomplet au niveau de la PROTEINE ===",
    ]
    inc = [l for l in leaves if prot_state[l] == "incomplet"]
    for l in sorted(inc, key=lambda l: (byname[l]["domaine"], byname[l]["accession"])):
        r = byname[l]
        pm = protein_motif(byacc[r["accession"]])
        org = r["organism_lineage"].split(" | ")[0] if r["organism_lineage"] else ""
        report.append(
            f"  {l:8s} {r['domaine'][:12]:12s} {r['accession']:12s} motif_proteine={pm:4s} "
            f"occ={r['seq_occupancy']} {org[:52]}"
        )
    report.append("")
    report.append(f"TOTAL feuilles a site incomplet (proteine) : {len(inc)}")

    with (OUTDIR / f"{prefix}_states.tsv").open("w") as fh:
        fh.write("leaf\tpfam_name\taccession\tdomaine\tmotif_segment\tmotif_proteine\t"
                 "etat_segment\tetat_proteine\torganism_lineage\n")
        for l in leaves:
            r = byname[l]
            fh.write("\t".join([
                l, r["pfam_name"], r["accession"], r["domaine"], r["motif"],
                protein_motif(byacc[r["accession"]]), seg_state[l], prot_state[l],
                r.get("organism_lineage", ""),
            ]) + "\n")

    txt = "\n".join(report) + "\n"
    (OUTDIR / f"{prefix}_report.txt").write_text(txt)
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
