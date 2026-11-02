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
STO = PROJ / "data" / "PF04417_full.sto"
NEAR_COMPLETE = 0.90   # seuil de quasi-complétude, identique à phase29 (90 % de la médiane)
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


def aligned_lengths() -> dict[str, int]:
    """Résidus alignés (hors gaps) de chaque segment de l'alignement Pfam."""
    out = {}
    for line in STO.open():
        if line.startswith(("#", "//")) or not line.strip():
            continue
        name, seq = line.split(None, 1)
        out[name] = sum(1 for ch in seq.strip() if ch not in "-.")
    return out


WILD = {"complet", "perdu"}   # état « manquant » : compatible avec les deux


def fitch(tree: dendropy.Tree, states: dict[str, str]) -> tuple[int, dict]:
    """Nombre minimal de changements (Fitch, arbre non racine via un passage ascendant
    sur l'arbre tel qu'il est lu ; le score ne depend pas du choix de la racine)."""
    changes = 0
    sets: dict[int, set[str]] = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            lab = node.taxon.label.replace(" ", "_")
            st = states[lab]
            sets[id(node)] = set(WILD) if st == "manquant" else {st}
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

    alen = aligned_lengths()
    tri = sorted(alen[r["pfam_name"]] for r in rows if r["triade_stricte"] == "oui")
    m = len(tri)
    mediane = tri[m // 2] if m % 2 else (tri[m // 2 - 1] + tri[m // 2]) / 2
    seuil = NEAR_COMPLETE * mediane

    # trois états au niveau de la PROTEINE : complet ; perdu = site incomplet sur une protéine
    # quasi complète (au moins un segment aligné >= seuil) ; manquant = site incomplet sur une
    # protéine dont le domaine n'est couvert que partiellement (fragment, donnée manquante).
    seg_state, prot2, prot_state = {}, {}, {}
    for l in leaves:
        r = byname[l]
        seg_state[l] = "complet" if r["motif"] == "ECH" else "incomplet"
        pm = protein_motif(byacc[r["accession"]])
        prot2[l] = "complet" if pm == "ECH" else "incomplet"
        if pm == "ECH":
            prot_state[l] = "complet"
        else:
            nc = max(alen[x["pfam_name"]] for x in byacc[r["accession"]]) >= seuil
            prot_state[l] = "perdu" if nc else "manquant"

    n_seg, _ = fitch(tree, seg_state)
    n_prot2, _ = fitch(tree, prot2)
    n_prot, sets = fitch(tree, prot_state)
    pdm = tree.phylogenetic_distance_matrix()
    lost = [l for l in leaves if prot_state[l] == "perdu"]
    tax = {n.taxon.label.replace(" ", "_"): n.taxon for n in tree.leaf_node_iter()}
    topo = {(a, b): pdm.path_edge_count(tax[a], tax[b]) for a in lost for b in lost if a < b}

    report = [
        f"Arbre : {treefile}",
        f"Feuilles : {len(leaves)}",
        "",
        "=== Caractere : etat du site ===",
        f"  par SEGMENT aligne   : {Counter(seg_state.values())}",
        f"  par PROTEINE, 2 etats (union des segments)  : {Counter(prot2.values())}",
        f"  par PROTEINE, 3 etats (complet/perdu/manquant) : {Counter(prot_state.values())}",
        f"  seuil de quasi-completude : {seuil:.1f} residus alignes ({NEAR_COMPLETE:.0%} de la mediane {mediane} des membres a triade)",
        "",
        "=== Nombre MINIMAL de changements (parsimonie de Fitch) ===",
        f"  segment, 2 etats                 : {n_seg} changements  (artefact de segmentation)",
        f"  proteine, 2 etats                : {n_prot2} changements  (compte les fragments comme des pertes)",
        f"  proteine, 3 etats, manquant = ?  : {n_prot} changements  <-- le seul defendable",
        "",
        "Lecture : les deux premiers comptages fabriquent des evenements, le premier par la",
        "segmentation du HMM, le second en confondant donnee manquante et perte. Le dernier ne",
        "retient que les proteines quasi completes dont le site est incomplet.",
        "",
        "=== Distance topologique (aretes) entre feuilles a perte confirmee ===",
        *[f"  {a} - {b} : {d} aretes" for (a, b), d in sorted(topo.items())],
        "",
        "=== Feuilles a perte confirmee (site incomplet, proteine quasi complete) ===",
    ]
    inc = lost
    for l in sorted(inc, key=lambda l: (byname[l]["domaine"], byname[l]["accession"])):
        r = byname[l]
        pm = protein_motif(byacc[r["accession"]])
        org = r["organism_lineage"].split(" | ")[0] if r["organism_lineage"] else ""
        report.append(
            f"  {l:8s} {r['domaine'][:12]:12s} {r['accession']:12s} motif_proteine={pm:4s} "
            f"occ={r['seq_occupancy']} {org[:52]}"
        )
    report.append("")
    report.append(f"TOTAL feuilles a perte confirmee : {len(inc)} ; manquantes (fragments) : "
                  f"{sum(1 for v in prot_state.values() if v == 'manquant')}")

    with (OUTDIR / f"{prefix}_states.tsv").open("w") as fh:
        fh.write("leaf\tpfam_name\taccession\tdomaine\tmotif_segment\tmotif_proteine\t"
                 "etat_segment\tetat_proteine\torganism_lineage\tetat_site\n")
        for l in leaves:
            r = byname[l]
            fh.write("\t".join([
                l, r["pfam_name"], r["accession"], r["domaine"], r["motif"],
                protein_motif(byacc[r["accession"]]), seg_state[l], prot2[l],
                r.get("organism_lineage", ""), prot_state[l],
            ]) + "\n")

    txt = "\n".join(report) + "\n"
    (OUTDIR / f"{prefix}_report.txt").write_text(txt)
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
