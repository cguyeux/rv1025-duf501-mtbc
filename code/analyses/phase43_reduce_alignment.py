#!/usr/bin/env python3
"""Phase 43 (piste AA.5) - jeu reduit non redondant pour un arbre ML interpretable.

Objet       : l'alignement du domaine ne porte que 129 colonnes pour 2142 sequences, soit
              4281 branches a estimer sur 129 sites : l'arbre complet est structurellement
              sous-determine (AICc de ModelFinder au-dela de 3e7, cf. duf501_core.log).
              On reduit la redondance par agregation gloutonne a seuil d'identite, en
              GARDANT d'office les sequences qui portent le signal a cartographier.
Entrees     : résultats/phase40_tree/align_tree.fasta (2142 x 129, toutes ECH)
              résultats/phase40_tree/align_full.fasta (2228 x 129, avec les lacunaires)
              résultats/phase40_tree/leaf_metadata.tsv
Sorties     : résultats/phase43_reduced/align_reduced.fasta
              résultats/phase43_reduced/clusters.tsv, report.txt
Reutilisable: oui (parametre IDENTITY, liste des sequences a garder d'office).
Projet      : mtbc/clos_soumis/Rv1025, piste AA.5
Date        : 2026-11-02
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

PROJ = Path(__file__).resolve().parent.parent
TREEDIR = PROJ / "résultats" / "phase40_tree"
OUTDIR = PROJ / "résultats" / "phase43_reduced"
IDENTITY = float(sys.argv[1]) if len(sys.argv) > 1 else 0.90  # seuil d'agregation
MIN_SHARED = 60          # colonnes non-gap communes exigees pour comparer deux sequences


def require(p: Path, label: str) -> Path:
    if not p.exists():
        raise SystemExit(f"entree introuvable ({label}) : {p}")
    return p


def read_fasta(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    name = None
    for line in require(path, path.name).read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].strip()
            seqs[name] = ""
        elif name:
            seqs[name] += line.strip()
    return seqs


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    seqs = read_fasta(TREEDIR / "align_full.fasta")
    meta = {
        r["newick_name"]: r
        for r in csv.DictReader(require(TREEDIR / "leaf_metadata.tsv", "meta").open(), delimiter="\t")
    }

    names = sorted(seqs)
    ncol = len(seqs[names[0]])
    arr = np.frombuffer("".join(seqs[n] for n in names).encode(), dtype="S1").reshape(len(names), ncol)
    gap = arr == b"-"

    # sequences gardees d'office : l'ancre, tous les eucaryotes, toutes les non-ECH,
    # et tout membre d'une accession a plusieurs segments (les cas a double detection).
    from collections import Counter
    acc_count = Counter(r["accession"] for r in meta.values())
    forced = set()
    for n, r in meta.items():
        if r["pfam_name"].startswith("P96375_MYCTU"):
            forced.add(n)
        if r["domaine"] == "Eukaryota" or r["triade_stricte"] == "non" or acc_count[r["accession"]] > 1:
            forced.add(n)

    order = [n for n in names if n in forced] + [n for n in names if n not in forced]
    idx = {n: i for i, n in enumerate(names)}

    reps: list[str] = []
    assigned: dict[str, str] = {}
    for n in order:
        if n in assigned:
            continue
        i = idx[n]
        if reps:
            ri = np.array([idx[r] for r in reps])
            shared = (~gap[i]) & (~gap[ri])
            nshared = shared.sum(axis=1)
            same = ((arr[i] == arr[ri]) & shared).sum(axis=1)
            with np.errstate(invalid="ignore", divide="ignore"):
                ident = np.where(nshared >= MIN_SHARED, same / np.maximum(nshared, 1), 0.0)
            best = int(np.argmax(ident))
            if ident[best] >= IDENTITY and n not in forced:
                assigned[n] = reps[best]
                continue
        reps.append(n)
        assigned[n] = n

    with (OUTDIR / "align_reduced.fasta").open("w") as fh:
        for n in reps:
            fh.write(f">{n}\n{seqs[n]}\n")
    with (OUTDIR / "clusters.tsv").open("w") as fh:
        fh.write("sequence\trepresentant\tpfam_name\tdomaine\tmotif\n")
        for n in names:
            r = meta[n]
            fh.write(f"{n}\t{assigned[n]}\t{r['pfam_name']}\t{r['domaine']}\t{r['motif']}\n")

    dom = Counter(meta[n]["domaine"] for n in reps)
    mot = Counter(meta[n]["motif"] for n in reps)
    report = [
        f"Entree : {len(names)} sequences x {ncol} colonnes",
        f"Seuil d'identite : {IDENTITY:.0%} ; colonnes communes minimales : {MIN_SHARED}",
        f"Sequences gardees d'office : {len(forced)} "
        f"(ancre Rv1025, eucaryotes, non-ECH, accessions a plusieurs segments)",
        f"Representants retenus : {len(reps)}",
        f"Branches a estimer : {2 * len(reps) - 3} pour {ncol} sites "
        f"(ratio sites/branches = {ncol / (2 * len(reps) - 3):.2f})",
        "",
        "Composition du jeu reduit :",
    ]
    for d, n in dom.most_common():
        report.append(f"  {d:28s} {n:5d}")
    report.append("")
    report.append("Motifs de la triade dans le jeu reduit :")
    for m, n in mot.most_common():
        report.append(f"  {m:5s} {n:5d}")
    txt = "\n".join(report) + "\n"
    (OUTDIR / "report.txt").write_text(txt)
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
