#!/usr/bin/env python3
"""Phase 46 (piste AA.6) - mise a jour de l'Online Resource 10 (membres eucaryotes de DUF501).

Objet       : ajouter a la table supplementaire les colonnes qui portent la correction d'AA.5 :
              nombre de segments de l'accession dans l'alignement Pfam, motif du site au niveau
              de la PROTEINE (union des segments) et verdict. L'ancienne table ne donnait que le
              motif du SEGMENT, ce qui comptait huit pertes inexistantes.
Entrees     : article/supplementary_materials/table_S10_duf501_eukaryotes.tsv (version a corriger,
              dont la colonne d'identite a Rv1025 est reprise telle quelle)
              résultats/phase40_tree/leaf_metadata.tsv
              résultats/phase42_mafft/verdict.tsv (controle sur sequences entieres)
Sorties     : article/supplementary_materials/table_S10_duf501_eukaryotes.tsv (reecrite)
Reutilisable: non (table d'un article donne).
Projet      : mtbc/clos_soumis/Rv1025, piste AA.6
Date        : 2026-11-02
"""
from __future__ import annotations

import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
TABLE = PROJ / "article" / "supplementary_materials" / "table_S10_duf501_eukaryotes.tsv"
META = PROJ / "résultats" / "phase40_tree" / "leaf_metadata.tsv"
MAFFT = PROJ / "résultats" / "phase42_mafft" / "verdict.tsv"


def require(p: Path, label: str) -> Path:
    if not p.exists():
        raise SystemExit(f"entree introuvable ({label}) : {p}")
    return p


def main() -> int:
    old = list(csv.DictReader(require(TABLE, "table S10").open(), delimiter="\t"))
    meta = list(csv.DictReader(require(META, "table des feuilles").open(), delimiter="\t"))
    mafft = {
        r["accession"]: r
        for r in csv.DictReader(require(MAFFT, "verdict MAFFT").open(), delimiter="\t")
    }

    byacc: dict[str, list[dict]] = defaultdict(list)
    for r in meta:
        byacc[r["accession"]].append(r)

    def pmotif(segs: list[dict]) -> str:
        out = []
        for i in range(3):
            aas = {s["motif"][i] for s in segs if s["motif"][i] != "-"}
            out.append(sorted(aas)[0] if aas else "-")
        return "".join(out)

    fields = ["sequence", "accession", "domaine", "longueur", "motif_E59_C113_H115",
              "n_segments_alignement", "motif_proteine_E59_C113_H115", "verdict",
              "identite_Rv1025_%"]
    rows, counts = [], Counter()
    for r in old:
        acc = r["accession"]
        segs = byacc.get(acc, [])
        if not segs:
            raise SystemExit(f"accession {acc} de la table S10 absente de la table des feuilles")
        pm = pmotif(segs)
        if r["motif_E59_C113_H115"] == "ECH":
            verdict = "site complet"
        elif pm == "ECH":
            verdict = "site complet au niveau proteine ; segment partiel (double detection HMM)"
        else:
            m = mafft.get(acc)
            verdict = (
                "perte confirmee sur la sequence entiere realignee"
                if m and m["motif_proteine_complete"] != "ECH"
                else "perte non verifiee (sequence UniProt indisponible)"
            )
        counts[verdict] += 1
        rows.append({
            "sequence": r["sequence"], "accession": acc, "domaine": r["domaine"],
            "longueur": r["longueur"], "motif_E59_C113_H115": r["motif_E59_C113_H115"],
            "n_segments_alignement": str(len(segs)),
            "motif_proteine_E59_C113_H115": pm, "verdict": verdict,
            "identite_Rv1025_%": r["identite_Rv1025_%"],
        })

    with TABLE.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"{TABLE} reecrite, {len(rows)} lignes")
    for v, n in counts.most_common():
        print(f"  {n:4d}  {v}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
