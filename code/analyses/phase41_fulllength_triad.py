#!/usr/bin/env python3
"""Phase 41 (piste AA.5) - la triade est-elle vraiment perdue, ou seulement hors du segment HMM ?

Objet       : pour chaque accession dont l'alignement Pfam ne porte pas la triade complete,
              recuperer la sequence UniProt ENTIERE et chercher les residus homologues de
              E59/C113/H115 de Rv1025 par alignement par paire. Un hit HMM partiel coupe une
              protéine longue : l'absence dans l'alignement Pfam n'est pas une perte du site.
              Controle positif : des accessions dont l'alignement Pfam porte deja ECH.
Entrees     : résultats/phase40_tree/leaf_metadata.tsv (etat par segment aligne)
              data/Rv1025.faa (ancre, 155 aa)
              UniProt REST (sequences completes, telechargees dans un repertoire dedie)
Sorties     : résultats/phase41_fulllength/uniprot/*.fasta (donnees telechargees, isolees)
              résultats/phase41_fulllength/triad_fulllength.tsv
              résultats/phase41_fulllength/report.txt
Reutilisable: oui (changer ANCRE et POSITIONS) pour tout site cartographie sur un alignement Pfam.
Projet      : mtbc/clos_soumis/Rv1025, piste AA.5
Date        : 2026-11-02
"""
from __future__ import annotations

import csv
import sys
import time
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path

from Bio import Align
from Bio.Align import substitution_matrices

PROJ = Path(__file__).resolve().parent.parent
META = PROJ / "résultats" / "phase40_tree" / "leaf_metadata.tsv"
ANCHOR_FAA = PROJ / "data" / "Rv1025.faa"
OUTDIR = PROJ / "résultats" / "phase41_fulllength"
DLDIR = OUTDIR / "uniprot"

POSITIONS = {"E59": (59, "E"), "C113": (113, "C"), "H115": (115, "H")}
N_CONTROLS = 12


def require(path: Path, label: str) -> Path:
    if not path.exists():
        raise SystemExit(f"entree introuvable ({label}) : {path}")
    return path


def read_anchor() -> str:
    seq = "".join(
        l.strip() for l in require(ANCHOR_FAA, "ancre Rv1025").read_text().splitlines()
        if not l.startswith(">")
    )
    for label, (pos, aa) in POSITIONS.items():
        if seq[pos - 1] != aa:
            raise SystemExit(
                f"ANCRE FAUSSE : Rv1025[{pos}] = {seq[pos - 1]!r}, attendu {aa!r} pour {label}"
            )
    return seq


def fetch(acc: str) -> str | None:
    path = DLDIR / f"{acc}.fasta"
    if path.exists() and path.stat().st_size > 0:
        return path.read_text()
    url = f"https://rest.uniprot.org/uniprotkb/{acc}.fasta"
    for attempt in range(3):
        try:
            with urllib.request.urlopen(url, timeout=30) as fh:
                txt = fh.read().decode()
            if txt.strip():
                path.write_text(txt)
                return txt
            return None
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return None
            time.sleep(2 * (attempt + 1))
        except Exception:
            time.sleep(2 * (attempt + 1))
    return None


def seq_of(fasta: str) -> str:
    return "".join(l.strip() for l in fasta.splitlines() if l and not l.startswith(">"))


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    DLDIR.mkdir(parents=True, exist_ok=True)
    anchor = read_anchor()

    rows = list(csv.DictReader(require(META, "table des feuilles").open(), delimiter="\t"))
    byacc: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        byacc[r["accession"]].append(r)

    def union_complete(segs: list[dict]) -> bool:
        return all(any(s["motif"][i] != "-" for s in segs) for i in range(3))

    targets = [a for a, v in byacc.items() if not union_complete(v)]
    controls = [
        a for a, v in byacc.items()
        if len(v) == 1 and v[0]["motif"] == "ECH" and v[0]["domaine"] == "Eukaryota"
    ][:N_CONTROLS]

    aligner = Align.PairwiseAligner()
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aligner.mode = "local"

    report: list[str] = []
    report.append(f"Accessions a triade incomplete apres union des segments : {len(targets)}")
    report.append(f"Controles positifs (eucaryotes ECH, 1 segment) : {len(controls)}")

    out_rows = []
    for kind, accs in (("cible", targets), ("controle", controls)):
        for acc in accs:
            segs = byacc[acc]
            fasta = fetch(acc)
            if fasta is None:
                out_rows.append(
                    dict(kind=kind, accession=acc, domaine=segs[0]["domaine"],
                         n_segments=len(segs), motifs_pfam=";".join(s["motif"] for s in segs),
                         len_full="", aligned_span="", res59="?", res113="?", res115="?",
                         motif_full="?", verdict="sequence UniProt indisponible")
                )
                continue
            full = seq_of(fasta)
            aln = aligner.align(full, anchor)[0]
            # map : position de l'ancre (1-based) -> residu de la sequence complete
            amap: dict[int, str] = {}
            ta, tb = aln.aligned  # blocs apparies (full, anchor)
            for (s1, e1), (s2, e2) in zip(ta, tb):
                for k in range(e1 - s1):
                    amap[s2 + k + 1] = full[s1 + k]
            got = {}
            for label, (pos, aa) in POSITIONS.items():
                got[label] = amap.get(pos, "-")
            motif_full = got["E59"] + got["C113"] + got["H115"]
            if motif_full == "ECH":
                verdict = "TRIADE PRESENTE dans la proteine complete"
            elif all(c != "-" for c in motif_full):
                verdict = "residus presents mais substitues : " + motif_full
            else:
                verdict = "position non couverte par l'alignement par paire"
            out_rows.append(
                dict(kind=kind, accession=acc, domaine=segs[0]["domaine"], n_segments=len(segs),
                     motifs_pfam=";".join(s["motif"] for s in segs), len_full=str(len(full)),
                     aligned_span=f"{tb[0][0] + 1}-{tb[-1][1]}",
                     res59=got["E59"], res113=got["C113"], res115=got["H115"],
                     motif_full=motif_full, verdict=verdict)
            )

    fields = list(out_rows[0].keys())
    with (OUTDIR / "triad_fulllength.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    from collections import Counter
    report.append("")
    report.append("=== CONTROLE POSITIF (ces accessions portent ECH dans l'alignement Pfam) ===")
    cc = Counter(r["motif_full"] for r in out_rows if r["kind"] == "controle")
    for m, n in cc.most_common():
        report.append(f"  motif retrouve par la methode : {m}  n={n}")
    report.append("  -> la methode est valide si elle rend ECH pour (presque) tous les controles.")
    report.append("")
    report.append("=== CIBLES ===")
    tc = Counter((r["domaine"], r["motif_full"]) for r in out_rows if r["kind"] == "cible")
    for (d, m), n in sorted(tc.items()):
        report.append(f"  {d[:22]:22s} motif de la proteine complete = {m:4s}  n={n}")
    report.append("")
    nech = sum(1 for r in out_rows if r["kind"] == "cible" and r["motif_full"] == "ECH")
    ncib = sum(1 for r in out_rows if r["kind"] == "cible")
    report.append(
        f"BILAN : {nech}/{ncib} des accessions declarees sans triade portent en realite "
        f"les trois residus dans leur sequence complete."
    )
    report.append("")
    report.append("=== DETAIL DES CIBLES ===")
    for r in sorted((r for r in out_rows if r["kind"] == "cible"),
                    key=lambda r: (r["domaine"], r["motif_full"])):
        report.append(
            f"  {r['domaine'][:10]:10s} {r['accession']:12s} len={r['len_full']:>5s} "
            f"pfam={r['motifs_pfam']:12s} complete={r['motif_full']:4s} "
            f"span_ancre={r['aligned_span']:>9s}  {r['verdict']}"
        )

    txt = "\n".join(report) + "\n"
    (OUTDIR / "report.txt").write_text(txt)
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
