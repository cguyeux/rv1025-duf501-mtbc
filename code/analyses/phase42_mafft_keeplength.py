#!/usr/bin/env python3
"""Phase 42 (piste AA.5) - test decisif : la triade manque-t-elle a la PROTEINE, ou au hit HMM ?

Objet       : ajouter les sequences UniProt ENTIERES au profil Pfam par MAFFT --add
              --keeplength (les 471 colonnes de reference sont conservees, donc les colonnes
              172/395/397 restent lisibles), puis relire l'etat de la triade. Les cibles sont
              les accessions dont l'alignement Pfam ne porte pas la triade complete ; les
              controles sont des accessions dont le hit Pfam porte deja ECH : si la methode
              leur rend ECH, elle est valide.
Entrees     : data/PF04417_full.sto ; résultats/phase40_tree/leaf_metadata.tsv
              résultats/phase41_fulllength/uniprot/*.fasta (cache UniProt de la phase 41)
Sorties     : résultats/phase42_mafft/ref471.fasta, add.fasta (envoyes sur mp)
              résultats/phase42_mafft/added.fasta (retour de MAFFT)
              résultats/phase42_mafft/verdict.tsv, report.txt
Modes       : prep (ecrit les entrees) ; read (lit le retour de MAFFT)
Reutilisable: oui, pour tout site cartographie sur un alignement Pfam dont les membres sont
              des hits partiels de proteines plus longues.
Projet      : mtbc/clos_soumis/Rv1025, piste AA.5
Date        : 2026-11-02
"""
from __future__ import annotations

import csv
import sys
import time
import urllib.error
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
STO = PROJ / "data" / "PF04417_full.sto"
META = PROJ / "résultats" / "phase40_tree" / "leaf_metadata.tsv"
CACHE = PROJ / "résultats" / "phase41_fulllength" / "uniprot"
OUTDIR = PROJ / "résultats" / "phase42_mafft"

COLS = {"E59": (172, "E"), "C113": (395, "C"), "H115": (397, "H")}  # 0-based, 471 colonnes
ANCHOR = "P96375_MYCTU/16-137"
N_CONTROLS = 30


def require(p: Path, label: str) -> Path:
    if not p.exists():
        raise SystemExit(f"entree introuvable ({label}) : {p}")
    return p


def read_sto() -> dict[str, str]:
    seqs: dict[str, str] = {}
    for line in require(STO, "alignement Pfam").read_text().splitlines():
        if not line or line[0] == "#" or line.startswith("//"):
            continue
        n, s = line.split()
        seqs[n] = seqs.get(n, "") + s
    return seqs


def load_meta() -> tuple[list[dict], dict[str, list[dict]]]:
    rows = list(csv.DictReader(require(META, "table des feuilles").open(), delimiter="\t"))
    byacc: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        byacc[r["accession"]].append(r)
    return rows, byacc


def union_complete(segs: list[dict]) -> bool:
    return all(any(s["motif"][i] != "-" for s in segs) for i in range(3))


def fetch(acc: str) -> str | None:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{acc}.fasta"
    if path.exists() and path.stat().st_size > 0:
        return path.read_text()
    for attempt in range(3):
        try:
            with urllib.request.urlopen(
                f"https://rest.uniprot.org/uniprotkb/{acc}.fasta", timeout=30
            ) as fh:
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


def do_prep() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    seqs = read_sto()
    ncol = len(seqs[ANCHOR])
    for label, (c, aa) in COLS.items():
        got = seqs[ANCHOR][c].upper()
        if got != aa:
            raise SystemExit(f"ANCRAGE FAUX : col {c} de l'ancre = {got!r}, attendu {aa!r} ({label})")

    with (OUTDIR / "ref471.fasta").open("w") as fh:
        for i, (n, s) in enumerate(sorted(seqs.items()), start=1):
            if len(s) != ncol:
                raise SystemExit("alignement de reference non rectangulaire")
            fh.write(f">ref{i:04d}\n{s.replace('.', '-').upper()}\n")
    (OUTDIR / "ref_names.tsv").write_text(
        "ref_name\tpfam_name\n"
        + "".join(f"ref{i:04d}\t{n}\n" for i, n in enumerate(sorted(seqs), start=1))
    )

    _, byacc = load_meta()
    targets = sorted(a for a, v in byacc.items() if not union_complete(v))
    controls = sorted(
        a for a, v in byacc.items()
        if len(v) == 1 and v[0]["motif"] == "ECH" and v[0]["domaine"] == "Eukaryota"
    )[:N_CONTROLS]

    kept, missing = [], []
    with (OUTDIR / "add.fasta").open("w") as fh:
        for kind, accs in (("cible", targets), ("controle", controls)):
            for acc in accs:
                fa = fetch(acc)
                if fa is None:
                    missing.append((kind, acc))
                    continue
                s = seq_of(fa)
                fh.write(f">{kind}__{acc}\n{s}\n")
                kept.append((kind, acc, len(s)))

    report = [
        f"Reference : {len(seqs)} sequences x {ncol} colonnes (ref471.fasta)",
        f"Cibles (triade incomplete apres union des segments) : {len(targets)}",
        f"Controles (eucaryotes ECH, un seul segment) : {len(controls)}",
        f"Sequences completes ajoutees : {len(kept)}",
        f"Indisponibles chez UniProt (entrees supprimees) : {len(missing)} "
        f"-> {', '.join(a for _, a in missing)}",
    ]
    (OUTDIR / "prep_report.txt").write_text("\n".join(report) + "\n")
    print("\n".join(report))
    return 0


def do_read() -> int:
    added = require(OUTDIR / "added.fasta", "retour de MAFFT")
    seqs: dict[str, str] = {}
    name = None
    for line in added.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line.strip()

    ref_len = {len(s) for n, s in seqs.items() if n.startswith("ref")}
    if len(ref_len) != 1:
        raise SystemExit(f"longueurs de reference incoherentes : {ref_len}")
    ncol = ref_len.pop()
    if ncol != 471:
        raise SystemExit(f"--keeplength n'a pas ete respecte : {ncol} colonnes au lieu de 471")

    _, byacc = load_meta()
    rows, report = [], []
    counts: Counter = Counter()
    for n, s in seqs.items():
        if n.startswith("ref"):
            continue
        kind, acc = n.split("__", 1)
        motif = "".join(s[c].upper() if s[c] not in "-." else "-" for c, _ in COLS.values())
        segs = byacc.get(acc, [])
        pfam_motifs = ";".join(x["motif"] for x in segs)
        dom = segs[0]["domaine"] if segs else "?"
        counts[(kind, dom, motif)] += 1
        rows.append(
            dict(kind=kind, accession=acc, domaine=dom, n_segments=len(segs),
                 motifs_pfam=pfam_motifs, motif_proteine_complete=motif,
                 longueur=str(len(s.replace("-", ""))))
        )

    with (OUTDIR / "verdict.tsv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(sorted(rows, key=lambda r: (r["kind"], r["domaine"], r["accession"])))

    ctrl = [r for r in rows if r["kind"] == "controle"]
    cib = [r for r in rows if r["kind"] == "cible"]
    n_ok = sum(1 for r in ctrl if r["motif_proteine_complete"] == "ECH")
    report.append("=== CONTROLE DE LA METHODE (accessions dont le hit Pfam porte ECH) ===")
    report.append(f"  ECH retrouve : {n_ok}/{len(ctrl)}")
    for m, k in Counter(r["motif_proteine_complete"] for r in ctrl).most_common():
        report.append(f"    motif {m} : {k}")
    report.append("")
    report.append("=== CIBLES : motif de la PROTEINE COMPLETE realignee au profil ===")
    for (kind, dom, motif), k in sorted(counts.items()):
        if kind == "cible":
            report.append(f"  {dom[:22]:22s} {motif:4s} n={k}")
    n_rec = sum(1 for r in cib if r["motif_proteine_complete"] == "ECH")
    report.append("")
    report.append(
        f"BILAN : {n_rec}/{len(cib)} des accessions declarees sans triade portent les trois "
        f"residus une fois la proteine entiere realignee au profil."
    )
    report.append("")
    report.append("=== DETAIL ===")
    for r in sorted(cib, key=lambda r: (r["domaine"], r["motif_proteine_complete"], r["accession"])):
        report.append(
            f"  {r['domaine'][:10]:10s} {r['accession']:12s} len={r['longueur']:>5s} "
            f"pfam={r['motifs_pfam']:12s} complete={r['motif_proteine_complete']}"
        )
    txt = "\n".join(report) + "\n"
    (OUTDIR / "report.txt").write_text(txt)
    print(txt)
    return 0


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "prep"
    sys.exit(do_prep() if mode == "prep" else do_read())
