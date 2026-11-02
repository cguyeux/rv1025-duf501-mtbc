#!/usr/bin/env python3
"""Phase 40 (piste AA.5) - preparation de l'alignement PF04417 pour l'arbre ML.

Objet       : filtrer l'alignement Pfam full de DUF501 et produire la table des feuilles
              (etat des trois positions de la triade, domaine du vivant, lignee) qui sera
              cartographiee sur l'arbre ML.
Entrees     : data/PF04417_full.sto (alignement Pfam full, 2228 x 471)
              résultats/duf501_family_triad.tsv, duf501_taxonomy.tsv,
              duf501_eukaryote_uniprot.tsv (phases 28-29, pour controle et lignees)
Sorties     : résultats/phase40_tree/align_full.fasta   (toutes les sequences, colonnes filtrees)
              résultats/phase40_tree/align_tree.fasta   (sequences retenues pour l'arbre)
              résultats/phase40_tree/leaf_metadata.tsv  (etat de la triade, domaine, lignee)
              résultats/phase40_tree/name_map.tsv       (nom sur Newick <-> nom Pfam)
              résultats/phase40_tree/prep_report.txt
Reutilisable: oui pour toute famille Pfam dont on veut cartographier un site sur un arbre
              (changer STO, ANCHOR et TRIAD_COLS_1B).
Projet      : mtbc/clos_soumis/Rv1025, piste AA.5
Date        : 2026-11-02

Choix documentes :
 - colonnes gardees : occupation (fraction de residus non-gap) >= COL_OCC ;
   les trois colonnes de la triade sont verifiees presentes apres filtrage.
 - sequences gardees pour l'arbre : occupation des colonnes retenues >= SEQ_OCC.
 - renommage : les noms Pfam portent '/' et des bornes de domaine, illegaux en Newick ;
   chaque sequence recoit un identifiant court, la correspondance est ecrite.
"""
from __future__ import annotations

import sys
from collections import Counter
from pathlib import Path

PROJ = Path(__file__).resolve().parent.parent
STO = PROJ / "data" / "PF04417_full.sto"
OUTDIR = PROJ / "résultats" / "phase40_tree"

COL_OCC = 0.50   # occupation minimale d'une colonne
SEQ_OCC = 0.70   # occupation minimale d'une sequence sur les colonnes retenues

# Ancre : P96375_MYCTU (Rv1025), domaine 16-137. Les colonnes rapportees par phase28
# (prep_report "colonne d'alignement") sont des index 0-BASED dans la ligne Stockholm.
ANCHOR = "P96375_MYCTU/16-137"
TRIAD_COLS_0B = {"E59": (172, "E"), "C113": (395, "C"), "H115": (397, "H")}


def read_stockholm(path: Path) -> dict[str, str]:
    seqs: dict[str, str] = {}
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line.startswith("//"):
                continue
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(f"ligne Stockholm inattendue : {line[:80]}")
            name, seq = parts
            if name in seqs:
                seqs[name] += seq
            else:
                seqs[name] = seq
    return seqs


def main() -> int:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    seqs = read_stockholm(STO)
    names = list(seqs)
    ncol = len(seqs[names[0]])
    if any(len(s) != ncol for s in seqs.values()):
        raise ValueError("alignement non rectangulaire")
    report: list[str] = []
    report.append(f"Alignement lu : {len(names)} sequences x {ncol} colonnes")

    if ANCHOR not in seqs:
        raise ValueError(f"ancre absente : {ANCHOR}")

    # --- occupation par colonne -------------------------------------------------
    gaps = set("-.")
    occ = []
    for j in range(ncol):
        n = sum(1 for nm in names if seqs[nm][j] not in gaps)
        occ.append(n / len(names))
    keep = [j for j in range(ncol) if occ[j] >= COL_OCC]
    report.append(f"Colonnes gardees (occupation >= {COL_OCC:.0%}) : {len(keep)}/{ncol}")
    for thr in (0.1, 0.25, 0.5, 0.75, 0.9):
        report.append(f"  colonnes a occupation >= {thr:.0%} : {sum(1 for o in occ if o >= thr)}")

    # positions de la triade dans le nouvel indexage
    triad_new = {}
    for label, (j, expected) in TRIAD_COLS_0B.items():
        anchor_res = seqs[ANCHOR][j].upper()
        if anchor_res != expected:
            raise SystemExit(
                f"ANCRAGE FAUX : colonne {j} de {ANCHOR} = {anchor_res!r}, attendu {expected!r} "
                f"pour {label}. L'alignement ou la convention d'index a change : ne pas continuer."
            )
        if j not in keep:
            raise ValueError(f"colonne de la triade {label} (col {j}) perdue au filtrage")
        triad_new[label] = keep.index(j) + 1
        report.append(
            f"  triade {label} : colonne Pfam {j} (0-based) -> colonne filtree "
            f"{triad_new[label]} ; residu de l'ancre = {anchor_res} (attendu {expected}, OK)"
        )

    # --- occupation par sequence sur les colonnes retenues ---------------------
    filt = {nm: "".join(seqs[nm][j] for j in keep).replace(".", "-").upper() for nm in names}
    seq_occ = {nm: sum(1 for c in filt[nm] if c != "-") / len(keep) for nm in names}
    tree_names = [nm for nm in names if seq_occ[nm] >= SEQ_OCC]
    report.append(
        f"Sequences retenues pour l'arbre (occupation >= {SEQ_OCC:.0%}) : "
        f"{len(tree_names)}/{len(names)} ; ecartees : {len(names) - len(tree_names)}"
    )

    # --- metadonnees : triade, domaine, lignee ---------------------------------
    triad_tsv = PROJ / "résultats" / "duf501_family_triad.tsv"
    taxo_tsv = PROJ / "résultats" / "duf501_taxonomy.tsv"
    euk_tsv = PROJ / "résultats" / "duf501_eukaryote_uniprot.tsv"

    triad_prev: dict[str, dict[str, str]] = {}
    with triad_tsv.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) != len(hdr):
                continue
            triad_prev[f[0]] = dict(zip(hdr, f))

    domain: dict[str, str] = {}
    with taxo_tsv.open() as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) >= 2:
                domain[f[0]] = f[1]

    lineage: dict[str, str] = {}
    with euk_tsv.open() as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        icol = hdr.index("Taxonomic lineage") if "Taxonomic lineage" in hdr else 4
        ocol = hdr.index("Organism") if "Organism" in hdr else 3
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(icol, ocol):
                lineage[f[0]] = f[ocol] + " | " + f[icol]

    def acc_of(name: str) -> str:
        return name.split("/")[0].split("_")[0]

    # --- renommage sur pour Newick --------------------------------------------
    name_map = {}
    width = len(str(len(names)))
    for i, nm in enumerate(sorted(names), start=1):
        name_map[nm] = f"s{i:0{width}d}"

    with (OUTDIR / "name_map.tsv").open("w") as fh:
        fh.write("newick_name\tpfam_name\n")
        for nm in sorted(names):
            fh.write(f"{name_map[nm]}\t{nm}\n")

    def write_fasta(path: Path, subset: list[str]) -> None:
        with path.open("w") as fh:
            for nm in subset:
                fh.write(f">{name_map[nm]}\n")
                s = filt[nm]
                for k in range(0, len(s), 60):
                    fh.write(s[k:k + 60] + "\n")

    write_fasta(OUTDIR / "align_full.fasta", sorted(names))
    write_fasta(OUTDIR / "align_tree.fasta", sorted(tree_names))

    # --- table des feuilles ---------------------------------------------------
    counts = Counter()
    with (OUTDIR / "leaf_metadata.tsv").open("w") as fh:
        fh.write(
            "newick_name\tpfam_name\taccession\tdomaine\tgroupe\t"
            "res59\tres113\tres115\tmotif\ttriade_stricte\tin_tree\tseq_occupancy\torganism_lineage\n"
        )
        for nm in sorted(names):
            acc = acc_of(nm)
            prev = triad_prev.get(nm, {})
            r59 = filt[nm][triad_new["E59"] - 1]
            r113 = filt[nm][triad_new["C113"] - 1]
            r115 = filt[nm][triad_new["H115"] - 1]
            motif = f"{r59}{r113}{r115}".replace("-", "-")
            strict = "oui" if (r59 == "E" and r113 == "C" and r115 == "H") else "non"
            # controle de coherence avec phase28
            if prev and prev.get("triade_stricte") and prev["triade_stricte"] != strict:
                counts["desaccord_phase28"] += 1
            dom = domain.get(acc, "non renseigne")
            counts[(dom, strict)] += 1
            fh.write(
                "\t".join(
                    [
                        name_map[nm], nm, acc, dom, prev.get("groupe", "?"),
                        r59, r113, r115, motif, strict,
                        "oui" if nm in set(tree_names) else "non",
                        f"{seq_occ[nm]:.3f}", lineage.get(acc, ""),
                    ]
                )
                + "\n"
            )

    report.append(f"Desaccords avec phase28 sur la triade stricte : {counts['desaccord_phase28']}")
    report.append("Etat de la triade par domaine (toutes sequences) :")
    doms = sorted({k[0] for k in counts if isinstance(k, tuple)})
    for d in doms:
        n_oui = counts.get((d, "oui"), 0)
        n_non = counts.get((d, "non"), 0)
        report.append(f"  {d:28s} triade stricte {n_oui:5d} / {n_oui + n_non:5d}")

    txt = "\n".join(report) + "\n"
    (OUTDIR / "prep_report.txt").write_text(txt)
    print(txt)
    return 0


if __name__ == "__main__":
    sys.exit(main())
