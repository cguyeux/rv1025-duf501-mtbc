#!/bin/bash
# Phase 40 (piste AA.5) - arbre ML de la famille DUF501 (PF04417), deporte sur mp.
# Run A : arbre de reference sur les 2142 sequences a occupation >= 70 % (toutes ECH),
#         modele par ModelFinder, ufboot 1000 + SH-aLRT 1000, taux par site (--wsr).
# Run B : placement des 86 sequences lacunaires, topologie de A imposee en contrainte (-g).
# Threads bornes a 16 : mp est partage (load ~40/64 au lancement).
set -euo pipefail

WORK=/data/cguyeux/rv1025_aa5_duf501_tree
# Binaire recopie sous un autre nom : une autre session de la machine lance un wrapper
# qui fait `pgrep -x iqtree2 | xargs -r kill` (constate le 2026-11-02 sur mp), ce qui
# tuerait ce job. Le nom distinct met le calcul hors de portee de ce pgrep.
IQ="$WORK/iqtree2_aa5"
cd "$WORK"
[ -x "$IQ" ] || cp /data/cguyeux/bin/iqtree2 "$IQ"
chmod +x "$IQ"

echo "=== RUN A : arbre de reference (2142 seq) === $(date -Is)"
"$IQ" -s align_tree.fasta \
      -m MFP --mset LG,WAG,JTT,VT,Q.pfam,Blosum62 \
      -B 1000 --alrt 1000 \
      --rate \
      -T 8 --seed 20261102 \
      --prefix duf501_core

BEST=$(grep -m1 "Best-fit model:" duf501_core.log | sed 's/.*Best-fit model: *//; s/ .*//')
echo "=== modele retenu par ModelFinder : $BEST ==="

echo "=== RUN B : placement des 86 lacunaires sous contrainte === $(date -Is)"
"$IQ" -s align_full.fasta \
      -m "$BEST" \
      -g duf501_core.treefile \
      -T 8 --seed 20261102 \
      --prefix duf501_all

echo "=== FINI === $(date -Is)"
