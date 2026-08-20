# Registre de vérification des claims

**Article :** A conserved cysteine–histidine–glutamate metal site identifies DUF501 (Rv1025)… as a candidate metalloenzyme and drug target
**Dernière vérification :** 2026-07-05
**Méthode dominante :** re-dérivation des chiffres depuis les fichiers de résultats du projet (les claims P1 sont nos propres résultats computationnels) + sources primaires pour les faits de littérature (voir /bib-check, 16 réfs vérifiées CrossRef).

## Claims

| # | Prio | Section | Claim (abrégé) | Statut | Vérifié le | Source / Preuve | Notes |
|---|------|---------|----------------|--------|------------|-----------------|-------|
| 1 | P1 | Results | Rv1025 essentiel (TnSeq) + CRISPRi VI −5,7 (IC −6,2..−5,2) | confirmé | 2026-07-05 | data/Rv1025_atlas_fiche.json (vi −5.717, ci −6.214/−5.202) ; DeJesus 2017, Bosch 2021 | arrondi correct |
| 2 | P1 | Results | DUF501 = 4370 protéines / 5484 taxons, 0 GO, 0 PDB | confirmé | 2026-07-05 | InterPro API PF04417/IPR007511 (vérif P1) | — |
| 3 | P1 | Results | Gène DUF501 essentiel dans 3 genres (Rv1025, cgp_1113, DIP0919) | confirmé | 2026-07-05 | Goodall 2023 (texte lu mot pour mot, P1) | — |
| 4 | P1 | Results | Repli nouveau : Foldseek complet = 0 homologue significatif | confirmé | 2026-07-05 | résultats/structure/foldseek_out (swissprot 0 ; pdb100 1 hit 24%/26%cov ; cath50 5 @13-19%, e≥0.01) | re-parsé foldseek_parse.py |
| 5 | P1 | Results | Opéron conservé 61/63 génomes, divIC voisin invariant, absent E.coli/B.subtilis, M. leprae gardé | confirmé | 2026-07-05 | résultats/synteny_operon.tsv (61 présents/intacts ; 48/48 Myco ; leprae+lepromatosis 4/4 ; absents = Bsubtilis, Ecoli) | re-compté |
| 6 | P1 | Results | Pas de complexe Rv1025-divIC : ipTM 0,25→0,08 / PAE 14,7 ; contrôle+ divIC-ftsQ 0,38 / 5,1 ; eno 0,19 ; ppx2 0,17 | confirmé | 2026-07-05 | résultats/af3_out/*summary_confidences* (re-dérivé) | valeurs exactes |
| 7 | P1 | Results | Conservation 8700 homologues ; 23 quasi-invariants ; C113/H115 100 %, E59 99 % ; donneurs convergents ~2,3 Å | confirmé | 2026-07-05 | analyses/phase5_conservation_site.py (8700 séq, 23 invariants, convergence 2,30 Å) | « family »→« homologous » corrigé (imprécision) |
| 8 | P1 | Results | Holo AF3 : Zn/Fe/Mn coordonnés 2,25-2,47 Å, pLDDT ion 94-98, Fe le plus confiant, coordination incomplète | confirmé | 2026-07-05 | analyses/phase7_metal_parse.py + résultats/af3_metal_out (Table 1 re-dérivée exacte) | — |
| 9 | P2 | Results | 5 sites polymorphes / 145 209 génomes, aucun disruptif | confirmé | 2026-07-05 | fiche atlas (snp_sites 5, n 145209, disrupt 0) | — |
| 10 | P2 | Results | pLDDT moyen du modèle AF = 95,3 | confirmé | 2026-07-05 | fiche atlas (plddt_af 95.3) ; AF-P96375-F1 (CA mean 95.3) | — |
| 11 | P3 | Intro | ppx2 (Rv1026) = exopolyphosphatase (poly(P), dormance/paroi) | confirmé | 2026-07-05 | Chuang 2015, Choi 2012, Tiwari 2019 (réfs vérifiées) | — |
| 12 | P3 | Discussion | PPK2 = cible antitub. validée de la voie poly(P) | confirmé | 2026-07-05 | Singh 2016 (réf vérifiée) ; formulation corrigée (PPK2 hors opéron, même voie que Ppx2) | corrigé en bib-check |
| 13 | P4 | Intro/Méth | Rv1025 155 aa, PF04417, IPR007511, COG1507, K09009, UniProt P96375 | confirmé | 2026-07-05 | fiche atlas + InterPro/KEGG/UniProt (P1) | — |

## Audit de cohérence numérique interne

- 4370 / 5484 : cohérent abstract ↔ results.
- 61/63, 48/48 : cohérent abstract ↔ results ↔ fichier.
- 8700 : cohérent abstract ↔ results (désormais « homologous sequences »).
- Table 1 métal (2,25-2,47 Å ; pLDDT 94-98) : cohérent avec le texte (« 2.25--2.47 », « 94--98 »).
- ipTM 0,25 / 0,38 : cohérent abstract ↔ results.
- Aucune incohérence interne détectée ; aucune somme à contrôler (pas de décomposition additive).

## Résumé

| Prio | Total | Confirmé | À corriger | Partiel | Non vérifiable |
|------|-------|----------|------------|---------|----------------|
| P1 | 8 | 8 | 0 | 0 | 0 |
| P2 | 2 | 2 | 0 | 0 | 0 |
| P3 | 2 | 2 | 0 | 0 | 0 |
| P4 | 1 | 1 | 0 | 0 | 0 |

**0 claim infirmé.** 2 imprécisions mineures corrigées dans le manuscrit lors du pipeline qualité (« family sequences »→« homologous sequences » ; formulation PPK2/voie polyphosphate). Tous les chiffres du manuscrit re-dérivent exactement depuis les fichiers de résultats du projet.
