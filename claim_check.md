# Registre de vérification des claims

**Article :** A conserved cysteine–histidine–glutamate metal site identifies DUF501 (Rv1025)… as a candidate metalloenzyme and drug target
**Dernière vérification :** 2026-11-02 (relecture de fond AA.7 du manuscrit recadré, section ci-dessous ; le tableau du 2026-07-05 porte sur l'ancienne version)
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


## Relecture de fond AA.7 du manuscrit recadré (2026-11-02)

Portée : les treize affirmations du résumé et la conclusion entière, recoupées avec les fichiers de résultats du projet (`résultats/`), l'audit AA.5 et l'API InterPro interrogée le jour même. Le tableau du 2026-07-05 ci-dessus porte sur la version avant recadrage et reste valide pour ses lignes.

| # | Section | Claim | Statut | Preuve |
|---|---------|-------|--------|--------|
| 14 | Résumé | 4 370 protéines, 0 GO, 0 structure résolue, 0 parent de repli | confirmé | `conservation_profile.tsv`, InterPro 2026-11-02 : 4 375 protéines et 0 structure pour PF04417 (le décompte dérive de quelques unités entre versions) |
| 15 | Résumé | essentiel dans trois genres actinobactériens | confirmé | claim 3 (Goodall 2023, texte lu) |
| 16 | Résumé | 1 472 séquences bactériennes quasi complètes toutes porteuses, aucune substitution non conservative sur 2 228 | confirmé | `duf501_family_summary.txt` (0 substitution, 86 gaps) ; `AA5_audit_perte_eucaryote.md` (1 472/1 472 au niveau protéine). Précision ajoutée : 2 228 est un nombre de segments (2 204 protéines) et couvre environ la moitié des 4 370 membres InterPro |
| 17 | Résumé | trois chaînes latérales convergent à distance de coordination ; AF3 holo place un métal ; prédicteur de géométrie du squelette retrouve le site ; la mutation le déplace | confirmé | claim 7 et 8 ; `metal_orthogonal/` (BioMetAll) ; `af3_mutant_out/` (déplacement de 16 à 22 Å) |
| 18 | Résumé | bloc de gènes synténique à l'échelle du phylum « au-dessus d'une base aléatoire » | **corrigé** | `synteny_null/null_summary.json` : médiane 0,55 contre 1,00, mais p empirique 0,2045 (8 blocs aléatoires sur 43 aussi synténiques). Reformulé « dans le cinquième supérieur des blocs aléatoires » (résumé, introduction) ; « reflète une sélection » devenu « est compatible avec » (discussion) |
| 19 | Résumé | seconde cavité conservée indépendamment | confirmé | `dca_interpretation.txt` : 80 paires inter-poches, p = 0,23 ; conservation 74 % contre fond 54 %, p = 0,005 |
| 20 | Résumé | la coévolution ne retrouve que le repli | confirmé | `dca_summary.txt` : 58 % de contacts (×20,6) ; 16 paires éloignées sur 155 contre 112 attendues, p = 6,0 × 10⁻⁵⁹ |
| 21 | Résumé | trois eucaryotes ont perdu le site, paire Cys-His perdue, glutamate gardé | confirmé | `AA5_audit_perte_eucaryote.md` (32 pertes candidates testées sur séquence entière, 3 confirmées, Fisher 2,7 × 10⁻⁴) |
| 22 | Résumé | première hypothèse fonctionnelle de la famille | confirmé | claim d'antériorité de P6 : PubMed entier, 0 hit pour Rv1025 |
| 23 | Résumé | calculatoire, soutient sans prouver l'occupation | confirmé | formulation conservée ; le texte dit que la confiance par ion ne prouve pas l'occupation (relocalisation du mutant) |
| 24 | Résultats §3.3 | 3 830 sur 4 378 (87 %) Actinomycetota | **corrigé** | 4 378 contre 4 370 ailleurs dans le texte sans explication : relevés à deux dates. Le texte le dit désormais ; le reste ne se décrivait plus comme tel (3 830 + 264 + 179 + 24 = 4 297, il manque 81) et passe à « suivies de » |
| 25 | Résultats §3.7 | 576 séquences sans lignée, toutes porteuses de la triade | **corrigé** | recalcul sur `duf501_family_triad.tsv` et `duf501_taxonomy.tsv` : 585 segments sans lignée, dont 576 porteurs et 9 non porteurs (fragments de 37 à 99 résidus alignés, deux au nom d'entrée d'apicomplexe) |
| 26 | Résultats §3.7 | 75 des 86 exceptions sont des fragments (76 contre 124 résidus, p = 10⁻⁴⁶) | confirmé | recalcul : 75/86 sous le seuil de 111,6 résidus, Mann-Whitney p = 2,6 × 10⁻⁴⁶ |
| 27 | Résultats §3.5 | accessibilité 19,8 % ouvert / 14,3 % vers le solvant ; témoins 10,1-13,9 % et 0,5-4,9 % ; Zif268 0,0 % | confirmé | `metal_accessibility.txt` |
| 28 | Résultats §3.5 | poche druggable 361 Å3, score 0,65 ; cavité du métal 255 Å3, 22 sphères, score 0,24 | confirmé | `druggability/pocket_metalsite_map.tsv` |
| 29 | Résultats §3.5 | poche à 74 % (67 % pondéré, p = 0,011) ; patch 92 % pondéré | confirmé | `conservation_weighted.tsv` (74,1 ; 66,5 ; 0,0111 ; 91,8) |
| 30 | Conclusion | synténie conservée sur un phylum, pas de complexe Rv1025-DivIC, site strictement conservé chez toute séquence bactérienne quasi complète, convergence holo et prédicteur sans évolution | confirmé après la correction 18 | claims 5, 6, 16, 17 |

Non recoupé à la source dans cette relecture : l'effectif efficace de 2 599 séquences pondérées par Henikoff (le script `phase25_conservation_weighted.py` l'affiche mais aucune sortie conservée ne le porte ; à régénérer si une revue le demande).

Quatre défauts corrigés (18, 24, 25 et la nature du nombre 2 228), dans `main.tex` et `main_fr.tex`. Les fichiers `submission/abstract_plain.txt`, `title_plain.txt` et `rv1025_source.zip` datent de la soumission du 2026-08-25 et sont à régénérer au dépôt effectif.

## AA.5 : arbre de la famille et cartographie du site (2026-11-02, soir)

| # | Section | Claim | Statut | Preuve |
|---|---------|-------|--------|--------|
| 31 | Méthodes, §3.7 | 363 représentants à 70 % d'identité, 129 colonnes, 4 281 branches pour le jeu complet, 0,18 site par branche | confirmé | `résultats/phase43_reduced/report.txt` (363 ; 129 ; ratio 0,18 ; 723 branches pour le jeu réduit) ; 4 281 branches : `AA5_audit_perte_eucaryote.md` |
| 32 | Méthodes | modèle Q.pfam+I+R7 par BIC parmi LG, WAG, JTT, VT, Q.pfam ; 1 000 ufboot et 1 000 SH-aLRT | confirmé | `duf501_red.iqtree` (Best-fit model BIC) ; ligne de commande `-m MFP --mset LG,WAG,JTT,VT,Q.pfam -B 1000 --alrt 1000` relevée sur `mp` |
| 33 | §3.7 | trois pertes dans un clade de cinq feuilles, avec deux fragments d'état inconnu, ufboot 62 %, SH-aLRT 0 % | confirmé | `phase44_triad_on_tree.py` et recalcul du MRCA ; étiquette « 0/62 » du nœud, ordre SH-aLRT / ufboot lu dans `duf501_red.iqtree` |
| 34 | §3.7 | minimum de 1 perte (fragments inconnus) à 2 (fragments complets) ; 40 changements par segment, 22 par protéine incomplète | confirmé | `résultats/phase44_tree_states/duf501_red_report.txt` et bornes recalculées |
| 35 | Méthodes | seuil de quasi-complétude 112 résidus (90 % de la médiane 124) | confirmé | même seuil que `phase29` (111,6) ; médiane lue sur l'alignement |
| 36 | Références | Minh 2020, Kalyaanamoorthy 2017, Hoang 2018, Guindon 2010, Minh 2021, Fitch 1971 | confirmé | métadonnées lues sur Crossref le 2026-11-02 (titres, volumes, pages, DOI) |

Limite assumée : le jeu réduit garde d'office tous les eucaryotes, toutes les séquences sans triade complète et l'ancre, il n'est donc pas un échantillon de la famille ; le nombre d'événements décrit cet échantillon, il n'estime pas un taux de perte. L'arbre complet (2 142 séquences) n'a pas été inféré (sous-déterminé, 129 sites).
| 37 | Résumé | « in one or two events on a family tree » (trois pertes eucaryotes) | confirmé | lignes 33 et 34 ci-dessus ; résumé de 246 mots, borne 250 |
| 38 | Discussion | 0,18 site par branche, clade de cinq feuilles, ufboot 62 %, SH-aLRT 0 %, pas de datation ni de test de sélection relâchée | confirmé | lignes 31 et 33 ci-dessus |
