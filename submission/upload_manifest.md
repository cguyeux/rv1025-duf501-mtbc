# Manifeste de téléversement — Rv1025

Préparé le 2026-08-25. Ordre de dépôt et libellé à saisir pour chaque fichier.

Revue cible : **à confirmer**. Archives of Microbiology a été écartée le 2026-08-25 :
le brouillon Snapp `530bf3c9-e9b3-4911-94f5-77b916fec16e` ouvert sur ce titre porte le
manuscrit Rv3222c, et envoyer deux dossiers du même auteur à la même revue reproduirait
la concentration qui a valu deux refus de bureau au Journal of Bacteriology en août.

Le portail demande, pour chaque fichier, un type et une légende. Les voici, prêts
à coller, pour ne pas les improviser à quatorze reprises.

## 1. Manuscrit principal

| Fichier | Type | Légende |
|---|---|---|
| `main.pdf` | Manuscript | A conserved cysteine–histidine–glutamate metal site identifies DUF501 (Rv1025), an essential uncharacterised protein family of *Mycobacterium tuberculosis*, as a candidate metalloenzyme and drug target |
| `submission/rv1025_source.zip` | LaTeX source | main.tex, references.bib, main.bbl et les trois figures, 2,8 Mo — à déposer si le portail demande les sources compilables plutôt qu'un PDF |

17 pages, figures incluses, recompilée le 2026-08-25 après ajout des quatre renvois
manquants aux tables supplémentaires. Zéro erreur, zéro citation non résolue.

## 2. Figures supplémentaires

Déposer les versions PDF, vectorielles. Les PNG du même répertoire sont des copies
de travail et n'ont pas à être téléversées.

| Fichier | Type | Légende |
|---|---|---|
| `supplementary_materials/figureS_conservation.pdf` | Supplementary Figure S1 | Conservation of Cys113, His115 and Glu59 across the DUF501 family, independent of alignment redundancy |
| `supplementary_materials/figureS_druggability.pdf` | Supplementary Figure S2 | Pocket detection on the AlphaFold model: druggability score and cavity volume |
| `supplementary_materials/figureS_synteny_phylo.pdf` | Supplementary Figure S3 | Schematic of the eno-divIC-Rv1025-ppx2 block along the Actinobacterial taxonomy |

## 3. Tables supplémentaires citées dans le manuscrit

| Fichier | Type | Légende |
|---|---|---|
| `supplementary_materials/table_S1_synteny.csv` | Supplementary Table S1 | Synteny of the eno-divIC-Rv1025-ppx2 block across 63 genomes |
| `supplementary_materials/table_S2_orthogonal_metal.tsv` | Supplementary Table S2 | Orthogonal predictions of the metal-binding site |
| `supplementary_materials/table_S3_druggability.tsv` | Supplementary Table S3 | Druggability metrics of the detected pockets |
| `supplementary_materials/table_S4_synteny_null.tsv` | Supplementary Table S4 | Null model for the synteny signal |
| `supplementary_materials/table_S8_dca.tsv` | Supplementary Table S8 | Direct coupling analysis: top-ranked coevolving pairs |
| `supplementary_materials/table_S9_duf501_family.tsv` | Supplementary Table S9 | Conservation of the three metal-site positions across 2 228 DUF501 sequences |

## 4. Tables supplémentaires S5, S6, S7 et S10 : point résolu

Ces quatre tables n'étaient citées nulle part dans le manuscrit, alors que les
chiffres qu'elles portent y figuraient déjà en toutes lettres. Les renvois ont été
ajoutés le 2026-08-25, en anglais et en français, aux endroits où ces chiffres
apparaissent :

| Table | Contenu | Où elle est désormais appelée |
|---|---|---|
| S5 | conservation des résidus de poche contre tirage aléatoire apparié en taille | après « above 99.5% of size-matched random residue sets (p=0.005) » |
| S6 | les deux clusters de résidus quasi invariants, enfoui et exposé | après « a three-residue N-terminal patch » |
| S7 | même mesure sous pondération Henikoff, correction de redondance | après « 67% against a 51% background (p=0.011) » |
| S10 | les 102 membres eucaryotes de la famille | après « We report this eukaryotic branch as an observation of sequence » |

Les dix tables et les trois figures supplémentaires sont donc toutes appelées, la
numérotation ne saute plus, et le contrôle technique de l'éditeur ne trouvera pas
de fichier orphelin. Les deux versions recompilent proprement, 17 pages en anglais
et 19 en français.

Le dépôt comprend donc l'ensemble des treize fichiers supplémentaires.

## 5. Autres documents

| Fichier | Type | Remarque |
|---|---|---|
| `submission/cover_letter.txt` | Cover letter | 429 mots. ATTENTION : elle est adressée à Archives of Microbiology et cite la revue deux fois. À réadresser à la revue retenue avant dépôt |
| `main_fr.tex` et `main_fr.pdf` | ne pas déposer | Version française, pour l'archive du projet, sans destinataire chez l'éditeur |
