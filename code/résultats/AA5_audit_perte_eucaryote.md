# AA.5 — Audit de la « perte polarisée du site chez les eucaryotes »

Séance du 2026-11-02, piste AA.5. Scripts : `analyses/phase40_tree_prep.py`,
`phase41_fulllength_triad.py`, `phase42_mafft_keeplength.py`, `phase43_reduce_alignment.py`.
Données : `data/PF04417_full.sto` (PF04417.18, 2 228 segments × 471 colonnes), UniProt REST
(sequences complètes, 2026-11-02), MAFFT 7.526 sur `mp`.

## Ce que le manuscrit affirme aujourd'hui

`article/main.tex` §3.7 (et le résumé, et l'annonce de l'introduction) :

> Eleven of them have lost part of the site, almost always the Cys--His pair with retention of
> the glutamate ($p=5\times10^{-14}$ against bacteria), and this is not explained by divergence,
> since those eleven are no less similar to Rv1025 than the ninety-one that retain the complete
> triad ($p=0.36$).

Résumé : « Where eukaryotic members have lost the site, the loss is polarised: the
cysteine--histidine pair goes, the glutamate stays. »

## Le défaut : le décompte porte sur le SEGMENT aligné, non sur la PROTÉINE

Les 11 « pertes » sont des **segments** de l'alignement Pfam dont une ou deux positions de la
triade sont des **gaps** (aucune substitution : les motifs observés sont `E--`, `-CH`, `-C-`,
`EC-`). Or 24 accessions UniProt de l'alignement y figurent par **deux segments distincts**, et
pour 23 d'entre elles ces deux segments couvrent des colonnes **complémentaires** du profil
(indice de Jaccard des colonnes occupées : 0,00 à 0,12 ; un seul cas intermédiaire à 0,18) :
c'est le même domaine détecté en deux morceaux par le HMM, la moitié N portant Glu59 et la
moitié C portant Cys113-His115. Sur plusieurs de ces accessions les bornes des deux segments se
**chevauchent** sur la protéine (de 5 à 22 aa), ce qui exclut deux copies en tandem.

**Huit des onze pertes eucaryotes sont de ce type** : l'autre segment de la même protéine porte
les résidus déclarés manquants.

| accession | motif du segment compté | motifs des deux segments | protéine |
|---|---|---|---|
| A0A2A9M7T6 | `E--` | `E--` + `-CH` | triade complète |
| A0A448ZKR0 | `E--` | `-CH` + `E--` | triade complète |
| A0A6A5CBS0 | `E--` | `-CH` + `E--` | triade complète |
| A0A9K3KWR1 | `E--` | `E--` + `-CH` | triade complète |
| A0A9K3PP58 | `E--` | `E--` + `-CH` | triade complète |
| A0A9W5T8J3 | `-CH` | `E--` + `-CH` | triade complète |
| A0AA88KMG9 | `E--` | `E--` + `-CH` | triade complète |
| F0VAI8 | `E--` | `E--` + `-CH` | triade complète |

**L'argument décisif qui qualifie l'artefact : le même motif existe chez 12 bactéries** (et
2 entrées depuis supprimées d'UniProt), dont la triade est elle aussi reconstituée par union des
deux segments. La « dégénérescence complémentaire » n'est donc pas un trait de la branche
eucaryote : c'est une signature de la segmentation du HMM, présente dans les deux domaines du
vivant. La phase 38 avait vu le phénomène sur 3 apicomplexes et l'avait lu comme un tandem à
dégénérescence complémentaire, sans étendre le contrôle au reste de la famille.

## Contrôle : les 3 pertes restantes sont réelles

Les 32 accessions dont la triade reste incomplète après union des segments ont été testées sur
leur **séquence UniProt entière**, réalignée au profil Pfam par `mafft --add --keeplength`
(les 471 colonnes de référence sont préservées, donc les colonnes 172/395/397 restent lisibles).

- **Contrôle de la méthode** : 10 accessions dont le hit Pfam porte déjà `ECH` → `ECH` retrouvé
  **10/10**.
- **Cibles** : **0/32** ne retrouve la triade. Les protéines complètes confirment `E--` (22),
  `-CH` (5+1), `-C-` (2), `EC-` (1).

La perte est donc réelle pour ces accessions : elle ne se réduit pas à un hit HMM partiel.

Une méthode préalable par alignement **local** par paire contre Rv1025 (phase 41) a été
**écartée** : elle tronque l'alignement (spans de 30 à 50 colonnes) et son contrôle positif
n'atteint que 9/12. Ses sorties ne doivent pas être citées.

## Les chiffres à porter dans le manuscrit

Au niveau de la **protéine**, sur les accessions ayant au moins un segment quasi complet
(≥ 90 % de la longueur médiane 124 des membres à triade complète, critère de la phase 29) :

| domaine | triade complète | incomplète | taux |
|---|---|---|---|
| Bacteria | 1 472 / 1 472 | 0 | 0 % |
| Eukaryota | 99 / 102 | **3** | 2,9 % |
| supprimé d'UniProt | 526 / 526 | 0 | 0 % |
| non retourné par UniProt | 50 / 50 | 0 | 0 % |

- Fisher eucaryotes contre bactéries, niveau protéine : table `[[3, 99], [0, 1472]]`,
  **p = 2,7 × 10⁻⁴** (et non 5 × 10⁻¹⁴, qui est la valeur au niveau segment).
- L'**invariance bactérienne stricte est confirmée** : 1 472 / 1 472, aucune substitution non
  conservative. Ce résultat central du manuscrit n'est pas touché.
- La **polarité est confirmée mais sur n = 3** : les trois pertes sont `E--` (glutamate retenu,
  paire Cys-His perdue). « almost always » devient « in all three ».
- Le test de divergence (p = 0,36 sur onze contre quatre-vingt-onze) n'a plus de puissance sur
  n = 3 : à refaire ou à retirer.

## Les trois pertes confirmées

| accession | organisme | position taxonomique | motif |
|---|---|---|---|
| A0A1D3CV76 | *Cyclospora cayetanensis* | Eucoccidiorida / Eimeriorina / Eimeriidae | `E--` |
| U6GA88 | *Eimeria acervulina* | Eucoccidiorida / Eimeriorina / Eimeriidae | `E--` |
| A0A2C6L2K9 | *Cystoisospora suis* | Eucoccidiorida / Eimeriorina / Sarcocystidae | `E--` |

Les trois appartiennent au même sous-ordre d'apicomplexes (Eimeriorina), mais à deux familles
distinctes, et d'autres Eimeriorina de l'alignement (*Besnoitia besnoiti*, *Neospora caninum*)
portent la triade complète — par deux segments complémentaires, précisément. La perte n'est donc
pas ancestrale au sous-ordre : au moins deux événements indépendants, à confirmer sur l'arbre.

## Réserve de méthode, à écrire si ces chiffres sont publiés

Un gap après réalignement de la protéine entière au profil signifie « région non alignable au
profil », ce qui est le meilleur test disponible mais n'est pas une preuve d'absence du résidu
dans la protéine. Le contrôle positif à 10/10 borne ce risque sans l'annuler.
