# Code and intermediate results for the Rv1025 (DUF501) manuscript

This directory mirrors the layout of the project from which the manuscript was produced, so that the scripts run unchanged with `code/` as their root.

- `analyses/` : the analysis scripts, numbered by phase. Synteny: `phase1`, `phase10`, `phase17` (random-block null), `phase18`. AlphaFold-Multimer complex tests: `phase2`, `phase3`, `phase11`. Structural search: `phase4`, `foldseek_parse.py`. Conservation of the site: `phase5`, `phase12`, `phase16`, `phase22` to `phase25`. Holo metal and mutant control: `phase6`, `phase7`, `phase8`, `phase9`, `metal_geometry.py`. Orthogonal metal-site prediction, pockets and accessibility: `phase13`, `phase14`, `phase15`, `phase31`, `phase32`. Homodimer and divisome panel: `phase19` to `phase21`, `phase33`. Coevolution: `phase26`, `phase27`. Family-wide conservation of the triad and the eukaryotic branch: `phase28`, `phase29`, `phase30`, `phase34`, `phase35`, `phase38`, `phase39`, `phase41` to `phase43`, `phase46`. Family tree: `phase40`, `phase44`, `phase45`.
- `résultats/` : tabular and text outputs under 2 MB that carry the numbers quoted in the text. Large model outputs (AlphaFold3 and Boltz-2 archives, multiple-sequence-alignment caches) are not duplicated here.
- `data/` : input data under 2 MB, including the Pfam full alignment of PF04417 (release PF04417.18) used without realignment.

The numbers of the manuscript are traced to their files in `../claim_check.md`. A snapshot of the whole project (scripts, data, results including model outputs) is archived on Zenodo, DOI 10.5281/zenodo.22025970.
