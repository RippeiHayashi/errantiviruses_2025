# Host-taxonomic concordance permutation analyses for Lepidoptera and Annelida POL RT clades

This README accompanies two lightweight Python scripts:

```text
lepidoptera_permutation_analysis.py
annelida_permutation_analysis.py
```

The scripts test whether highly supported small clades in fixed errantivirus POL RT trees are enriched for host-taxonomic concordance. They are intended as simple, readable analyses rather than full gene-tree/species-tree reconciliation tests.

No third-party Python packages are required.

---

## Overview

Both scripts:

1. parse an IQ-TREE-style Newick tree;
2. identify internal nodes with bootstrap support at or above a selected threshold;
3. retain nodes containing a selected number of descendant tips;
4. score host-taxonomic coherence in the observed tree;
5. keep the POL RT topology fixed;
6. randomly shuffle host-taxonomic labels across tips for a specified number of permutations;
7. rescore the same supported nodes after each shuffle;
8. compare observed values with the random-label null distribution.

The default settings are:

```text
bootstrap support >= 95
3 <= number of descendant tips <= 10
10,000 permutations
random seed = 123
```

The random seed makes the Monte Carlo label shuffles reproducible.

---

# Part 1: Lepidoptera analysis

## Purpose

`lepidoptera_permutation_analysis.py` asks whether family and superfamily labels are distributed non-randomly across a fixed Lepidoptera POL RT tree.

The analysis is suitable for statements such as:

> Highly supported small POL RT clades are enriched for elements from related Lepidoptera host groups, although strict family-level monophyly is rare.

## Input tree

Example input:

```text
Tip_Lepidoptera_subclade.txt
```

Expected tip-label format:

```text
Lepidoptera_Family_SpeciesCode_errantivirus_Number
```

Example:

```text
Lepidoptera_Erebidae_Emi_errantivirus_17
```

Bootstrap values should be stored as IQ-TREE-style square-bracket comments after branch lengths, for example:

```text
):0.032[100]
```

## Taxonomic mapping

The family-to-superfamily mapping is included in the script as:

```python
FAMILY_TO_SUPERFAMILY
```

Edit this dictionary if the taxonomic scheme changes.

## Run command

Basic run:

```bash
python lepidoptera_permutation_analysis.py Tip_Lepidoptera_subclade.txt
```

Recommended reproducible run:

```bash
python lepidoptera_permutation_analysis.py Tip_Lepidoptera_subclade.txt \
  --nperm 10000 \
  --seed 123 \
  --observed-out observed_supported_clades_lepidoptera.tsv \
  --summary-out permutation_summary_lepidoptera.tsv
```

## Metrics

For each supported 3-10-tip clade, the script records:

- number of descendant tips;
- family composition;
- superfamily composition;
- whether all tips are from the same family;
- whether all tips are from the same superfamily;
- whether at least 80% of tips are from one superfamily.

## Permutation null model

The Lepidoptera script shuffles **family labels** across terminal tips while preserving the exact observed number of tips per family. Superfamily labels are then derived from the permuted family labels.

The null model asks:

> If family labels were randomly distributed across this fixed POL RT topology, how many highly supported small clades would appear family- or superfamily-coherent by chance?

---

# Part 2: Annelida analysis

## Purpose

`annelida_permutation_analysis.py` asks whether host-species and family labels are distributed non-randomly across a fixed Annelida POL RT tree.

The analysis is suitable for statements such as:

> Highly supported small POL RT clades within an Annelida-associated radiation are enriched for elements from the same host species or family relative to randomized host labels.

## Input tree

Example input:

```text
Errantivirus_Annelida_1.txt
```

Expected tip-label formats:

```text
Annelida_Polychaeta_SpeciesCode_errantivirus_Number
Annelida_Sipuncula_SpeciesCode_errantivirus_Number
```

Examples:

```text
Annelida_Polychaeta_Pdum_errantivirus_1
Annelida_Sipuncula_Snud_errantivirus_1
```

## Host mapping

The Annelida script includes the following host-code mapping:

| Tree code | Host species | Host family |
|---|---|---|
| `Apac` | `Amphiduros_pacificus` | Hesionidae |
| `Avir` | `Alitta_virens` | Nereididae |
| `Himp` | `Harmothoe_impar` | Polynoidae |
| `Lcla` | `Lepidonotus_clava` | Polynoidae |
| `Pdum` | `Platynereis_dumerilii` | Nereididae |
| `Slim` | `Sthenelais_limicola` | Sigalionidae |
| `Snud` | `Sipunculus_nudus` | Sipunculidae |

The mapping is stored in:

```python
HOST_CODE_TO_TAXON
```

Add entries to this dictionary if additional Annelida host codes are included in a future tree.

## Run command

Basic run:

```bash
python annelida_permutation_analysis.py Errantivirus_Annelida_1.txt
```

Recommended reproducible run:

```bash
python annelida_permutation_analysis.py Errantivirus_Annelida_1.txt \
  --nperm 10000 \
  --seed 123
```

## Metrics

For each supported 3-10-tip clade, the Annelida script records:

- number of descendant tips;
- host-code composition;
- species composition;
- family composition;
- whether all tips are from the same host species;
- whether all tips are from the same host family;
- whether at least 80% of tips are from one host species;
- whether at least 80% of tips are from one host family.

## Permutation null model

The Annelida script shuffles **host-species codes** across terminal tips while preserving the exact observed number of tips per host species. Family labels are then derived from the shuffled host codes.

This approach preserves both:

```text
copy number per host species
copy number per host family
```

The null model asks:

> If host-species labels were randomly distributed across this fixed POL RT topology, how many highly supported small clades would appear species- or family-coherent by chance?

## Verified output for the supplied Annelida tree

Using:

```bash
python annelida_permutation_analysis.py Errantivirus_Annelida_1.txt \
  --nperm 10000 \
  --seed 123
```

the supplied 36-tip Annelida tree produces:

| Metric | Observed | Null mean | Empirical P |
|---|---:|---:|---:|
| Same host species | 4 / 6 | 0.1282 | 0.00009999 |
| Same host family | 6 / 6 | 0.2654 | 0.00009999 |
| At least 80% one host species | 5 / 6 | 0.1608 | 0.00009999 |
| At least 80% one host family | 6 / 6 | 0.3501 | 0.00009999 |

The empirical P value uses the standard +1 correction:

```text
(1 + number of permutations >= observed) / (number of permutations + 1)
```

With 10,000 permutations, the smallest reportable value is therefore approximately `1 x 10^-4`.

---

# Command-line options

Both scripts support:

```text
--bootstrap 95       Minimum bootstrap support for scored clades
--min-tips 3         Minimum descendant tips per scored clade
--max-tips 10        Maximum descendant tips per scored clade
--nperm 10000        Number of random label permutations
--seed 123           Random seed for reproducibility
--observed-out FILE  Output table listing observed supported clades
--summary-out FILE   Output table summarizing the permutation test
```

---

# Output files

## Observed-clade table

The observed-clade table lists each supported 3-10-tip clade in the real tree, including bootstrap support, taxonomic composition, concordance metrics, and descendant tip labels.

Suggested use:

- supplementary table;
- figure annotation;
- manual inspection of representative clades.

## Permutation-summary table

The permutation-summary table includes:

```text
metric
observed
n_scored_clades
null_mean
null_median
null_2.5pct
null_97.5pct
empirical_p_ge_observed
n_permutations
random_seed
```

---

# Suggested Methods wording

> We scored internal nodes with bootstrap support >=95 and 3-10 descendant tips. To test whether supported small POL RT clades were enriched for related host taxa, we kept the tree topology fixed and randomly permuted host-taxonomic labels across terminal tips 10,000 times while preserving the observed label counts. For each permutation, we rescored the same set of supported nodes. In the Lepidoptera analysis, family labels were permuted and superfamily assignments were derived from the permuted labels. In the Annelida analysis, host-species labels were permuted and host-family assignments were derived from the permuted labels.

# Interpretation and limitations

These analyses are intentionally lightweight.

They do not:

- estimate duplication, loss, or horizontal-transfer events;
- identify a unique transmission mechanism;
- test strict host-element co-speciation;
- replace formal gene-tree/species-tree reconciliation.

They ask a narrower question:

> Are host-taxonomic labels more clustered on the observed POL RT topology than expected under random label placement?

This makes the analyses useful for illustrating non-random host-element concordance alongside fine-scale discordance and lineage-specific amplification.
