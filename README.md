# 6PP_Ryegrass_Metatranscriptomics

This repository contains bioinformatic workflows, statistical analyses, and figure-generation scripts used in a perennial ryegrass metatranscriptomic study investigating the effects of 6-pentyl-α-pyrone (6-PP) seed coating under drought and post-drought recovery conditions.

## Repository structure

```text
6PP_Ryegrass_Metatranscriptomics/

├── 1.bioinformatics/
├── 2.germination_endophyte/
├── 3.plant_biomass/
├── 4.plant_transcriptome/
├── 5.microbial_taxa/
├── 6.microbial_function/
├── 7.manuscript_figure/
├── 8.trait_function_linkage/
├── docs/
├── README.md
├── LICENSE
└── .gitignore
```

## Folder descriptions

### `1.bioinformatics`

Contains shell scripts for preprocessing and analysis of metatranscriptomic reads.

Scripts should be run sequentially.

```text
01_fastp.sh
02_cutadapt.sh
03_host_removal.sh
04_kraken2.sh
05_bracken.sh
06_metaphlan4.sh
07_humann4.sh
08_merge_outputs.sh
```

### `2.germination_endophyte`

Generates manuscript Table 1.

Input files:

```text
data/
```

Outputs:

```text
results/
7.manuscript_figure/1.table1/
```

Run:

```r
Table1_Germination_Endophyte.R
```

### `3.plant_biomass`

Generates manuscript Figure 1.

Run:

```r
Figure1_Plant_Biomass.R
```

Outputs:

```text
results/
7.manuscript_figure/2.figure1/
```

### `4.plant_transcriptome`

Generates manuscript Figure 2.

Run:

```r
Figure2_Plant_Transcriptome.R
```

Outputs:

```text
results/
7.manuscript_figure/3.figure2/
```

### `5.microbial_taxa`

Generates manuscript Figure 3.

Run:

```r
Figure3_Microbial_Taxa.R
```

Outputs:

```text
results/
7.manuscript_figure/4.figure3/
```

### `6.microbial_function`

Generates manuscript Figure 4.

Run:

```r
Figure4_Microbial_Function.R
```

Outputs:

```text
results/
7.manuscript_figure/5.figure4/
```

### `8.trait_function_linkage`

Generates manuscript Figure 5.

Run:

```r
Figure5_Trait_Function_Linkage.R
```

Outputs:

```text
results/
7.manuscript_figure/6.figure5/
```

## Notes

The repository structure should be preserved after cloning.

R scripts were written assuming the following hierarchy:

```text
6PP_Ryegrass_Metatranscriptomics/
├── 1.bioinformatics/
├── 2.germination_endophyte/
├── 3.plant_biomass/
├── 4.plant_transcriptome/
├── 5.microbial_taxa/
├── 6.microbial_function/
├── 7.manuscript_figure/
└── 8.trait_function_linkage/
```

No modifications to file paths should be required if this directory structure is maintained.

Software versions and database information are provided in:

```text
docs/software_versions.txt
```

## License

This repository is distributed under the MIT License.
