# Metatranscriptomics bioinformatics scripts

This folder contains cleaned `.sh` versions of the eight scripts used for the perennial ryegrass metatranscriptomic workflow.

## Workflow order

1. `01_fastp_quality_trim.sh` — quality trimming with fastp.
2. `02_cutadapt_fastqc_multiqc.sh` — cutadapt trimming, FastQC, and optional MultiQC mode.
3. `03_merge_two_flowcells.sh` — merge final-trimmed reads across flowcells/lanes.
4. `04_remove_plant_host.sh` — remove perennial ryegrass reads with Bowtie2/SAMtools.
5. `05_remove_human_reads.sh` — remove possible human reads with Bowtie2.
6. `06_kraken2_taxonomic_classification.sh` — classify nonhuman reads with Kraken2.
7. `07_bracken_multirank_abundance.sh` — extract multirank taxonomic summaries and estimate genus abundance with Bracken.
8. `08_humann4_functional_profiling.sh` — run MetaPhlAn/HUMAnN4 functional profiling.

## Important edits made

- Added executable shebangs and `set -euo pipefail`.
- Added `mkdir -p` for output/log directories.
- Standardised file naming between cutadapt and flowcell merging.
- Made flowcells, lanes, database paths, and base directories clear user-editable variables.
- Added sanity checks for input files, database/index presence, and empty outputs.
- Added skip logic where final outputs already exist.
- Kept FastQC inside the cutadapt script as requested.
- Added a `MODE=multiqc` option to run MultiQC once after trimming is complete.


## Software and database versions

The exact software and database versions are listed in `software_versions.txt`. In brief, the workflow used the default available NeSI module versions for fastp through Bracken, and the project Conda environment for HUMAnN4/MetaPhlAn4/DIAMOND.

Key versions include fastp/1.1.0, cutadapt/5.2, FastQC/0.12.1, MultiQC/1.24.1, Bowtie2/2.5.4, SAMtools/1.23.1, Kraken2/2.1.6, Bracken/2.7, HUMAnN4 v4.0.0.alpha.1, MetaPhlAn4 v4.1.1, and DIAMOND v2.1.14.

## Suggested run pattern

```bash
mkdir -p logs

sbatch scripts/01_fastp_quality_trim.sh
sbatch --dependency=afterok:<JOBID_FASTP> scripts/02_cutadapt_fastqc_multiqc.sh
sbatch --dependency=afterok:<JOBID_CUTADAPT> --export=MODE=multiqc scripts/02_cutadapt_fastqc_multiqc.sh
sbatch --dependency=afterok:<JOBID_CUTADAPT> scripts/03_merge_two_flowcells.sh
sbatch --dependency=afterok:<JOBID_MERGE> scripts/04_remove_plant_host.sh
sbatch --dependency=afterok:<JOBID_PLANT_HOST> scripts/05_remove_human_reads.sh
sbatch --dependency=afterok:<JOBID_HUMAN_HOST> scripts/06_kraken2_taxonomic_classification.sh
sbatch --dependency=afterok:<JOBID_KRAKEN2> scripts/07_bracken_multirank_abundance.sh
sbatch scripts/08_humann4_functional_profiling.sh
```

Replace `<JOBID_...>` with the SLURM job IDs returned by `sbatch`.

## Notes before GitHub release

Please check the account names, module names, database paths, and the base project path before committing. These scripts are cleaned for reproducibility but still contain project-specific NeSI paths that should either be retained as a record of the original run or moved into a separate config file for public reuse.

Software and database versions have been added in `software_versions.txt`. Please still confirm reference genome/database file paths before public release, especially the Lolium perenne reference genome and GRCh38 index paths.
