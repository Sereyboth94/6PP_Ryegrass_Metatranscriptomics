# 6PP Ryegrass Metatranscriptomics

Bioinformatic workflows, statistical analyses and manuscript figures for a perennial ryegrass study of **6-pentyl-α-pyrone (6-PP) seed coating during drought and post-drought recovery**. The analyses cover plant biomass and physiology, host gene expression, microbial taxonomy, microbial functions, and links between plant traits and microbial pathways.

Repository: [Sereyboth94/6PP_Ryegrass_Metatranscriptomics](https://github.com/Sereyboth94/6PP_Ryegrass_Metatranscriptomics)

Documentation updated: **13 October 2026**. This README incorporates the larger-label figure instructions.

## Repository folders and entry points

Paths below are relative to the project root.

| Folder | Purpose | Entry point |
| --- | --- | --- |
| `1.bioinformatics/` | Read preprocessing, host removal, taxonomy and functional profiling | Shell workflow listed below |
| `2.germination_endophyte/` | Germination and endophyte status; manuscript Table 1 | `Table1_Germination_Endophyte.R` in this analysis folder |
| `3.plant_biomass/` | Plant biomass and physiology; Figure 1 | `3.plant_biomass/script/Fig.1_Generated_Script.R` |
| `4.plant_transcriptome/` | Plant gene expression; Figure 2 | `4.plant_transcriptome/script/Figure2_Plant_Transcriptome.R` |
| `5.microbial_taxa/` | Microbial taxonomy; Figure 3 and genus composition Figure S4 | `5.microbial_taxa/script/Figure3_Taxonomy.R` |
| `6.microbial_function/` | Microbial functional profiles; Figure 4 | `6.microbial_function/scripts/Figure4_Function.R` |
| `7.manuscript_figure/` | Assembled manuscript figures | Output folders listed below |
| `8.trait_function_linkage/` | Plant trait–microbial pathway associations; Figure 5 | `8.trait_function_linkage/script/Figure5_Trait_Function_Linkage.R` |
| `docs/` | Software versions and database information | `docs/software_versions.txt` |
| `Run_All_Figures.R` | Runs the five figure scripts in separate R sessions | Project-root runner |
| `Update_GitHub_Repo.ps1` | Copies the selected updated files, commits them and pushes to GitHub | PowerShell update script |

Keep this folder structure when cloning the repository. Figure 4 uses **`scripts/`**, whereas the other figure entry points use **`script/`**.

## Bioinformatics workflow

Run the preprocessing and profiling steps sequentially, with the input paths, databases and compute settings configured for the analysis environment:

| Order | Script | Main task |
| --- | --- | --- |
| 1 | `01_fastp.sh` | Initial read quality control |
| 2 | `02_cutadapt.sh` | Adapter/sequence trimming |
| 3 | `03_host_removal.sh` | Remove host reads for microbial profiling |
| 4 | `04_kraken2.sh` | Taxonomic classification |
| 5 | `05_bracken.sh` | Taxonomic abundance estimation |
| 6 | `06_metaphlan4.sh` | Microbial marker profiling |
| 7 | `07_humann4.sh` | Gene-family and pathway profiling |
| 8 | `08_merge_outputs.sh` | Merge sample-level outputs |

The figure runner starts from processed analysis inputs; it does not run these shell workflows. Software and database versions belong in `docs/software_versions.txt`.

## Run the manuscript figures

Place `Run_All_Figures.R` at the project root. It resolves the five scripts from the paths in the table above.

For the current Windows project, run this in RStudio:

```r
Sys.setenv(
  MANUSCRIPT_PROJECT_ROOT =
    "C:/Users/soths/OneDrive - Lincoln University/Writing/NatMicro/github"
)
source(file.path(Sys.getenv("MANUSCRIPT_PROJECT_ROOT"), "Run_All_Figures.R"))
```

For a different checkout, set `MANUSCRIPT_PROJECT_ROOT` to that checkout's root. If the variable is unset, the runner uses the folder containing `Run_All_Figures.R`.

At the top of the runner:

```r
RUN_FIGURES <- 1:5       # Example: c(2, 3, 4) for selected figures.
EXPORT_WIDTH_IN <- 8.0
FONT_SCALE <- 1.0        # Example: 1.10 for 10% larger text.
SYMBOL_SCALE <- 0.80
DISPLAY_WIDTH_IN <- 6.27
```

The runner parses all selected scripts, checks their required files and packages, and launches each script in a separate R session. It checks the PDF, PNG and TIFF outputs for each successful run, including the genus composition supplement when Figure 3 is selected. Individual script revisions are reported rather than required to be identical.

These are complete analysis scripts, so rerunning can take time. Figure 2 uses AnnotationHub and may require its annotation cache or internet access.

### Required processed inputs

Supply these files under the corresponding analysis folders before running. The runner checks their presence; package and analysis-specific checks continue inside each figure script.

| Figure | Required input location and files |
| --- | --- |
| 1 | `3.plant_biomass/data/`: SPAD and fresh/dry-weight `.xlsx` workbooks; `Biotron_Experiment.png` |
| 2 | `4.plant_transcriptome/data/`: `featureCounts_ryegrass_clean.txt`, `Metadata.csv` |
| 3 | `5.microbial_taxa/data/`: `Genus.tsv`, `Metadata_tidy.csv`, and a supported phylum source described below |
| 4 | `6.microbial_function/data/`: `Metadata_tidy.csv`, `humann4_pathabundance_cpm.cleanIDs.tsv`, `pathabundance_merged.unstratified.clean.tsv` |
| 5 | The plant workbooks from Figure 1; functional metadata and CPM pathway abundances from Figure 4; `8.trait_function_linkage/data/`: `icon_seed.png`, `icon_root.png`, `icon_microbes.png`, `icon_membrane.png`, `icon_recovery.png` |

Figure 3 uses a complete `Genus_Phylum_mapping.csv` when supplied. Otherwise it uses `FIGURE3_PHYLUM_FILE`, then `Bracken_Phylum_long_with_metadata.csv`, then `Kraken2_Phylum_long_with_metadata.csv`. The mapping/table schema and QC are checked in the taxonomy script; the actual source and denominator are recorded with the outputs.

Required R packages are listed in `Run_All_Figures.R`. It reports missing packages before launching any selected figure. Bioconductor packages include DESeq2, ALDEx2, Maaslin2, AnnotationHub, AnnotationDbi, clusterProfiler and enrichplot; install these with `BiocManager::install()`.

## Larger labels and figure exports

The five main figures default to **8.0 inches wide**. PDF, PNG and TIFF use the same layout and white background; filenames end in `_LargeLabels`.

- **PDF:** vector output; Cairo embeds fonts when available, with the standard PDF device used as a fallback.
- **PNG:** 600 dpi.
- **TIFF:** 600 dpi with lossless LZW compression.

DPI applies to raster formats. An 8-inch raster export is 4800 pixels wide. Labels also scale when a figure is resized in Word; use `FONT_SCALE` to adjust the native text size independently of the export width.

| Figure | Default main export, inches | Current presentation |
| --- | --- | --- |
| 1 | 8.0 × 5.0 | Wider canvas for biomass, physiology and the experiment illustration |
| 2 | 8.0 × 8.05 | Shared PC1 scale in A/B with explicit −50, 0 and 50 ticks and padded limits; gene names removed from C/D |
| 3 | 8.0 × calculated height | A/B/D retain genus-based analyses; main C presents phylum composition; height accounts for the measured legends |
| 4 | 8.0 × 7.4 | Taller A/B and extra room for C's pathway labels |
| 5 | 8.0 × 6.0 | Wider network/data areas and smaller conceptual-model text with larger icons |

The Figure 3 PERMANOVA axis has room beyond 0.6 and the largest bar-end annotation. The data, model terms and significance thresholds remain defined by the analysis scripts.

### Standalone genus composition: Figure S4

Figure S4 uses **the top 50 genera plus Other**, selected from the QC-retained genus profiles. Other is the summed abundance of genera outside the selected set; the top genera are not renormalised to 100% after dropping Other.

The current layout has:

- Panels **A: Control, B: Polymer, C: 6-PP**, with equal plotting widths and one shared y-axis title.
- **Irrigated/Drought** columns and **During/Post Drought** rows.
- **Two columns of italic genus names on the right**, with Other in plain text.
- Smaller **Nil/AR37** labels: 6.0 points in the manuscript display calibration, reduced from 7.5.
- A separate **0.18-inch white margin above A–C**.
- A default width of **10.5 inches** and a minimum content height of **4.5 inches**; the top margin is added to the final height. The legend can increase the content height when necessary.

Set optional controls before sourcing the runner or taxonomy script:

```r
Sys.setenv(
  FIGURE3_TOP_N_GENERA = "50",
  FIGURE3_GENUS_LABEL_PT = "5.75",
  FIGURE3_GENUS_X_LABEL_PT = "6.0",
  FIGURE3_GENUS_WIDTH_IN = "10.50",
  FIGURE3_GENUS_BODY_HEIGHT_IN = "4.50",
  FIGURE3_GENUS_TOP_MARGIN_IN = "0.18"
)
```

The main phylum panel defaults to the top 10 named phyla; its controls are `FIGURE3_TOP_N_PHYLA` and `FIGURE3_PHYLUM_LABEL_PT`. The genus supplementary canvas is controlled separately from the shared main-figure export width.

Use the updated supplementary files named **`Figure3_Genus_Composition_Supplement_LargeLabels`** in `7.manuscript_figure/4.figure3/`.

### Output locations

| Output | Folder | Main filename stem |
| --- | --- | --- |
| Table 1 | `7.manuscript_figure/1.table1/` | Produced by the germination/endophyte script |
| Figure 1 | `7.manuscript_figure/2.figure1/` | `Figure1_NatCom_LargeLabels` |
| Figure 2 | `7.manuscript_figure/3.figure2/` | `Figure2_NatCom_LargeLabels` |
| Figure 3 | `7.manuscript_figure/4.figure3/` | `Figure3_NatCom_taxonomy_LargeLabels` |
| Figure S4 | `7.manuscript_figure/4.figure3/` | `Figure3_Genus_Composition_Supplement_LargeLabels` |
| Figure 4 | `7.manuscript_figure/5.figure4/` | `Figure4_NatCom_function_LargeLabels` |
| Figure 5 | `7.manuscript_figure/6.figure5/` | `Figure5_NatCom_linking_LargeLabels` |

Statistical tables, QC files and intermediate plots continue to be written under the corresponding analysis `results/` folders.

## Update this GitHub repository from PowerShell

Keep `README.md`, `Run_All_Figures.R` and `Update_GitHub_Repo.ps1` together in the extracted update folder. The five figure scripts are read from their current locations under the Windows project root.

Install [Git for Windows](https://gitforwindows.org/) and use your existing GitHub authentication. If Git identity is not already configured, set it once with your own email:

```powershell
git config --global user.name "Sereyboth Soth"
git config --global user.email "YOUR_GITHUB_EMAIL_OR_NOREPLY_EMAIL"
```

From the extracted update folder, run:

```powershell
powershell.exe -NoProfile -ExecutionPolicy Bypass -File ".\Update_GitHub_Repo.ps1"
```

This command **commits and pushes the selected updated files** to the existing GitHub repository. It includes the merged README, the runner, this PowerShell script and the five figure scripts. Tracked copies of the obsolete larger-label README at its known locations are removed as part of consolidation.

The script uses the project root if it is already a Git checkout. Otherwise it clones the existing GitHub repository into a sibling folder named `6PP_Ryegrass_Metatranscriptomics`, preserving the repository history, and copies the selected files there. It also installs the merged README and runner at the local project root.

For a different existing checkout, provide its path:

```powershell
.\Update_GitHub_Repo.ps1 -RepositoryPath "C:/path/to/6PP_Ryegrass_Metatranscriptomics"
```

Optional modes:

```powershell
.\Update_GitHub_Repo.ps1 -PreviewOnly  # Show the selected files without changing them.
.\Update_GitHub_Repo.ps1 -NoPush       # Make the commit locally, then inspect it before pushing.
```

The script verifies the fetch/push remote and branch, stages only its explicit file list, and uses a path-scoped commit. Data, result folders and figure images are outside that update list. It uses a normal push and stops on unresolved Git history or push errors; it does not force-push or discard local changes. If it finds pre-existing unpushed commits, reconcile them first or use a fresh checkout with `-RepositoryPath`.

Git Credential Manager may open a browser for GitHub sign-in. See [GitHub's authentication guidance](https://docs.github.com/en/get-started/git-basics/caching-your-github-credentials-in-git).

## License

MIT License; see `LICENSE`.
