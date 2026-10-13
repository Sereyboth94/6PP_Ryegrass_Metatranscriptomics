# Run the five manuscript figure scripts from their analysis folders.
# Runner revision: 2026-10-13 (project-relative paths and Figure S4 output checks).
# Place this file at the project root; MANUSCRIPT_PROJECT_ROOT can override that root.

find_script_directory <- function() {
  frames <- rev(sys.frames())
  for (frame in frames) {
    if (!is.null(frame$ofile)) {
      return(dirname(normalizePath(frame$ofile, mustWork = TRUE)))
    }
  }
  args <- commandArgs(trailingOnly = FALSE)
  arg <- args[grepl("^--file=", args)]
  if (length(arg)) {
    return(dirname(normalizePath(sub("^--file=", "", arg[1]), mustWork = TRUE)))
  }
  normalizePath(getwd(), mustWork = TRUE)
}

RUNNER_DIR <- find_script_directory()
PROJECT_ROOT <- Sys.getenv("MANUSCRIPT_PROJECT_ROOT", RUNNER_DIR)
if (!nzchar(PROJECT_ROOT)) PROJECT_ROOT <- RUNNER_DIR

FONT_SCALE <- 1.0       # Use 1.10 for a further 10% increase in manuscript labels.
SYMBOL_SCALE <- 0.80
DISPLAY_WIDTH_IN <- 6.27
EXPORT_WIDTH_IN <- 8.0
RUN_FIGURES <- 1:5
RUNNER_REVISION <- "2026-10-13 project-folder runner update"

if (!dir.exists(PROJECT_ROOT)) {
  stop("Project folder does not exist. Set MANUSCRIPT_PROJECT_ROOT: ", PROJECT_ROOT)
}
PROJECT_ROOT <- normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE)
if (!is.numeric(FONT_SCALE) || length(FONT_SCALE) != 1L ||
    !is.finite(FONT_SCALE) || FONT_SCALE <= 0 ||
    !is.numeric(DISPLAY_WIDTH_IN) || length(DISPLAY_WIDTH_IN) != 1L ||
    !is.finite(DISPLAY_WIDTH_IN) || DISPLAY_WIDTH_IN <= 0 ||
    !is.numeric(SYMBOL_SCALE) || length(SYMBOL_SCALE) != 1L ||
    !is.finite(SYMBOL_SCALE) || SYMBOL_SCALE <= 0 ||
    !is.numeric(EXPORT_WIDTH_IN) || length(EXPORT_WIDTH_IN) != 1L ||
    !is.finite(EXPORT_WIDTH_IN) || EXPORT_WIDTH_IN <= 0) {
  stop("FONT_SCALE, DISPLAY_WIDTH_IN, EXPORT_WIDTH_IN and SYMBOL_SCALE must be positive numbers.")
}
if (!is.numeric(RUN_FIGURES) || length(RUN_FIGURES) == 0L || anyNA(RUN_FIGURES) ||
    any(!RUN_FIGURES %in% 1:5)) {
  stop("RUN_FIGURES must contain one or more figure numbers from 1 to 5.")
}
RUN_FIGURES <- unique(as.integer(RUN_FIGURES))

script_relative_paths <- c(
  "3.plant_biomass/script/Fig.1_Generated_Script.R",
  "4.plant_transcriptome/script/Figure2_Plant_Transcriptome.R",
  "5.microbial_taxa/script/Figure3_Taxonomy.R",
  "6.microbial_function/scripts/Figure4_Function.R",
  "8.trait_function_linkage/script/Figure5_Trait_Function_Linkage.R"
)
script_paths <- file.path(PROJECT_ROOT, script_relative_paths)
script_names <- basename(script_relative_paths)
message("Manuscript runner revision: ", RUNNER_REVISION)
message("Project root: ", PROJECT_ROOT)
missing_scripts <- script_paths[RUN_FIGURES][!file.exists(script_paths[RUN_FIGURES])]
if (length(missing_scripts)) {
  stop("Required figure scripts are missing from the analysis folders:\n",
       paste(missing_scripts, collapse = "\n"))
}

# Parse every selected script before starting any analysis.
# Individual figures can have different layout revision dates.
for (i in RUN_FIGURES) {
  revision_lines <- grep(
    "^FIGURE_REVISION", readLines(script_paths[i], warn = FALSE, encoding = "UTF-8"),
    value = TRUE
  )
  if (length(revision_lines)) message("Figure ", i, ": ", revision_lines[1L])
  tryCatch(parse(file = script_paths[i], encoding = "UTF-8"), error = function(e) {
    stop("Syntax check failed for ", script_names[i], ": ", conditionMessage(e))
  })
}

fixed_inputs <- list(
  file.path("3.plant_biomass", "data", "Biotron_Experiment.png"),
  file.path("4.plant_transcriptome", "data", c("featureCounts_ryegrass_clean.txt", "Metadata.csv")),
  file.path("5.microbial_taxa", "data", c("Genus.tsv", "Metadata_tidy.csv")),
  file.path("6.microbial_function", "data", c(
    "Metadata_tidy.csv", "humann4_pathabundance_cpm.cleanIDs.tsv",
    "pathabundance_merged.unstratified.clean.tsv"
  )),
  c(
    file.path("6.microbial_function", "data", c(
      "Metadata_tidy.csv", "humann4_pathabundance_cpm.cleanIDs.tsv"
    )),
    file.path("8.trait_function_linkage", "data", c(
      "icon_seed.png", "icon_root.png", "icon_microbes.png", "icon_membrane.png", "icon_recovery.png"
    ))
  )
)
input_paths <- unique(file.path(PROJECT_ROOT, unlist(fixed_inputs[RUN_FIGURES])))
missing_inputs <- input_paths[!file.exists(input_paths)]
if (any(RUN_FIGURES %in% c(1, 5))) {
  plant_dir <- file.path(PROJECT_ROOT, "3.plant_biomass", "data")
  xlsx <- list.files(plant_dir, pattern = "\\.xlsx$", ignore.case = TRUE)
  spad <- xlsx[
    grepl("spad|before|during|post", xlsx, ignore.case = TRUE) &
    !grepl("fresh|dry|weight|biomass", xlsx, ignore.case = TRUE)
  ]
  biomass <- xlsx[grepl("fresh|dry|weight", xlsx, ignore.case = TRUE)]
  if (!length(spad)) missing_inputs <- c(missing_inputs, paste0(plant_dir, "/[SPAD .xlsx files]"))
  if (!length(biomass)) missing_inputs <- c(missing_inputs, paste0(plant_dir, "/[fresh/dry-weight .xlsx files]"))
}
if (3L %in% RUN_FIGURES) {
  taxa_dir <- file.path(PROJECT_ROOT, "5.microbial_taxa", "data")
  phylum_mapping <- file.path(taxa_dir, "Genus_Phylum_mapping.csv")
  phylum_override <- Sys.getenv("FIGURE3_PHYLUM_FILE", "")
  if (!file.exists(phylum_mapping)) {
    phylum_candidates <- if (nzchar(phylum_override)) {
      phylum_override
    } else {
      file.path(taxa_dir, c("Bracken_Phylum_long_with_metadata.csv",
                           "Kraken2_Phylum_long_with_metadata.csv"))
    }
    if (!any(file.exists(phylum_candidates))) {
      missing_inputs <- c(missing_inputs, paste0(
        taxa_dir, "/[Genus_Phylum_mapping.csv or a supported phylum table; see README]"
      ))
    }
  }
}

if (length(missing_inputs)) {
  stop("Required input data are missing:\n", paste(missing_inputs, collapse = "\n"))
}

required_packages <- list(
  c("readxl", "tidyverse", "stringr", "emmeans", "multcomp", "patchwork", "cowplot", "png", "magick"),
  c("DESeq2", "tidyverse", "pheatmap", "RColorBrewer", "matrixStats", "ggplot2",
    "ggrepel", "AnnotationHub", "AnnotationDbi", "tibble", "stringr", "here",
    "patchwork", "vegan", "BiocManager", "clusterProfiler", "cowplot", "enrichplot"),
  c("dplyr", "tidyr", "readr", "stringr", "tibble", "forcats", "ggplot2", "ggrepel",
    "patchwork", "scales", "vegan", "ape", "ALDEx2", "pheatmap", "RColorBrewer", "ggsci", "here", "cowplot"),
  c("tidyverse", "vegan", "ape", "pheatmap", "matrixStats", "ggrepel", "here",
    "BiocManager", "Maaslin2", "cowplot", "forcats"),
  c("tidyverse", "readxl", "stringr", "matrixStats", "cowplot", "patchwork", "ggimage", "magick")
)
required <- unique(unlist(required_packages[RUN_FIGURES]))
missing_packages <- setdiff(required, rownames(installed.packages()))
if (length(missing_packages)) {
  stop("Install the missing R packages before rerunning:\n",
       paste(missing_packages, collapse = ", "),
       "\nDESeq2, ALDEx2, Maaslin2, AnnotationHub, AnnotationDbi, clusterProfiler and enrichplot use Bioconductor.")
}

Sys.setenv(
  MANUSCRIPT_PROJECT_ROOT = normalizePath(PROJECT_ROOT, winslash = "/", mustWork = TRUE),
  MANUSCRIPT_FONT_SCALE = as.character(FONT_SCALE),
  MANUSCRIPT_SYMBOL_SCALE = as.character(SYMBOL_SCALE),
  MANUSCRIPT_FIGURE_WIDTH_IN = as.character(DISPLAY_WIDTH_IN),
  MANUSCRIPT_EXPORT_WIDTH_IN = as.character(EXPORT_WIDTH_IN)
)

rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
if (!file.exists(rscript)) stop("Could not find the Rscript executable: ", rscript)

output_dirs <- file.path(PROJECT_ROOT, "7.manuscript_figure", c(
  "2.figure1", "3.figure2", "4.figure3", "5.figure4", "6.figure5"
))
output_stems <- paste0(c(
  "Figure1_NatCom", "Figure2_NatCom", "Figure3_NatCom_taxonomy",
  "Figure4_NatCom_function", "Figure5_NatCom_linking"
), "_LargeLabels")

failed <- integer()
for (i in RUN_FIGURES) {
  message("\nRunning Figure ", i, ": ", script_names[i])
  flush.console()
  status <- tryCatch(
    system2(rscript, args = shQuote(script_paths[i]), stdout = "", stderr = ""),
    error = function(e) { message(conditionMessage(e)); 1L }
  )
  stems <- output_stems[i]
  if (i == 3L) {
    stems <- c(stems, "Figure3_Genus_Composition_Supplement_LargeLabels")
  }
  expected <- file.path(output_dirs[i], unlist(lapply(stems, function(stem) {
    paste0(stem, c(".pdf", ".png", ".tiff"))
  })))
  if (status != 0L || !all(file.exists(expected))) {
    failed <- c(failed, i)
    message("Figure ", i, " failed. See its error above; continuing with the remaining figures.")
    missing_outputs <- expected[!file.exists(expected)]
    if (length(missing_outputs)) {
      message("Missing output files:\n", paste(missing_outputs, collapse = "\n"))
    }
  } else {
    message("Figure ", i, " complete: PDF, PNG and TIFF in ", output_dirs[i])
  }
}
if (length(failed)) stop("The rerun failed for Figure(s): ", paste(failed, collapse = ", "))
message("\nAll selected figures regenerated successfully.")
