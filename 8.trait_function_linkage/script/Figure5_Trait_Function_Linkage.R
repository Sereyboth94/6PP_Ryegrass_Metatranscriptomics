###############################################################################
# Figure 5 — Nature Communications
# Linking microbial pathways with plant performance
#
# Panels:
#   A = correlation network between plant traits and microbial pathways
#   B = Spearman correlation heatmap
#   C = conceptual model linking 6-PP seed coating to drought recovery
#
# Outputs:
#   Figure5_NatCom_linking.pdf
#   Figure5_NatCom_linking.png
#   Figure5_trait_pathway_correlations.csv
#   Figure5_trait_table_used.csv
#   Figure5_pathway_matrix_used.csv
#
# Notes:
#   - This script is independent. It reads the same source files used for
#     Figures 1 and 4.
#   - Fresh and dry biomass are merged as a single z-scored Plant biomass trait.
#   - By default, correlations are calculated for Drought samples from both
#     During Drought and Post Drought harvests. Change analysis_water or
#     analysis_harvest below if needed.
###############################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(stringr)
  library(matrixStats)
  library(cowplot)
  library(patchwork)
  library(grid)
  library(ggimage)
  library(magick)
})

# ─────────────────────────────────────────────────────────────────────────────
# 1. Paths and analysis settings
# ─────────────────────────────────────────────────────────────────────────────

project_root <- "C:/Users/soths/OneDrive - Lincoln University/Writing/NatMicro/github"

# Input files used by the HUMAnN4 microbial-function pipeline
function_data_dir <- file.path(project_root, "6.microbial_function", "data")

# Physiology files used by Figure 1
plant_data_dir <- file.path(project_root, "3.plant_biomass", "data")

# Figure 5-specific outputs
results_dir <- file.path(project_root, "8.trait_function_linkage", "results")
fig5_dir <- file.path(project_root, "7.manuscript_figure", "6.figure5")

# Optional icon files for conceptual model panel C
icon_dir <- file.path(project_root, "8.trait_function_linkage", "data")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig5_dir, recursive = TRUE, showWarnings = FALSE)

meta_fp <- file.path(function_data_dir, "Metadata_tidy.csv")
pathway_cpm_fp <- file.path(function_data_dir, "humann4_pathabundance_cpm.cleanIDs.tsv")

spad_dir <- plant_data_dir
biomass_dir <- plant_data_dir

# Correlation subset. For recovery-focused interpretation, Drought samples are used.
analysis_water <- "Drought"
analysis_harvest <- c("During Drought", "Post Drought")

# Number of microbial pathways shown in panels A and B.
top_n_pathways <- 16

# Network edge filter for panel A.
network_fdr_cutoff <- 0.05
network_abs_rho_cutoff <- 0.35

# ─────────────────────────────────────────────────────────────────────────────
# 2. Theme and colours
# ─────────────────────────────────────────────────────────────────────────────

theme_natcom <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 3),
      plot.subtitle = element_text(face = "bold", hjust = 0, size = base_size - 1),
      plot.margin = margin(8, 8, 8, 8)
    )
}

treat_cols <- c(
  "Control" = "#F8766D",
  "Polymer" = "#00BA38",
  "6-PP" = "#619CFF"
)

# ─────────────────────────────────────────────────────────────────────────────
# 3. Helper functions
# ─────────────────────────────────────────────────────────────────────────────

standardise_names <- function(df) {
  names(df) <- names(df) %>%
    tolower() %>%
    gsub("[^a-z0-9]+", "_", .) %>%
    gsub("^_|_$", "", .)
  df
}

clean_treatment <- function(x) {
  x <- stringr::str_squish(as.character(x))
  dplyr::recode(
    x,
    "VOC" = "6-PP",
    "6PP" = "6-PP",
    "6PP plus polymer" = "6-PP",
    "6-PP plus polymer" = "6-PP",
    "6PP+polymer" = "6-PP",
    "6-PP+polymer" = "6-PP",
    .default = x
  )
}

clean_harvest <- function(x) {
  x <- stringr::str_squish(as.character(x))
  dplyr::recode(
    x,
    "During drought" = "During Drought",
    "during drought" = "During Drought",
    "During_Drought" = "During Drought",
    "Post drought" = "Post Drought",
    "post drought" = "Post Drought",
    "Post_Drought" = "Post Drought",
    .default = x
  )
}

pathway_label <- function(x) {
  x %>%
    as.character() %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_replace_all("\\.", " ") %>%
    stringr::str_replace_all("−|–", "-") %>%
    stringr::str_squish() %>%
    
    # specific messy pathway fixes
    stringr::str_replace_all(".*LIPASYN[- ]?PWY.*", "Fatty acid biosyn. II") %>%
    stringr::str_replace_all(".*dodecenoate.*", "Dodecenoate biosyn. II") %>%
    stringr::str_replace_all(".*430.*cis.*vaccenate.*", "Cis-vaccenate biosyn.") %>%
    stringr::str_replace_all(".*cytosolic glycolysis.*pyruvate dehydrogenase.*TCA cycle.*",
                             "Glycolysis-PDH-TCA cycle") %>%
    stringr::str_replace_all(".*TCA[: ]+tca cycle i.*", "TCA: TCA cycle I") %>%
    stringr::str_replace_all(".*TCA cycle V.*|.*TCA cycle v.*", "TCA cycle V") %>%
    
    # known biological labels
    stringr::str_replace_all(".*fatty acid beta oxidation V.*", "FA beta-oxidation V") %>%
    stringr::str_replace_all(".*fatty acid beta oxidation IV.*", "FA beta-oxidation IV") %>%
    stringr::str_replace_all(".*fatty acid beta oxidation I.*", "FA beta-oxidation I") %>%
    stringr::str_replace_all(".*409.*purine.*salvage.*", "Purine nt. salvage") %>%
    stringr::str_replace_all(".*REDCITCYC.*TCA cycle VI.*", "TCA cycle VI") %>%
    stringr::str_replace_all(".*PWY[- ]?4661.*inositol.*", "Inositol phosphate biosyn.") %>%
    stringr::str_replace_all(".*PWY[- ]?6435.*", "4-hydroxybenzoate biosyn.") %>%
    stringr::str_replace_all(".*PWY[- ]?735.*", "L-tyrosine biosyn.") %>%
    stringr::str_replace_all(".*PWY[- ]?6628.*", "L-phenylalanine biosyn.") %>%
    stringr::str_replace_all(".*PWY[- ]?1422.*", "Vitamin E biosyn. (tocopherols)") %>%
    stringr::str_replace_all(".*PWY[- ]?7436.*", "Vitamin E biosyn. (tocotrienols)") %>%
    stringr::str_replace_all(".*PWY[- ]?8143.*", "Ascorbate biosyn.") %>%
    stringr::str_replace_all(".*PWY[- ]?7392.*", "Taxadiene biosyn.") %>%
    stringr::str_replace_all(".*PWY[- ]?2201.*", "Folate transformations") %>%
    stringr::str_replace_all(".*PWY[- ]?821.*", "Sulfur AA biosyn.") %>%
    stringr::str_replace_all(".*PWY[- ]?6609.*", "Adenosine salvage") %>%
    stringr::str_replace_all(".*PWY[- ]?6353.*", "Purine nt. degrad.") %>%
    stringr::str_replace_all(".*PWY[- ]?7254.*", "TCA cycle VII acetate prod.") %>%
    stringr::str_replace_all(".*PWY[- ]?241.*", "C4 photosyn. assim.") %>%
    stringr::str_replace_all(
      ".*palmitoleate.*dodec.*enoate.*",
      "Palmitoleate biosyn. I"
    ) %>%
    
    # remove remaining IDs/prefixes
    stringr::str_replace_all("^PWY[- ]?[0-9]+[: ]*", "") %>%
    stringr::str_replace_all("^P[0-9]+[- ]?PWY[: ]*", "") %>%
    stringr::str_replace_all("^[- ]?[0-9]+[: ]*[- ]*", "") %>%
    stringr::str_replace_all("^[- ]+", "") %>%
    
    # general shortening
    stringr::str_replace_all("superpathway of ", "") %>%
    stringr::str_replace_all("biosynthesis", "biosyn.") %>%
    stringr::str_replace_all("degradation", "degrad.") %>%
    stringr::str_replace_all("fermentation", "ferment.") %>%
    stringr::str_replace_all("nucleotide", "nt.") %>%
    stringr::str_replace_all("amino acid", "AA") %>%
    stringr::str_replace_all("\\(.*?\\)", "") %>%
    stringr::str_squish() %>%
    
    # sentence case, then restore acronyms and Roman numerals
    stringr::str_to_sentence() %>%
    stringr::str_replace_all("^[-–−]\\s*", "") %>%
    stringr::str_replace_all("\\bTca\\b", "TCA") %>%
    stringr::str_replace_all("\\btca\\b", "TCA") %>%
    stringr::str_replace_all("\\bPdh\\b", "PDH") %>%
    stringr::str_replace_all("\\bPwy\\b", "PWY") %>%
    stringr::str_replace_all("\\bFa\\b", "FA") %>%
    stringr::str_replace_all("\\bAa\\b", "AA") %>%
    stringr::str_replace_all("\\bNt\\b", "nt.") %>%
    stringr::str_replace_all("\\bI\\b$", "I") %>%
    stringr::str_replace_all("\\bIi\\b", "II") %>%
    stringr::str_replace_all("\\bIii\\b", "III") %>%
    stringr::str_replace_all("\\bIv\\b", "IV") %>%
    stringr::str_replace_all("\\bV\\b$", "V") %>%
    stringr::str_replace_all("\\bVi\\b", "VI") %>%
    stringr::str_replace_all("\\bVii\\b", "VII") %>%
    stringr::str_replace_all("^TCA: TCA cycle i$", "TCA: TCA cycle I") %>%
    stringr::str_replace_all("^TCA cycle v$", "TCA cycle V") %>%
    stringr::str_squish()
}

sig_stars <- function(q) {
  dplyr::case_when(
    is.na(q) ~ "",
    q < 0.001 ~ "***",
    q < 0.01  ~ "**",
    q < 0.05  ~ "*",
    TRUE ~ ""
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 4. Read HUMAnN4 pathway matrix and metadata
# ─────────────────────────────────────────────────────────────────────────────

stopifnot(file.exists(meta_fp), file.exists(pathway_cpm_fp))

meta2 <- read.csv(meta_fp, stringsAsFactors = FALSE) %>%
  dplyr::mutate(
    SampleID = as.character(SampleID),
    Treatment = clean_treatment(Treatment),
    Harvest = clean_harvest(Harvest),
    Endophyte = factor(Endophyte, levels = c("Nil", "AR37")),
    Treatment = factor(Treatment, levels = c("Control", "Polymer", "6-PP")),
    Water = factor(Water, levels = c("Irrigated", "Drought")),
    Harvest = factor(Harvest, levels = c("During Drought", "Post Drought"))
  )

path_cpm <- read.table(
  pathway_cpm_fp,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

drop_rows <- intersect(rownames(path_cpm), c("UNMAPPED", "UNINTEGRATED", "UNINTEGRATED|unclassified"))
if (length(drop_rows) > 0) path_cpm <- path_cpm[setdiff(rownames(path_cpm), drop_rows), , drop = FALSE]

meta2 <- meta2 %>% dplyr::filter(SampleID %in% colnames(path_cpm))
path_cpm <- path_cpm[, meta2$SampleID, drop = FALSE]
stopifnot(identical(colnames(path_cpm), meta2$SampleID))

X_cpm <- t(as.matrix(path_cpm))
X_cpm[is.na(X_cpm)] <- 0

# Add replicate if metadata lacks it. This matches the replicate structure used in
# the physiology files by ordering within each treatment/endophyte/water/harvest group.
if (!"Replicate" %in% names(meta2)) {
  meta2 <- meta2 %>%
    dplyr::group_by(Treatment, Endophyte, Water, Harvest) %>%
    dplyr::arrange(SampleID, .by_group = TRUE) %>%
    dplyr::mutate(Replicate = dplyr::row_number()) %>%
    dplyr::ungroup()
}

# ─────────────────────────────────────────────────────────────────────────────
# 5. Read plant physiology data and create Plant biomass trait
# ─────────────────────────────────────────────────────────────────────────────

get_timepoint <- function(path) {
  fn <- basename(path)
  dplyr::case_when(
    str_detect(fn, regex("before[_ ]?drought", ignore_case = TRUE)) ~ NA_character_,
    str_detect(fn, regex("during[_ ]?drought", ignore_case = TRUE)) ~ "During Drought",
    str_detect(fn, regex("post[_ ]?drought", ignore_case = TRUE)) ~ "Post Drought",
    TRUE ~ NA_character_
  )
}

read_spad_file <- function(path) {
  harvest_level <- get_timepoint(path)
  if (is.na(harvest_level)) return(NULL)
  df <- readxl::read_excel(path) %>% standardise_names()
  required <- c("treatment", "endophyte", "water", "replicate", "spad")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) stop("Missing columns in ", basename(path), ": ", paste(missing, collapse = ", "))
  df %>%
    transmute(
      Treatment = clean_treatment(treatment),
      Endophyte = stringr::str_squish(as.character(endophyte)),
      Water = stringr::str_squish(as.character(water)),
      Replicate = as.integer(replicate),
      Harvest = harvest_level,
      SPAD = as.numeric(spad)
    )
}

parse_biomass_file <- function(path) {
  fn <- basename(path)
  tibble(
    Metric = dplyr::case_when(
      str_detect(fn, regex("^fresh[_ ]?weight", ignore_case = TRUE)) ~ "Fresh_biomass",
      str_detect(fn, regex("^dry[_ ]?weight", ignore_case = TRUE)) ~ "Dry_biomass",
      TRUE ~ NA_character_
    ),
    Harvest = dplyr::case_when(
      str_detect(fn, regex("during[_ ]?drought", ignore_case = TRUE)) ~ "During Drought",
      str_detect(fn, regex("post[_ ]?drought", ignore_case = TRUE)) ~ "Post Drought",
      TRUE ~ NA_character_
    )
  )
}

read_biomass_file <- function(path) {
  meta_file <- parse_biomass_file(path)
  if (is.na(meta_file$Metric) || is.na(meta_file$Harvest)) return(NULL)
  df <- readxl::read_excel(path) %>% standardise_names()
  value_col <- dplyr::case_when(
    "freshweight" %in% names(df) ~ "freshweight",
    "fresh_weight" %in% names(df) ~ "fresh_weight",
    "dryweight" %in% names(df) ~ "dryweight",
    "dry_weight" %in% names(df) ~ "dry_weight",
    TRUE ~ NA_character_
  )
  if (is.na(value_col)) stop("Could not find biomass value column in: ", basename(path))
  df %>%
    transmute(
      Treatment = clean_treatment(treatment),
      Endophyte = stringr::str_squish(as.character(endophyte)),
      Water = stringr::str_squish(as.character(water)),
      Replicate = as.integer(replicate),
      Harvest = meta_file$Harvest,
      Metric = meta_file$Metric,
      Value = as.numeric(.data[[value_col]])
    )
}

spad_files <- list.files(spad_dir, pattern = "\\.xlsx$", full.names = TRUE)
spad_files <- spad_files[
  grepl("spad|before|during|post", basename(spad_files), ignore.case = TRUE) &
    !grepl("fresh|dry|weight|biomass", basename(spad_files), ignore.case = TRUE)
]

biomass_files <- list.files(biomass_dir, pattern = "\\.xlsx$", full.names = TRUE)
biomass_files <- biomass_files[
  grepl("fresh|dry|weight", basename(biomass_files), ignore.case = TRUE)
]

if (length(spad_files) == 0) stop("No SPAD .xlsx files found in: ", spad_dir)
if (length(biomass_files) == 0) stop("No biomass .xlsx files found in: ", biomass_dir)

spad_raw <- purrr::map_dfr(spad_files, read_spad_file) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Control", "Polymer", "6-PP")),
    Endophyte = factor(Endophyte, levels = c("Nil", "AR37")),
    Water = factor(Water, levels = c("Irrigated", "Drought")),
    Harvest = factor(Harvest, levels = c("During Drought", "Post Drought"))
  )

biomass_raw <- purrr::map_dfr(biomass_files, read_biomass_file) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("Control", "Polymer", "6-PP")),
    Endophyte = factor(Endophyte, levels = c("Nil", "AR37")),
    Water = factor(Water, levels = c("Irrigated", "Drought")),
    Harvest = factor(Harvest, levels = c("During Drought", "Post Drought"))
  )

biomass_wide <- biomass_raw %>%
  group_by(Treatment, Endophyte, Water, Harvest, Replicate, Metric) %>%
  summarise(Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Metric, values_from = Value)

plant_traits_design <- spad_raw %>%
  group_by(Treatment, Endophyte, Water, Harvest, Replicate) %>%
  summarise(SPAD = mean(SPAD, na.rm = TRUE), .groups = "drop") %>%
  left_join(
    biomass_wide,
    by = c("Treatment", "Endophyte", "Water", "Harvest", "Replicate")
  ) %>%
  mutate(
    Plant_biomass = as.numeric(scale(Fresh_biomass)) + as.numeric(scale(Dry_biomass))
  )

plant_traits <- meta2 %>%
  dplyr::mutate(
    Treatment = factor(Treatment, levels = c("Control", "Polymer", "6-PP")),
    Endophyte = factor(Endophyte, levels = c("Nil", "AR37")),
    Water = factor(Water, levels = c("Irrigated", "Drought")),
    Harvest = factor(Harvest, levels = c("During Drought", "Post Drought"))
  ) %>%
  dplyr::left_join(
    plant_traits_design %>%
      dplyr::group_by(Treatment, Endophyte, Water, Harvest) %>%
      dplyr::summarise(
        SPAD = mean(SPAD, na.rm = TRUE),
        Plant_biomass = mean(Plant_biomass, na.rm = TRUE),
        .groups = "drop"
      ),
    by = c("Treatment", "Endophyte", "Water", "Harvest")
  )

if (any(is.na(plant_traits$SPAD)) || any(is.na(plant_traits$Plant_biomass))) {
  warning("Some samples did not receive SPAD or Plant_biomass values. Check replicate matching.")
}

write_csv(plant_traits, file.path(results_dir, "Figure5_trait_table_used.csv"))

# ─────────────────────────────────────────────────────────────────────────────
# 6. Build correlation table
# ─────────────────────────────────────────────────────────────────────────────

analysis_samples <- plant_traits %>%
  filter(Water %in% analysis_water, Harvest %in% analysis_harvest) %>%
  filter(!is.na(SPAD), !is.na(Plant_biomass))

if (nrow(analysis_samples) < 8) stop("Too few matched samples for correlation analysis.")

X_analysis <- X_cpm[analysis_samples$SampleID, , drop = FALSE]
# Remove uninformative pathways
keep_pathways <- matrixStats::colVars(X_analysis) > 0
X_analysis <- X_analysis[, keep_pathways, drop = FALSE]

# Z-score pathways across samples for correlation and display.
X_z <- scale(X_analysis)
X_z[is.na(X_z)] <- 0

plant_mat <- analysis_samples %>%
  select(SampleID, SPAD, Plant_biomass) %>%
  column_to_rownames("SampleID")

cor_one <- function(trait_name, pathway_name) {
  x <- plant_mat[[trait_name]]
  y <- X_z[, pathway_name]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 6 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(tibble(Trait = trait_name, Pathway = pathway_name, rho = NA_real_, p = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
  tibble(Trait = trait_name, Pathway = pathway_name, rho = unname(ct$estimate), p = ct$p.value)
}

cor_tbl <- expand_grid(
  Trait = c("SPAD", "Plant_biomass"),
  Pathway = colnames(X_z)
) %>%
  pmap_dfr(~ cor_one(..1, ..2)) %>%
  mutate(
    q = p.adjust(p, method = "BH"),
    abs_rho = abs(rho),
    Pathway_label_raw = pathway_label(Pathway)
  )

# Keep one display label per pathway. If several features clean to the same name,
# preserve the strongest association to avoid duplicated heatmap labels.
pathway_rank <- cor_tbl %>%
  group_by(Pathway, Pathway_label_raw) %>%
  summarise(score = max(abs_rho, na.rm = TRUE), min_q = min(q, na.rm = TRUE), .groups = "drop") %>%
  arrange(min_q, desc(score))

selected_pathways <- pathway_rank %>%
  group_by(Pathway_label_raw) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  arrange(min_q, desc(score)) %>%
  slice_head(n = top_n_pathways) %>%
  pull(Pathway)

cor_plot <- cor_tbl %>%
  filter(Pathway %in% selected_pathways) %>%
  mutate(
    Trait_label = recode(Trait, "SPAD" = "SPAD", "Plant_biomass" = "Plant biomass"),
    Pathway_label = pathway_label(Pathway),
    stars = sig_stars(q)
  )

# Unique pathway label order based on strongest association.
path_order <- cor_plot %>%
  group_by(Pathway_label) %>%
  summarise(score = max(abs_rho, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(score)) %>%
  pull(Pathway_label)

cor_plot <- cor_plot %>%
  mutate(
    Trait_label = factor(Trait_label, levels = c("SPAD", "Plant biomass")),
    Pathway_label = factor(Pathway_label, levels = rev(path_order))
  )

write_csv(cor_tbl, file.path(results_dir, "Figure5_trait_pathway_correlations.csv"))
write_csv(as.data.frame(X_z) %>% rownames_to_column("SampleID"), file.path(results_dir, "Figure5_pathway_matrix_used.csv"))

# ─────────────────────────────────────────────────────────────────────────────
# 7. Panel A — correlation network
# ─────────────────────────────────────────────────────────────────────────────

nodes_trait <- tibble(
  name = c("SPAD", "Plant biomass"),
  label = c("SPAD", "Plant\nbiomass"),
  type = "Plant trait",
  x = 0,
  y = c(0.70, 0.30)
)

nodes_path <- tibble(
  name = rev(path_order),
  label = name,
  type = "Microbial pathway",
  x = 1.40,
  y = seq(0.95, 0.05, length.out = length(path_order))
)

edges <- cor_plot %>%
  filter(!is.na(q), q < network_fdr_cutoff, abs_rho >= network_abs_rho_cutoff) %>%
  mutate(
    trait_name = as.character(Trait_label),
    pathway_name = as.character(Pathway_label)
  ) %>%
  left_join(nodes_trait %>% select(trait_name = name, x1 = x, y1 = y), by = "trait_name") %>%
  left_join(nodes_path %>% select(pathway_name = name, x2 = x, y2 = y), by = "pathway_name")

if (nrow(edges) == 0) {
  edges <- cor_plot %>%
    group_by(Trait_label) %>%
    slice_max(abs_rho, n = 6, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      trait_name = as.character(Trait_label),
      pathway_name = as.character(Pathway_label)
    ) %>%
    left_join(nodes_trait %>% select(trait_name = name, x1 = x, y1 = y), by = "trait_name") %>%
    left_join(nodes_path %>% select(pathway_name = name, x2 = x, y2 = y), by = "pathway_name")
}

panel_a <- ggplot() +
  geom_segment(
    data = edges,
    aes(x = x1 + 0.05, y = y1, xend = x2 - 0.03, yend = y2,
        colour = rho, linewidth = abs_rho),
    alpha = 0.78
  ) +
  geom_label(
    data = nodes_trait,
    aes(x = x, y = y, label = label),
    fill = "#E5F5E0",
    colour = "#00441B",
    linewidth = 0.4,
    fontface = "bold",
    size = 5.2,
    label.padding = unit(0.35, "lines")
  ) +
  geom_point(
    data = nodes_path,
    aes(x = x, y = y),
    shape = 21,
    fill = "#EFF3FF",
    colour = "#08519C",
    size = 4.5,
    stroke = 0.7
  ) +
  geom_text(
    data = nodes_path,
    aes(x = x + 0.1, y = y, label = label),
    hjust = 0,
    size = 4.1,
    colour = "black"
  ) +
  annotate("text", x = 0, y = 1.04, label = "Plant traits", colour = "#00441B", fontface = "bold", size = 5.2) +
  annotate("text", x = 1.68, y = 1.04, label = "Microbial pathways", colour = "#08519C", fontface = "bold", size = 5.2) +
  scale_colour_gradient2(
    low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
    limits = c(-0.8, 0.8), name = "Spearman ρ"
  ) +
  scale_linewidth_continuous(name = "|Spearman ρ|", range = c(0.25, 1.8)) +
  labs(
    title = "A",
    subtitle = "Correlation network linking microbial pathways with plant performance"
  ) +
  coord_cartesian(xlim = c(-0.10, 2.10), ylim = c(-0.03, 1.08), clip = "off") +
  theme_void(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0),
    plot.subtitle = element_text(face = "bold", size = 14, hjust = 0),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11),
    plot.margin = margin(8, 8, 8, 8)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 8. Panel B — Spearman heatmap
# ─────────────────────────────────────────────────────────────────────────────
panel_b <- ggplot(
  cor_plot,
  aes(x = rho, y = Pathway_label, colour = Trait_label)
) +
  geom_vline(xintercept = 0, linewidth = 0.4, linetype = "dashed", colour = "grey50") +
  geom_point(
    position = position_dodge(width = 0.55),
    size = 3.8,
    alpha = 0.95
  ) +
  scale_colour_manual(
    values = c(
      "SPAD" = "#1F78B4",
      "Plant biomass" = "#33A02C"
    ),
    name = "Plant trait"
  ) +
  scale_x_continuous(
    limits = c(-0.9, 0.9),
    breaks = seq(-0.8, 0.8, 0.4)
  ) +
  labs(
    title = "B",
    subtitle = "Treatment-level associations between plant traits and microbial pathways",
    x = "Spearman rho",
    y = "Microbial pathways"
  ) +
  theme_natcom(15) +
  theme(
    panel.background = element_blank(),
    plot.background  = element_blank(),
    panel.grid       = element_blank(),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12),
    plot.subtitle = element_text(size = 13, face = "bold"),
    plot.margin = margin(8, 8, 8, 8)
  )


# ─────────────────────────────────────────────────────────────────────────────
# 9. Panel C — conceptual model
# ─────────────────────────────────────────────────────────────────────────────

concept_nodes <- tibble(
  step = 1:5,
  x = seq(0.08, 0.92, length.out = 5),
  y_icon = 0.55,
  y_box = 0.17,
  title = c(
    "6-PP seed coating",
    "Host signalling",
    "Microbial functional\nactivation",
    "Osmoprotection and\ncellular protection",
    "Improved drought\nrecovery"
  ),
  icon_path = c(
    file.path(icon_dir, "icon_seed.png"),
    file.path(icon_dir, "icon_root.png"),
    file.path(icon_dir, "icon_microbes.png"),
    file.path(icon_dir, "icon_membrane.png"),
    file.path(icon_dir, "icon_recovery.png")
  ),
  detail = c(
    "6-PP-treated seed establishes\nbeneficial rhizosphere signals",
    
    "Enhanced root exudation and\ncommunication",
    
    "Activation of key pathways:\nFA beta-oxidation, nucleotide salvage\nand antioxidant biosynthesis",
    
    "Osmolyte accumulation,\nmembrane stability and redox balance",
    
    "Higher SPAD, plant biomass and\npost-drought recovery"
  ),
  col = c("#6A3D9A", "#006D2C", "#08519C", "#D95F02", "#006D2C"),
  fill = c("#F3E5F5", "#E5F5E0", "#EFF3FF", "#FFF5EB", "#E5F5E0")
)

if (!all(file.exists(concept_nodes$icon_path))) {
  warning(
    "Some Panel C icon files are missing. Expected files:\n",
    paste(concept_nodes$icon_path, collapse = "\n")
  )
}

arrow_df <- tibble(
  x = concept_nodes$x[-nrow(concept_nodes)] + 0.085,
  xend = concept_nodes$x[-1] - 0.085,
  y = 0.55,
  yend = 0.55
)

panel_c <- ggplot() +
  geom_segment(
    data = arrow_df,
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.22, "inches")),
    linewidth = 0.9,
    colour = "black"
  ) +
  geom_text(
    data = concept_nodes,
    aes(x = x, y = 0.87, label = title, colour = col),
    fontface = "bold",
    size = 5.0,
    lineheight = 0.95,
    show.legend = FALSE
  ) +
  ggimage::geom_image(
    data = concept_nodes,
    aes(x = x, y = y_icon, image = icon_path),
    size = 0.40,
    asp = 1
  ) +
  geom_label(
    data = concept_nodes,
    aes(x = x, y = y_box, label = detail, colour = col, fill = fill),
    fontface = "bold",
    size = 4.2,
    lineheight = 0.90,
    linewidth = 0.35,
    label.padding = unit(0.55, "lines"),
    show.legend = FALSE
  ) +
  scale_colour_identity() +
  scale_fill_identity() +
  labs(
    title = "C",
    subtitle = "Conceptual model: 6-PP seed coating enhances drought recovery via microbiome functions"
  ) +
  coord_cartesian(
    xlim = c(0.01, 0.99),
    ylim = c(0.02, 0.98),
    clip = "off"
  ) +
  theme_void(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 24, hjust = 0),
    plot.subtitle = element_text(face = "bold", size = 15, hjust = 0),
    plot.margin = margin(8, 8, 8, 8)
  )

# ─────────────────────────────────────────────────────────────────────────────
# 10. Assemble Figure 5
# ─────────────────────────────────────────────────────────────────────────────

top <- cowplot::plot_grid(panel_a, panel_b, nrow = 1, rel_widths = c(1.20, 1.10))

fig5 <- cowplot::plot_grid(
  top,
  panel_c,
  ncol = 1,
  rel_heights = c(1.15, 0.75)
)

# Save

ggsave(
  file.path(fig5_dir, "Figure5_NatCom_linking.pdf"),
  fig5,
  width = 22,
  height = 14,
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  file.path(fig5_dir, "Figure5_NatCom_linking.png"),
  fig5,
  width = 22,
  height = 14,
  dpi = 600,
  bg = "white", 
  limitsize = FALSE
)

message("Figure 5 saved to: ", fig5_dir)
message("Correlation subset: Water = ", paste(analysis_water, collapse = ", "),
        "; Harvest = ", paste(analysis_harvest, collapse = ", "))
