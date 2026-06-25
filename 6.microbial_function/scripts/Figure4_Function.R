# ==========================================================
# HUMAnN4 Pathway Analysis (Publication-ready version)
# - Bray-Curtis + Aitchison beta diversity
# - PERMANOVA (main + 2way + pairwise)
# - Dispersion tests for all factors
# - PCoA + NMDS
# - Heatmap (top variable pathways)
# - Differential pathways with MaAsLin2
# - Bubble plots from MaAsLin2
# - Pathway class scoring for mechanistic interpretation
# ==========================================================

# =========================
# Packages
# =========================
pkgs <- c("tidyverse", "vegan", "ape", "pheatmap", "matrixStats", "ggrepel", "here")

to_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("Maaslin2", quietly = TRUE)) BiocManager::install("Maaslin2")

invisible(lapply(pkgs, library, character.only = TRUE))
library(Maaslin2)

# =========================
# Paths
# =========================
project_dir <- "C:/Users/soths/OneDrive - Lincoln University/Writing/NatMicro/github"

data_dir <- file.path(project_dir, "6.microbial_function", "data")
out_dir  <- file.path(project_dir, "6.microbial_function", "results")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

dist_dir    <- file.path(out_dir, "01_Distances")
perm_dir    <- file.path(out_dir, "02_PERMANOVA")
ord_dir     <- file.path(out_dir, "03_Ordination")
heat_dir    <- file.path(out_dir, "04_Heatmaps")
maaslin_dir <- file.path(out_dir, "05_MaAsLin2")
bubble_dir  <- file.path(out_dir, "06_BubblePlots")
focus_dir   <- file.path(out_dir, "07_FocusedPathways")
class_dir   <- file.path(out_dir, "08_PathwayClassScoring")
repr_dir    <- file.path(out_dir, "09_Reproducibility")

dir.create(dist_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(perm_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(ord_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(heat_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(maaslin_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(bubble_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(focus_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(class_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(repr_dir,    recursive = TRUE, showWarnings = FALSE)

meta_fp         <- file.path(data_dir, "Metadata_tidy.csv")
pathway_cpm_fp  <- file.path(data_dir, "humann4_pathabundance_cpm.cleanIDs.tsv")
pathway_raw_fp  <- file.path(data_dir, "pathabundance_merged.unstratified.clean.tsv")

stopifnot(
  file.exists(meta_fp),
  file.exists(pathway_cpm_fp),
  file.exists(pathway_raw_fp)
)

# =========================
# Helper functions
# =========================
sanitize_string <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9]+", "_") %>%
    stringr::str_replace_all("^_+|_+$", "")
}

save_gg <- function(p, filename_base, out_path = out_dir, w = 10, h = 6, dpi = 300) {
  ggsave(file.path(out_path, paste0(filename_base, ".pdf")), p, width = w, height = h, useDingbats = FALSE)
  ggsave(file.path(out_path, paste0(filename_base, ".png")), p, width = w, height = h, dpi = dpi)
}

# =========================
# Read data
# =========================
meta <- read.csv(meta_fp, stringsAsFactors = FALSE)
stopifnot("SampleID" %in% names(meta))

meta$SampleID <- as.character(meta$SampleID)
if (anyDuplicated(meta$SampleID)) {
  dup <- meta$SampleID[duplicated(meta$SampleID)]
  stop("Duplicate SampleID in metadata: ", paste(unique(dup), collapse = ", "))
}

path <- read.table(
  pathway_cpm_fp,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

path_raw <- read.table(
  pathway_raw_fp,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

message("✅ Data loaded successfully")
message("Samples: ", nrow(meta))
message("Pathways (CPM): ", nrow(path))
message("Pathways (RAW): ", nrow(path_raw))

# =========================
# Clean pathway table
# =========================
drop_rows <- intersect(
  rownames(path),
  c("UNMAPPED", "UNINTEGRATED", "UNINTEGRATED|unclassified")
)

if (length(drop_rows) > 0) {
  message("Dropping rows: ", paste(drop_rows, collapse = ", "))
  path <- path[setdiff(rownames(path), drop_rows), , drop = FALSE]
}

# Align samples
meta2 <- meta %>% dplyr::filter(SampleID %in% colnames(path))
path2 <- path[, meta2$SampleID, drop = FALSE]

stopifnot(identical(colnames(path2), meta2$SampleID))
message("N samples: ", ncol(path2))
message("N pathways before filtering: ", nrow(path2))

# Ensure factors
meta2 <- meta2 %>%
  dplyr::mutate(
    Endophyte = factor(Endophyte, levels = c("Nil", "AR37")),
    Treatment = factor(Treatment, levels = c("Control", "Polymer", "6-PP")),
    Water     = factor(Water, levels = c("Irrigated", "Drought")),
    Harvest   = factor(Harvest, levels = c("During Drought", "Post Drought"))
  )

# =========================
# Matrices
# =========================
# samples x pathways
X <- t(as.matrix(path2))

# remove all-zero pathways
keep_feat <- colSums(X, na.rm = TRUE) > 0
X <- X[, keep_feat, drop = FALSE]

message("N pathways after removing all-zero: ", ncol(X))

# ensure sample order
stopifnot(all(meta2$SampleID %in% rownames(X)))
X <- X[meta2$SampleID, , drop = FALSE]
stopifnot(identical(rownames(X), meta2$SampleID))

# CPM matrix for ecological analyses
X_cpm <- X
X_cpm[is.na(X_cpm)] <- 0

# Pseudocount + CLR for compositional analysis
pseudocount <- 1
X_clr <- log(X_cpm + pseudocount)
X_clr <- sweep(X_clr, 1, rowMeans(X_clr), FUN = "-")

# =========================
# Distances
# =========================
set.seed(1)
d_bray  <- vegdist(X_cpm, method = "bray")
d_aitch <- dist(X_clr, method = "euclidean")

saveRDS(d_bray,  file.path(dist_dir, "dist_bray.rds"))
saveRDS(d_aitch, file.path(dist_dir, "dist_aitchison.rds"))

# ==========================================================
# PERMANOVA
# ==========================================================
perm_formula_main <- as.formula("dist_obj ~ Endophyte + Treatment + Water + Harvest")
perm_formula_2way <- as.formula("dist_obj ~ (Endophyte + Treatment + Water + Harvest)^2")

run_permanova <- function(dist_obj, meta_df, label) {
  perm_main <- adonis2(perm_formula_main, data = meta_df, permutations = 999, by = "margin")
  perm_2way <- adonis2(perm_formula_2way, data = meta_df, permutations = 999, by = "margin")
  
  write.table(
    as.data.frame(perm_main),
    file = file.path(perm_dir, paste0("PERMANOVA_main_effects_", label, ".tsv")),
    sep = "\t", quote = FALSE, col.names = NA
  )
  
  write.table(
    as.data.frame(perm_2way),
    file = file.path(perm_dir, paste0("PERMANOVA_2way_interactions_", label, ".tsv")),
    sep = "\t", quote = FALSE, col.names = NA
  )
  
  list(main = perm_main, twoway = perm_2way)
}

perm_bray  <- run_permanova(d_bray,  meta2, "Bray")
perm_aitch <- run_permanova(d_aitch, meta2, "Aitchison")

# =========================
# Pairwise PERMANOVA helper
# =========================
pairwise_adonis2 <- function(dist_obj, meta_df, factor_name, permutations = 999) {
  grp <- droplevels(factor(meta_df[[factor_name]]))
  levs <- levels(grp)
  if (length(levs) < 2) return(tibble())
  
  combs <- combn(levs, 2, simplify = FALSE)
  
  res <- lapply(combs, function(x) {
    keep <- grp %in% x
    dsub <- as.dist(as.matrix(dist_obj)[keep, keep])
    msub <- meta_df[keep, , drop = FALSE]
    msub$GroupTmp <- droplevels(factor(msub[[factor_name]], levels = x))
    
    a <- adonis2(dsub ~ GroupTmp, data = msub, permutations = permutations)
    
    tibble(
      Factor = factor_name,
      Group1 = x[1],
      Group2 = x[2],
      Df = a$Df[1],
      SumOfSqs = a$SumOfSqs[1],
      R2 = a$R2[1],
      F = a$F[1],
      p = a$`Pr(>F)`[1]
    )
  }) %>% bind_rows()
  
  res %>% dplyr::mutate(p_adj_BH = p.adjust(p, method = "BH"))
}

pw_all_bray <- bind_rows(
  pairwise_adonis2(d_bray,  meta2, "Treatment"),
  pairwise_adonis2(d_bray,  meta2, "Water"),
  pairwise_adonis2(d_bray,  meta2, "Endophyte"),
  pairwise_adonis2(d_bray,  meta2, "Harvest")
)

pw_all_aitch <- bind_rows(
  pairwise_adonis2(d_aitch, meta2, "Treatment"),
  pairwise_adonis2(d_aitch, meta2, "Water"),
  pairwise_adonis2(d_aitch, meta2, "Endophyte"),
  pairwise_adonis2(d_aitch, meta2, "Harvest")
)

write.table(
  pw_all_bray,
  file = file.path(perm_dir, "PERMANOVA_pairwise_by_factor_Bray.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  pw_all_aitch,
  file = file.path(perm_dir, "PERMANOVA_pairwise_by_factor_Aitchison.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ==========================================================
# Dispersion tests for all factors
# ==========================================================
run_betadisper <- function(dist_obj, meta_df, factor_name, label) {
  bd <- betadisper(dist_obj, meta_df[[factor_name]])
  bd_aov  <- anova(bd)
  bd_perm <- permutest(bd, permutations = 999)
  
  capture.output(bd_aov,  file = file.path(perm_dir, paste0("Betadisper_", factor_name, "_", label, "_ANOVA.txt")))
  capture.output(bd_perm, file = file.path(perm_dir, paste0("Betadisper_", factor_name, "_", label, "_permutest.txt")))
  
  invisible(list(anova = bd_aov, perm = bd_perm))
}

for (fac in c("Treatment", "Water", "Endophyte", "Harvest")) {
  run_betadisper(d_bray,  meta2, fac, "Bray")
  run_betadisper(d_aitch, meta2, fac, "Aitchison")
}

# ==========================================================
# PCoA
# ==========================================================
make_pcoa_plot <- function(dist_obj, meta_df, label) {
  pcoa <- ape::pcoa(as.matrix(dist_obj))
  
  scores_pcoa <- as.data.frame(pcoa$vectors[, 1:2])
  scores_pcoa$SampleID <- rownames(scores_pcoa)
  
  scores_pcoa <- scores_pcoa %>%
    dplyr::left_join(meta_df, by = "SampleID") %>%
    dplyr::rename(PCoA1 = Axis.1, PCoA2 = Axis.2)
  
  write.table(
    scores_pcoa,
    file = file.path(ord_dir, paste0("PCoA_scores_", label, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  p <- ggplot(scores_pcoa, aes(PCoA1, PCoA2)) +
    geom_point(aes(color = Treatment, shape = Water), size = 3, alpha = 0.9) +
    facet_wrap(~ Harvest) +
    theme_bw() +
    labs(
      title = paste0("PCoA (", label, ") - HUMAnN4 Pathway Matrix"),
      x = paste0("PCoA1 (", round(pcoa$values$Relative_eig[1] * 100, 1), "%)"),
      y = paste0("PCoA2 (", round(pcoa$values$Relative_eig[2] * 100, 1), "%)")
    )
  
  save_gg(p, paste0("PCoA_", label, "_color_Treatment_shape_Water_facet_Harvest"), out_path = ord_dir, w = 10, h = 6)
  invisible(list(scores = scores_pcoa, plot = p, ord = pcoa))
}

pcoa_bray  <- make_pcoa_plot(d_bray,  meta2, "Bray")
pcoa_aitch <- make_pcoa_plot(d_aitch, meta2, "Aitchison")

# ==========================================================
# NMDS
# ==========================================================
make_nmds_plot <- function(Xmat, meta_df, label, distance = "bray") {
  set.seed(1)
  nmds <- metaMDS(Xmat, distance = distance, k = 2, trymax = 200, autotransform = FALSE, trace = FALSE)
  
  stress <- nmds$stress
  message("NMDS stress (", label, ") = ", round(stress, 4))
  
  scores_nmds <- as.data.frame(scores(nmds, display = "sites"))
  scores_nmds$SampleID <- rownames(scores_nmds)
  scores_nmds <- scores_nmds %>% dplyr::left_join(meta_df, by = "SampleID")
  
  write.table(
    scores_nmds,
    file = file.path(ord_dir, paste0("NMDS_scores_", label, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  
  p <- ggplot(scores_nmds, aes(NMDS1, NMDS2)) +
    geom_point(aes(color = Treatment, shape = Water), size = 3, alpha = 0.9) +
    facet_wrap(~ Harvest) +
    theme_bw() +
    labs(title = paste0("NMDS (", label, ") | stress = ", round(stress, 3)))
  
  save_gg(p, paste0("NMDS_", label, "_color_Treatment_shape_Water_facet_Harvest"), out_path = ord_dir, w = 10, h = 6)
  invisible(list(scores = scores_nmds, plot = p, ord = nmds))
}

nmds_bray <- make_nmds_plot(X_cpm, meta2, "Bray", distance = "bray")

# ==========================================================
# HEATMAPS (Top 30 variable pathways; z-scored)
# ==========================================================
N_top <- 30
vars <- matrixStats::colVars(X_cpm)
top_feats <- names(sort(vars, decreasing = TRUE))[1:min(N_top, length(vars))]

X_top <- X_cpm[, top_feats, drop = FALSE]
X_top_z <- scale(X_top)

make_heatmap_by_harvest_grouped <- function(harvest_level,
                                            width_pdf = 16,
                                            height_pdf = 10,
                                            cellw = 10,
                                            cellh = 10,
                                            fs_row = 7,
                                            fs_col = 9,
                                            cluster_cols = FALSE) {
  stopifnot(harvest_level %in% levels(meta2$Harvest))
  
  keep <- meta2$Harvest == harvest_level
  meta_h <- meta2[keep, , drop = FALSE]
  
  meta_h <- meta_h %>%
    dplyr::arrange(Treatment, Water, Endophyte, SampleID)
  
  mat_h <- X_top_z[meta_h$SampleID, , drop = FALSE]
  mat_plot <- t(mat_h)
  
  anno_col <- meta_h %>%
    dplyr::select(SampleID, Treatment, Endophyte, Water) %>%
    as.data.frame()
  
  rownames(anno_col) <- anno_col$SampleID
  anno_col$SampleID <- NULL
  
  slug <- gsub(" ", "_", harvest_level)
  
  pdf(file.path(heat_dir, paste0("Heatmap_Top", N_top, "_", slug, "_samples_groupedByTreatment.pdf")),
      width = width_pdf, height = height_pdf)
  
  pheatmap::pheatmap(
    mat_plot,
    annotation_col = anno_col,
    cluster_cols = cluster_cols,
    cluster_rows = TRUE,
    gaps_col = c(
      sum(meta_h$Treatment == "Control"),
      sum(meta_h$Treatment %in% c("Control", "Polymer"))
    ),
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize_col = fs_col,
    fontsize_row = fs_row,
    cellwidth = cellw,
    cellheight = cellh,
    main = paste0("Top ", N_top, " variable pathways (z-scored) — ", harvest_level)
  )
  dev.off()
  
  png(file.path(heat_dir, paste0("Heatmap_Top", N_top, "_", slug, "_samples_groupedByTreatment.png")),
      width = 5200, height = 3200, res = 300)
  
  pheatmap::pheatmap(
    mat_plot,
    annotation_col = anno_col,
    cluster_cols = cluster_cols,
    cluster_rows = TRUE,
    gaps_col = c(
      sum(meta_h$Treatment == "Control"),
      sum(meta_h$Treatment %in% c("Control", "Polymer"))
    ),
    show_colnames = TRUE,
    show_rownames = TRUE,
    fontsize_col = fs_col,
    fontsize_row = fs_row,
    cellwidth = cellw,
    cellheight = cellh,
    main = paste0("Top ", N_top, " variable pathways (z-scored) — ", harvest_level)
  )
  dev.off()
  
  invisible(TRUE)
}

make_heatmap_by_harvest_grouped("During Drought",
                                width_pdf = 16, height_pdf = 10,
                                cellw = 10, cellh = 10,
                                fs_row = 7, fs_col = 9,
                                cluster_cols = FALSE)

make_heatmap_by_harvest_grouped("Post Drought",
                                width_pdf = 16, height_pdf = 10,
                                cellw = 10, cellh = 10,
                                fs_row = 7, fs_col = 9,
                                cluster_cols = FALSE)

# ======================================================================
# DIFFERENTIAL PATHWAY MODELLING WITH MaAsLin2
# ======================================================================

path_raw <- read.table(
  pathway_raw_fp,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  quote = "",
  comment.char = ""
)

drop_rows_raw <- intersect(
  rownames(path_raw),
  c("UNMAPPED", "UNINTEGRATED", "UNINTEGRATED|unclassified")
)

if (length(drop_rows_raw) > 0) {
  message("Dropping rows from raw table: ", paste(drop_rows_raw, collapse = ", "))
  path_raw <- path_raw[setdiff(rownames(path_raw), drop_rows_raw), , drop = FALSE]
}

meta_maaslin <- meta %>%
  dplyr::filter(SampleID %in% colnames(path_raw)) %>%
  dplyr::mutate(
    Endophyte = factor(Endophyte, levels = c("Nil", "AR37")),
    Treatment = factor(Treatment, levels = c("Control", "Polymer", "6-PP")),
    Water     = factor(Water, levels = c("Irrigated", "Drought")),
    Harvest   = factor(Harvest, levels = c("During Drought", "Post Drought"))
  )

path_raw2 <- path_raw[, meta_maaslin$SampleID, drop = FALSE]
stopifnot(identical(colnames(path_raw2), meta_maaslin$SampleID))

X_raw <- t(as.matrix(path_raw2))
X_raw <- as.data.frame(X_raw)
X_raw$SampleID <- rownames(X_raw)
rownames(X_raw) <- rownames(X_raw)

feat_keep <- colSums(X_raw[, setdiff(colnames(X_raw), "SampleID"), drop = FALSE], na.rm = TRUE) > 0
X_raw_use <- X_raw[, c("SampleID", names(feat_keep)[feat_keep]), drop = FALSE]

meta_maaslin_df <- meta_maaslin %>%
  dplyr::select(SampleID, Endophyte, Treatment, Water, Harvest) %>%
  as.data.frame()

rownames(meta_maaslin_df) <- meta_maaslin_df$SampleID
meta_maaslin_df$SampleID <- NULL

rownames(X_raw_use) <- X_raw_use$SampleID
X_raw_use$SampleID <- NULL

X_raw_use <- X_raw_use[rownames(meta_maaslin_df), , drop = FALSE]
stopifnot(identical(rownames(X_raw_use), rownames(meta_maaslin_df)))

message("MaAsLin2 samples: ", nrow(X_raw_use))
message("MaAsLin2 pathways: ", ncol(X_raw_use))

run_maaslin <- function(input_data,
                        input_metadata,
                        fixed_effects,
                        random_effects = NULL,
                        reference = NULL,
                        analysis_method = "LM",
                        normalization = "TSS",
                        transform = "LOG",
                        standardize = FALSE,
                        max_significance = 0.25,
                        out_subdir = "Maaslin2_run") {
  
  out_path <- file.path(maaslin_dir, out_subdir)
  dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
  
  fit <- Maaslin2(
    input_data = input_data,
    input_metadata = input_metadata,
    output = out_path,
    fixed_effects = fixed_effects,
    random_effects = random_effects,
    reference = reference,
    normalization = normalization,
    transform = transform,
    analysis_method = analysis_method,
    standardize = standardize,
    max_significance = max_significance,
    plot_heatmap = FALSE,
    plot_scatter = FALSE
  )
  
  return(fit)
}

# ----------------------------------------------------------------------
# 1. GLOBAL MODEL
# ----------------------------------------------------------------------
fit_global <- run_maaslin(
  input_data = X_raw_use,
  input_metadata = meta_maaslin_df,
  fixed_effects = c("Endophyte", "Treatment", "Water", "Harvest"),
  reference = c(
    "Endophyte,Nil",
    "Treatment,Control",
    "Water,Irrigated",
    "Harvest,During Drought"
  ),
  normalization = "TSS",
  transform = "LOG",
  out_subdir = "Maaslin2_global"
)

# ----------------------------------------------------------------------
# 2. PAIRWISE CONTRASTS
# ----------------------------------------------------------------------
run_maaslin_pairwise <- function(factor_name, level_a, level_b, prefix) {
  
  keep <- meta_maaslin_df[[factor_name]] %in% c(level_a, level_b)
  
  meta_sub <- meta_maaslin_df[keep, , drop = FALSE]
  data_sub <- X_raw_use[rownames(meta_sub), , drop = FALSE]
  
  meta_sub[[factor_name]] <- factor(meta_sub[[factor_name]], levels = c(level_a, level_b))
  
  ref_string <- paste0(factor_name, ",", level_a)
  
  fit <- run_maaslin(
    input_data = data_sub,
    input_metadata = meta_sub,
    fixed_effects = factor_name,
    reference = ref_string,
    normalization = "TSS",
    transform = "LOG",
    out_subdir = paste0("Maaslin2_", prefix)
  )
  
  invisible(fit)
}

fit_treatment_control_polymer <- run_maaslin_pairwise(
  "Treatment", "Control", "Polymer", "Treatment_Control_vs_Polymer"
)

fit_treatment_control_6pp <- run_maaslin_pairwise(
  "Treatment", "Control", "6-PP", "Treatment_Control_vs_6PP_plus_polymer"
)

fit_treatment_polymer_6pp <- run_maaslin_pairwise(
  "Treatment", "Polymer", "6-PP", "Treatment_Polymer_vs_6PP_plus_polymer"
)

fit_water <- run_maaslin_pairwise("Water", "Irrigated", "Drought", "Water_Irrigated_vs_Drought")
fit_endo <- run_maaslin_pairwise("Endophyte", "Nil", "AR37", "Endophyte_Nil_vs_AR37")
fit_harvest <- run_maaslin_pairwise("Harvest", "During Drought", "Post Drought", "Harvest_During_vs_Post")

read_maaslin_results <- function(results_fp, label) {
  res <- read.delim(results_fp, stringsAsFactors = FALSE)
  res$Comparison <- label
  res
}

volcano_maaslin <- function(res, title, filename_base, q_cut = 0.25, label_top = 10) {
  
  if (!all(c("feature", "coef", "qval") %in% names(res))) {
    stop("Expected columns feature, coef, qval not found in MaAsLin2 results.")
  }
  
  df <- res %>%
    dplyr::mutate(
      neglog10q = -log10(pmax(qval, 1e-300)),
      sig = qval <= q_cut
    ) %>%
    dplyr::arrange(qval, dplyr::desc(abs(coef)))
  
  lab_df <- df %>%
    dplyr::filter(sig) %>%
    dplyr::slice_head(n = label_top)
  
  p <- ggplot(df, aes(x = coef, y = neglog10q)) +
    geom_hline(yintercept = -log10(q_cut), linetype = "dashed", linewidth = 0.35) +
    geom_point(aes(color = sig), alpha = 0.8, size = 2) +
    ggrepel::geom_text_repel(
      data = lab_df,
      aes(label = feature),
      size = 2.7,
      max.overlaps = 30,
      box.padding = 0.35,
      point.padding = 0.2,
      segment.color = "grey50"
    ) +
    scale_color_manual(values = c("FALSE" = "grey70", "TRUE" = "black")) +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = title,
      x = "MaAsLin2 coefficient",
      y = expression(-log[10](qvalue))
    )
  
  save_gg(p, filename_base, out_path = maaslin_dir, w = 8, h = 7)
}

res_treatment_control_polymer <- read_maaslin_results(
  file.path(maaslin_dir, "Maaslin2_Treatment_Control_vs_Polymer", "all_results.tsv"),
  "Control vs Polymer"
)

res_treatment_control_6pp <- read_maaslin_results(
  file.path(maaslin_dir, "Maaslin2_Treatment_Control_vs_6PP_plus_polymer", "all_results.tsv"),
  "Control vs 6-PP"
)

res_treatment_polymer_6pp <- read_maaslin_results(
  file.path(maaslin_dir, "Maaslin2_Treatment_Polymer_vs_6PP_plus_polymer", "all_results.tsv"),
  "Polymer vs 6-PP"
)

res_water <- read_maaslin_results(
  file.path(maaslin_dir, "Maaslin2_Water_Irrigated_vs_Drought", "all_results.tsv"),
  "Irrigated vs Drought"
)

res_endo <- read_maaslin_results(
  file.path(maaslin_dir, "Maaslin2_Endophyte_Nil_vs_AR37", "all_results.tsv"),
  "Nil vs AR37"
)

res_harvest <- read_maaslin_results(
  file.path(maaslin_dir, "Maaslin2_Harvest_During_vs_Post", "all_results.tsv"),
  "During vs Post"
)

volcano_maaslin(res_treatment_control_polymer, "MaAsLin2: Control vs Polymer", "Volcano_MaAsLin2_Control_vs_Polymer")
volcano_maaslin(res_treatment_control_6pp, "MaAsLin2: Control vs 6-PP", "Volcano_MaAsLin2_Control_vs_6PP_plus_polymer")
volcano_maaslin(res_treatment_polymer_6pp, "MaAsLin2: Polymer vs 6-PP", "Volcano_MaAsLin2_Polymer_vs_6PP_plus_polymer")
volcano_maaslin(res_water, "MaAsLin2: Irrigated vs Drought", "Volcano_MaAsLin2_Water")
volcano_maaslin(res_endo, "MaAsLin2: Nil vs AR37", "Volcano_MaAsLin2_Endophyte")
volcano_maaslin(res_harvest, "MaAsLin2: During vs Post", "Volcano_MaAsLin2_Harvest")

# ==========================================================
# Bubble plots from MaAsLin2 results
# ==========================================================
make_bubble_maaslin <- function(res_list,
                                factor_label,
                                out_prefix,
                                top_n = 40,
                                q_cut = 0.25) {
  
  bubble_input <- dplyr::bind_rows(res_list) %>%
    dplyr::filter(metadata != "(Intercept)") %>%
    dplyr::mutate(
      neglog10 = -log10(pmax(qval, 1e-300)),
      abs_coef = abs(coef)
    )
  
  top_show <- bubble_input %>%
    dplyr::group_by(feature) %>%
    dplyr::summarise(
      max_score = max(neglog10, na.rm = TRUE),
      max_abs_coef = max(abs_coef, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(max_score), dplyr::desc(max_abs_coef)) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::pull(feature)
  
  bubble_df <- bubble_input %>%
    dplyr::filter(feature %in% top_show) %>%
    dplyr::mutate(
      Pathway = forcats::fct_reorder(feature, coef, .fun = median, .desc = TRUE),
      Significant = qval <= q_cut
    )
  
  p <- ggplot(bubble_df, aes(x = Comparison, y = Pathway)) +
    geom_point(aes(size = neglog10, color = coef), alpha = 0.85) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      title = paste0("Differential pathways across ", factor_label, " (MaAsLin2)"),
      subtitle = "Color = MaAsLin2 coefficient; Size = -log10(q-value)",
      x = NULL,
      y = NULL,
      color = "Coefficient",
      size = expression(-log[10](q))
    )
  
  save_gg(p, out_prefix, out_path = bubble_dir, w = 12, h = 10)
  
  ranked <- bubble_input %>%
    dplyr::arrange(dplyr::desc(neglog10), dplyr::desc(abs_coef))
  
  write.table(
    ranked,
    file = file.path(bubble_dir, paste0(out_prefix, "_RANKED.tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  invisible(bubble_df)
}

bubble_treatment <- make_bubble_maaslin(
  res_list = list(
    res_treatment_control_polymer,
    res_treatment_control_6pp,
    res_treatment_polymer_6pp
  ),
  factor_label = "Treatment",
  out_prefix = "Bubble_MaAsLin2_Treatment_top40",
  top_n = 40,
  q_cut = 0.25
)

bubble_water <- make_bubble_maaslin(
  res_list = list(res_water),
  factor_label = "Water",
  out_prefix = "Bubble_MaAsLin2_Water_top40",
  top_n = 40,
  q_cut = 0.25
)

bubble_endo <- make_bubble_maaslin(
  res_list = list(res_endo),
  factor_label = "Endophyte",
  out_prefix = "Bubble_MaAsLin2_Endophyte_top40",
  top_n = 40,
  q_cut = 0.25
)

bubble_harvest <- make_bubble_maaslin(
  res_list = list(res_harvest),
  factor_label = "Harvest",
  out_prefix = "Bubble_MaAsLin2_Harvest_top40",
  top_n = 40,
  q_cut = 0.25
)

# ==========================================================
# Focused pathway dot plot
# ==========================================================
focused_patterns <- c(
  "glycolysis",
  "fermentation",
  "TCA cycle",
  "tricarboxylic",
  "respiration",
  "oxidative phosphorylation",
  "L-glutamine biosynthesis",
  "nucleotide",
  "UTP and CTP dephosphorylation"
)

focused_pathways <- unique(unlist(lapply(focused_patterns, function(p) {
  grep(p, colnames(X_cpm), value = TRUE, ignore.case = TRUE)
})))

focused_pathways <- intersect(focused_pathways, colnames(X_cpm))

if (length(focused_pathways) > 12) {
  vars_focus <- apply(X_cpm[, focused_pathways, drop = FALSE], 2, var, na.rm = TRUE)
  focused_pathways <- names(sort(vars_focus, decreasing = TRUE))[1:12]
}

message("Focused pathways selected:")
print(focused_pathways)

focus_df <- as.data.frame(X_cpm[, focused_pathways, drop = FALSE]) %>%
  tibble::rownames_to_column("SampleID") %>%
  tidyr::pivot_longer(-SampleID, names_to = "Pathway", values_to = "CPM") %>%
  dplyr::left_join(meta2, by = "SampleID")

focus_df <- focus_df %>%
  dplyr::group_by(Pathway) %>%
  dplyr::mutate(PathwayZ = as.numeric(scale(CPM))) %>%
  dplyr::ungroup()

focus_sum <- focus_df %>%
  dplyr::group_by(Harvest, Water, Treatment, Pathway) %>%
  dplyr::summarise(
    PathwayZ = mean(PathwayZ, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Pathway = factor(Pathway, levels = rev(focused_pathways)),
    Treatment = factor(Treatment, levels = c("Control", "Polymer", "6-PP")),
    Water = factor(Water, levels = c("Irrigated", "Drought")),
    Harvest = factor(Harvest, levels = c("During Drought", "Post Drought"))
  )

p_focus <- ggplot(focus_sum, aes(x = Treatment, y = Pathway)) +
  geom_point(aes(size = abs(PathwayZ), color = PathwayZ), alpha = 0.9) +
  facet_grid(Water ~ Harvest) +
  scale_color_gradient2(
    low = "#2166AC",
    mid = "grey95",
    high = "#B2182B",
    midpoint = 0
  ) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 35, hjust = 1),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  ) +
  labs(
    title = "Focused energy/stress pathways",
    subtitle = "Dot size = |z|, color = relative pathway activity",
    x = NULL,
    y = NULL,
    color = "Pathway Z",
    size = "abs(Pathway Z)"
  )

save_gg(p_focus, "Fig_FocusedPathways_Dotplot_Treatment_facet_Water_Harvest", out_path = focus_dir, w = 12, h = 8)

sig_bind <- dplyr::bind_rows(
  res_treatment_control_polymer,
  res_treatment_control_6pp,
  res_treatment_polymer_6pp,
  res_water,
  res_endo,
  res_harvest
) %>%
  dplyr::filter(qval <= 0.25) %>%
  dplyr::arrange(Comparison, qval, dplyr::desc(abs(coef)))

write.table(
  sig_bind,
  file = file.path(maaslin_dir, "Maaslin2_significant_all_contrasts.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ==========================================================
# Pathway class scoring for mechanistic interpretation
# ==========================================================
zscore_cols <- function(M) {
  M <- as.matrix(M)
  Z <- scale(M)
  rownames(Z) <- rownames(M)
  colnames(Z) <- colnames(M)
  Z
}

df_long <- as.data.frame(X_cpm) %>%
  tibble::rownames_to_column("SampleID") %>%
  tidyr::pivot_longer(-SampleID, names_to = "Pathway", values_to = "CPM") %>%
  dplyr::left_join(meta2, by = "SampleID")

col_raw <- colnames(X_cpm)

pick_cols <- function(pattern) {
  grep(pattern, col_raw, value = TRUE, ignore.case = TRUE)
}

pathway_sets_mapped <- list(
  Glycolysis   = pick_cols("\\bglycolysis\\b"),
  Fermentation = pick_cols("fermentation|homolactic|mixed[- ]acid|lactate|ethanol"),
  TCA          = pick_cols("TCA cycle|tricarboxylic"),
  Resp_Redox   = pick_cols("respiration|electron transport|oxidative phosphorylation|NADH|FADH2|ubiquinone|menaquinone|oxidation|dehydrogenase"),
  Nucleotide   = pick_cols("purine|pyrimidine|nucleotide|inosine|adenosine|guanosine|uracil|thymid"),
  AminoAcid    = pick_cols("amino acid biosynthesis|lysine biosynthesis|methionine biosynthesis|arginine biosynthesis|tryptophan biosynthesis|histidine biosynthesis|branched-chain amino acid|isoleucine biosynthesis|valine biosynthesis|leucine biosynthesis")
)

pathway_sets_mapped <- pathway_sets_mapped[sapply(pathway_sets_mapped, length) > 0]
message("Mapped pathway-set sizes:")
print(sapply(pathway_sets_mapped, length))

if (length(pathway_sets_mapped) == 0) {
  stop("All pathway classes are empty after mapping. Check regex patterns vs colnames(X_cpm).")
}

X_z <- zscore_cols(X_cpm)

class_scores <- purrr::imap_dfr(pathway_sets_mapped, function(pws_cols, set_name) {
  pws_cols <- intersect(pws_cols, colnames(X_z))
  if (length(pws_cols) == 0) return(NULL)
  
  tibble::tibble(
    SampleID = rownames(X_z),
    Class = set_name,
    Score = rowMeans(X_z[, pws_cols, drop = FALSE], na.rm = TRUE),
    N_pathways_in_class = length(pws_cols)
  )
}) %>%
  dplyr::left_join(meta2, by = "SampleID")

write.table(
  class_scores,
  file = file.path(class_dir, "PathwayClassScores_perSample.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

signature_pathways <- unique(unlist(pathway_sets_mapped, use.names = FALSE))
signature_pathways <- intersect(signature_pathways, colnames(X_cpm))

write.table(
  data.frame(Pathway = signature_pathways),
  file = file.path(class_dir, "Signature_pathways_used_for_class_scoring.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ==========================================================
# Reproducibility log
# ==========================================================
sink(file.path(repr_dir, "sessionInfo.txt"))
print(sessionInfo())
sink()

message("DONE. Outputs saved to: ", out_dir)

# ==========================================================
# Figure 4 — Nature Communications microbial function figure
# A = HUMAnN4 Bray-Curtis PCoA
# B = PERMANOVA R2 summary
# C = Focused energy/stress pathway activity
# D = MaAsLin2 treatment-associated pathways
# ==========================================================

suppressPackageStartupMessages({
  library(cowplot)
  library(forcats)
})

fig4_dir <- file.path(project_dir, "7.manuscript_figure", "5.figure4")
dir.create(fig4_dir, recursive = TRUE, showWarnings = FALSE)

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

# Panel A — HUMAnN4 Bray-Curtis PCoA
pcoa_scores <- pcoa_bray$scores
pcoa_var1 <- round(100 * pcoa_bray$ord$values$Relative_eig[1], 1)
pcoa_var2 <- round(100 * pcoa_bray$ord$values$Relative_eig[2], 1)

panel_a <- ggplot(
  pcoa_scores,
  aes(x = PCoA1, y = PCoA2, colour = Treatment, shape = Water)
) +
  geom_point(size = 3.8, alpha = 0.9) +
  facet_wrap(~ Harvest, nrow = 1) +
  scale_colour_manual(values = c(
    "Control" = "#F8766D",
    "Polymer" = "#00BA38",
    "6-PP" = "#619CFF"
  )) +
  labs(
    title = "A",
    subtitle = "HUMAnN4 pathway Bray-Curtis PCoA",
    x = paste0("PCoA1 (", pcoa_var1, "%)"),
    y = paste0("PCoA2 (", pcoa_var2, "%)"),
    colour = "Seed treatment",
    shape = "Water"
  ) +
  theme_natcom(15)

# Panel B — PERMANOVA variance explained including main effects, interactions, and residual

sig_code <- function(p) {
  dplyr::case_when(
    is.na(p)    ~ "",
    p <= 0.001 ~ "***",
    p <= 0.01  ~ "**",
    p <= 0.05  ~ "*",
    TRUE       ~ ""
  )
}

clean_perm_term <- function(x) {
  dplyr::case_when(
    x == "Endophyte" ~ "Endophyte",
    x == "Treatment" ~ "Seed treatment",
    x == "Water" ~ "Water regime",
    x == "Harvest" ~ "Harvest stage",
    x == "Treatment:Water" ~ "Seed treatment × Water regime",
    x == "Treatment:Harvest" ~ "Seed treatment × Harvest stage",
    x == "Water:Harvest" ~ "Water regime × Harvest stage",
    x == "Residual" ~ "Unexplained variance",
    TRUE ~ x
  )
}

perm_main_df <- as.data.frame(perm_bray$main) %>%
  tibble::rownames_to_column("Term")

perm_2way_df <- as.data.frame(perm_bray$twoway) %>%
  tibble::rownames_to_column("Term")

perm_function_df <- dplyr::bind_rows(
  perm_main_df %>%
    dplyr::filter(Term %in% c("Endophyte", "Treatment", "Water", "Harvest")),
  
  perm_2way_df %>%
    dplyr::filter(Term %in% c(
      "Treatment:Water",
      "Treatment:Harvest",
      "Water:Harvest",
      "Residual"
    ))
) %>%
  dplyr::rename(p_value = `Pr(>F)`) %>%
  dplyr::mutate(
    Term_clean = clean_perm_term(Term),
    sig = sig_code(p_value),
    significant = dplyr::case_when(
      Term == "Residual" ~ "Unexplained variance",
      p_value < 0.05 ~ "Significant",
      TRUE ~ "Not significant"
    ),
    label = dplyr::case_when(
      Term == "Residual" ~ sprintf("%.3f", R2),
      sig == "" ~ sprintf("%.3f", R2),
      TRUE ~ sprintf("%.3f %s", R2, sig)
    ),
    Term_clean = factor(
      Term_clean,
      levels = rev(c(
        "Harvest stage",
        "Water regime",
        "Seed treatment",
        "Endophyte",
        "Water regime × Harvest stage",
        "Seed treatment × Water regime",
        "Seed treatment × Harvest stage",
        "Unexplained variance"
      ))
    )
  )

panel_b <- ggplot(
  perm_function_df,
  aes(x = R2, y = Term_clean, fill = significant)
) +
  geom_col(width = 0.72) +
  geom_text(aes(label = label), hjust = -0.12, size = 3.8) +
  scale_fill_manual(
    values = c(
      "Significant" = "#8A3FFC",
      "Not significant" = "grey70",
      "Unexplained variance" = "grey35"
    ),
    drop = FALSE,
    name = NULL
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "B",
    subtitle = "Variance explained by PERMANOVA model terms",
    x = expression(R^2~"(variance explained)"),
    y = NULL
  ) +
  theme_natcom(14) +
  theme(
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.title.x = element_text(size = 14, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 10),
    plot.subtitle = element_text(size = 14, face = "bold"),
    plot.margin = margin(8, 30, 8, 8)
  )

# Panel C — Focused energy/stress pathway activity

shorten_pathway <- function(x, wrap_width = 32) {
  x %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_replace_all("\\.", " ") %>%
    
    # aggressive fixes first
    stringr::str_replace_all(".*fatty acid beta oxidation V.*", "FA beta-oxidation V") %>%
    stringr::str_replace_all(".*fatty acid beta oxidation IV.*", "FA beta-oxidation IV") %>%
    stringr::str_replace_all(".*fatty acid beta oxidation I.*", "FA beta-oxidation I") %>%
    stringr::str_replace_all(".*409.*purine.*salvage.*", "purine nt. salvage") %>%
    stringr::str_replace_all(".*REDCITCYC.*TCA cycle VI.*", "TCA cycle VI") %>%
    
    # remove pathway IDs/prefixes
    stringr::str_replace_all("^PWY[- ]?[0-9]*[: ]*", "") %>%
    stringr::str_replace_all("^[A-Z0-9-]+[- ]?PWY[: ]*", "") %>%
    stringr::str_replace_all("^P[0-9]+[- ]?PWY[: ]*", "") %>%
    stringr::str_replace_all("^P[0-9]+[: ]*", "") %>%
    stringr::str_replace_all("^TCA[−-]GLYOX[−-]BYPASS[: ]*", "TCA-glyoxylate bypass ") %>%
    stringr::str_replace_all("^RED[C]?T[C]?YC[: ]*", "") %>%
    stringr::str_replace_all("^COMPLETE[- ]ARO[- ]PWY[: ]*", "") %>%
    stringr::str_replace_all("^GLYCOLYSIS[-: ]*", "glycolysis ") %>%
    stringr::str_replace_all("^FERMENTATION[-: ]*", "") %>%
    stringr::str_replace_all("^ANAGLYCOLYSIS[-: ]*", "anaerobic glycolysis ") %>%
    stringr::str_replace_all("^CENTFERM[-: ]*", "") %>%
    stringr::str_replace_all("^FAO[-: ]*", "FA beta-oxidation ") %>%
    stringr::str_replace_all("^ARGDEG[- ]PWY[: ]*", "") %>%
    stringr::str_replace_all("^FASYN[- ]INITIAL[- ]PWY[: ]*", "") %>%
    stringr::str_replace_all("^BRANCHED[- ]CHAIN[- ]AA[- ]SYN[- ]PWY[: ]*", "") %>%
    
    # biological shortening
    stringr::str_replace_all("mixed acid ferment.*", "mixed-acid fermentation") %>%
    stringr::str_replace_all("homolactic ferment.*", "homolactic fermentation") %>%
    stringr::str_replace_all("pyruvate ferment.*butan.*", "pyruvate fermentation") %>%
    stringr::str_replace_all("pyruvate ferment.*propanoate.*", "pyruvate to propanoate") %>%
    stringr::str_replace_all("glycolysis.*Entner.*", "Entner-Doudoroff glycolysis") %>%
    stringr::str_replace_all("glycolysis.*glucose.*", "glycolysis") %>%
    stringr::str_replace_all("fatty acid biosyn.*", "fatty acid biosyn.") %>%
    stringr::str_replace_all("branched chain AA biosyn.*", "branched-chain AA biosyn.") %>%
    stringr::str_replace_all("L arginine.*putrescine.*", "arginine/putrescine/GABA") %>%
    stringr::str_replace_all("aromatic AA biosyn.*", "aromatic AA biosyn.") %>%
    stringr::str_replace_all("1D myo inositol hexakisphosphate.*", "inositol phosphate biosyn.") %>%
    stringr::str_replace_all("4 hydroxybenzoate biosyn.*", "4-hydroxybenzoate biosyn.") %>%
    stringr::str_replace_all("C4 photosynthetic carbon assimilation.*", "C4 photosyn. assim.") %>%
    stringr::str_replace_all("L ascorbate biosyn.*", "ascorbate biosyn.") %>%
    stringr::str_replace_all("taxadiene biosyn.*", "taxadiene biosyn.") %>%
    stringr::str_replace_all("vitamin E biosyn.*tocotrienols.*", "vitamin E biosyn. tocotrienols") %>%
    stringr::str_replace_all("vitamin E biosyn.*tocopherols.*", "vitamin E biosyn. tocopherols") %>%
    stringr::str_replace_all("vitamin E biosyn.*", "vitamin E biosyn.") %>%
    stringr::str_replace_all("purine nucleotides degrad.*", "purine nt. degrad.") %>%
    stringr::str_replace_all("adenine and adenosine salvage.*", "adenosine salvage") %>%
    stringr::str_replace_all("TCA cycle VII acetate producers", "TCA cycle acetate producers") %>%
    stringr::str_replace_all("sulfur amino acid biosyn.*", "sulfur AA biosyn.") %>%
    
    # general shortening
    stringr::str_replace_all("superpathway of ", "") %>%
    stringr::str_replace_all("biosynthesis", "biosyn.") %>%
    stringr::str_replace_all("degradation", "degrad.") %>%
    stringr::str_replace_all("fermentation", "ferment.") %>%
    stringr::str_replace_all("tricarboxylic acid", "TCA") %>%
    stringr::str_replace_all("amino acid", "AA") %>%
    stringr::str_replace_all("nucleotide", "nt.") %>%
    stringr::str_replace_all("photosynthetic", "photosyn.") %>%
    stringr::str_replace_all("assimilation", "assim.") %>%
    
    # remove species/detail
    stringr::str_replace_all("\\(.*?\\)", "") %>%
    stringr::str_replace_all("\\bE coli\\b", "") %>%
    stringr::str_replace_all("\\bplants\\b", "") %>%
    stringr::str_replace_all("\\bSpirodela polyrhiza\\b", "") %>%
    stringr::str_replace_all("\\bSaccharomyces cerevisiae\\b", "") %>%
    stringr::str_replace_all("\\bHelicobacter\\b", "") %>%
    stringr::str_replace_all("\\s+I{1,3}$", "") %>%
    stringr::str_replace_all("\\s+IV$", "") %>%
    stringr::str_replace_all("\\s+V$", "") %>%
    stringr::str_replace_all("\\s+VI$", "") %>%
    stringr::str_replace_all("\\s+VII$", "") %>%
    
    stringr::str_squish() %>%
    stringr::str_wrap(width = wrap_width)
}

focus_label <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "FERMENTATION-PWY") ~ "Mixed-acid fermentation",
    stringr::str_detect(x, "PWY-5913") ~ "Partial TCA cycle",
    stringr::str_detect(x, "PWY-6549") ~ "L-glutamine biosyn. III",
    stringr::str_detect(x, "GLYCOLYSIS: glycolysis I") ~ "Glycolysis I",
    stringr::str_detect(x, "PWY-7185") ~ "UTP/CTP dephosphorylation I",
    stringr::str_detect(x, "PWY-3781") ~ "Aerobic respiration I",
    stringr::str_detect(x, "PWY-5484") ~ "Glycolysis II",
    stringr::str_detect(x, "PWY-5690") ~ "TCA cycle II",
    stringr::str_detect(x, "PWY-1042") ~ "Glycolysis IV",
    stringr::str_detect(x, "ANAGLYCOLYSIS-PWY") ~ "Glycolysis III",
    stringr::str_detect(x, "PWY-7210") ~ "Pyrimidine dNTP biosyn. from CTP",
    stringr::str_detect(x, "PWY0-166") ~ "De novo pyrimidine dNTP biosyn.",
    TRUE ~ shorten_pathway(x, wrap_width = 28)
  )
}

panel_d_label <- function(x) {
  shorten_pathway(x, wrap_width = 30)
}

# Focused pathway identifiers used for manuscript Figure 4 panel C.
# These terms are matched against the HUMAnN4 pathway column names in X_cpm.
focused_pathway_terms <- c(
  "FERMENTATION-PWY",
  "PWY-5913",
  "PWY-6549",
  "GLYCOLYSIS: glycolysis I",
  "PWY-7185",
  "PWY-3781",
  "PWY-5484",
  "PWY-5690",
  "PWY-1042",
  "ANAGLYCOLYSIS-PWY",
  "PWY-7210",
  "PWY0-166"
)

focus_feats <- unlist(lapply(
  focused_pathway_terms,
  function(term) {
    colnames(X_cpm)[stringr::str_detect(colnames(X_cpm), fixed(term))]
  }
))

focus_feats <- unique(focus_feats)

if (length(focus_feats) == 0) {
  stop(
    "No focused pathway features matched focused_pathway_terms. ",
    "Check focused_pathway_terms against colnames(X_cpm)."
  )
}

focus_z <- scale(X_cpm[, focus_feats, drop = FALSE])

focus_long <- as.data.frame(focus_z) %>%
  tibble::rownames_to_column("SampleID") %>%
  tidyr::pivot_longer(
    cols = -SampleID,
    names_to = "Pathway",
    values_to = "Pathway_Z"
  ) %>%
  dplyr::left_join(meta2, by = "SampleID") %>%
  dplyr::group_by(Pathway, Treatment, Water, Harvest) %>%
  dplyr::summarise(
    mean_z = mean(Pathway_Z, na.rm = TRUE),
    abs_z = abs(mean_z),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    Pathway_clean = focus_label(Pathway),
    Pathway_clean = factor(
      Pathway_clean,
      levels = rev(c(
        "Mixed-acid fermentation",
        "Partial TCA cycle",
        "L-glutamine biosyn. III",
        "Glycolysis I",
        "UTP/CTP dephosphorylation I",
        "Aerobic respiration I",
        "Glycolysis II",
        "TCA cycle II",
        "Glycolysis IV",
        "Glycolysis III",
        "Pyrimidine dNTP biosyn. from CTP",
        "De novo pyrimidine dNTP biosyn."
      ))
    ),
    Treatment = factor(Treatment, levels = c("Control", "Polymer", "6-PP")),
    Harvest = factor(Harvest, levels = c("During Drought", "Post Drought")),
    Water = factor(Water, levels = c("Irrigated", "Drought"))
  )

panel_c <- ggplot(
  focus_long,
  aes(x = Treatment, y = Pathway_clean, size = abs_z, colour = mean_z)
) +
  geom_point(alpha = 0.9) +
  facet_grid(Water ~ Harvest) +
  scale_colour_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1.5, 1.5),
    name = "Pathway Z"
  ) +
  scale_size_continuous(
    name = "|Z|",
    range = c(2.8, 8.5),
    breaks = c(0.5, 1.0, 1.5)
  ) +
  labs(
    title = "C",
    subtitle = "Central carbon and stress-response pathways",
    x = "Seed treatment",
    y = NULL
  ) +
  theme_natcom(16) +
  
  theme(
    panel.grid.major.y = element_line(
      colour = "grey88",
      linewidth = 0.35
    ),
    panel.grid.major.x = element_blank(),
    
    axis.text.x = element_text(
      angle = 35,
      hjust = 1,
      size = 13
    ),
    
    axis.text.y = element_text(
      size = 14,
      colour = "black"
    ),
    
    axis.title.x = element_text(
      size = 15,
      face = "bold"
    ),
    
    strip.text = element_text(
      size = 15,
      face = "bold"
    ),
    
    legend.title = element_text(
      size = 15,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 14
    ),
    
    plot.subtitle = element_text(
      size = 14,
      face = "bold"
    ),
    
    legend.position = "right",
    
    plot.margin = margin(
      8, 8, 8, 55
    )
  )

# Panel D — MaAsLin2 treatment-associated pathways

panel_d_label <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, "6837|odd number") ~
      "FA β-oxidation V",
    
    stringr::str_detect(x, "5138|even number") ~
      "FA β-oxidation IV",
    
    stringr::str_detect(x, "FAO.PWY|generic") ~
      "FA β-oxidation I",
    
    stringr::str_detect(x, "4661|inositol") ~
      "Inositol phosphate biosyn.",
    
    stringr::str_detect(x, "6435") ~
      "4-hydroxybenzoate biosyn.",
    
    stringr::str_detect(x, "735") ~
      "L-tyrosine biosyn.",
    
    stringr::str_detect(x, "6628") ~
      "L-phenylalanine biosyn.",
    
    stringr::str_detect(x, "1422") ~
      "Vitamin E biosyn. (tocopherols)",
    
    stringr::str_detect(x, "7436") ~
      "Vitamin E biosyn. (tocotrienols)",
    
    stringr::str_detect(x, "8143") ~
      "Ascorbate biosyn.",
    
    stringr::str_detect(x, "7392") ~
      "Taxadiene biosyn.",
    
    stringr::str_detect(x, "2201") ~
      "Folate transformations",
    
    stringr::str_detect(x, "821") ~
      "Sulfur AA biosyn.",
    
    stringr::str_detect(x, "6609") ~
      "Adenosine salvage",
    
    stringr::str_detect(x, "6353") ~
      "Purine nt. degrad.",
    
    stringr::str_detect(x, "66.409|409") ~
      "Purine nt. salvage",
    
    stringr::str_detect(x, "REDCITCYC") ~
      "TCA cycle VI",
    
    stringr::str_detect(x, "7254") ~
      "TCA cycle VII acetate prod.",
    
    stringr::str_detect(x, "241") ~
      "C4 photosyn. assim.",
    
    TRUE ~ shorten_pathway(x, wrap_width = 26)
  )
}

maaslin_treatment_df <- bind_rows(
  res_treatment_control_polymer,
  res_treatment_control_6pp,
  res_treatment_polymer_6pp
) %>%
  dplyr::filter(metadata != "(Intercept)") %>%
  dplyr::mutate(
    neglog10q = -log10(pmax(qval, 1e-300)),
    abs_coef = abs(coef)
  )

top_pathways <- maaslin_treatment_df %>%
  dplyr::filter(Comparison == "Control vs 6-PP") %>%
  dplyr::arrange(dplyr::desc(neglog10q), dplyr::desc(abs_coef)) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::pull(feature)

bubble_df <- maaslin_treatment_df %>%
  dplyr::filter(feature %in% top_pathways) %>%
  dplyr::mutate(
    Comparison = factor(
      Comparison,
      levels = c("Control vs Polymer", "Control vs 6-PP", "Polymer vs 6-PP")
    ),
    Pathway = panel_d_label(feature),
    Pathway = forcats::fct_reorder(Pathway, neglog10q, .fun = max)
  )

panel_d <- ggplot(
  bubble_df,
  aes(x = Comparison, y = Pathway, size = neglog10q, colour = coef)
) +
  geom_point(alpha = 0.9) +
  scale_colour_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Coefficient"
  ) +
  scale_size_continuous(
    name = expression(-log[10](q)),
    range = c(1.8, 6.5)
  ) +
  labs(
    title = "D",
    subtitle = "Differential metabolic pathways associated with seed treatment",
    x = NULL,
    y = NULL
  ) +
  theme_natcom(15) +
  
  theme(
    
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      size = 13,
      colour = "black"
    ),
    
    axis.text.y = element_text(
      size = 15,
      colour = "black"
    ),
    
    axis.title = element_text(
      size = 15,
      face = "bold"
    ),
    
    strip.text = element_text(
      size = 16,
      face = "bold"
    ),
    
    legend.title = element_text(
      size = 15,
      face = "bold"
    ),
    
    legend.text = element_text(
      size = 14
    ),
    
    plot.subtitle = element_text(
      size = 14,
      face = "bold"
    ),
    
    plot.title = element_text(
      size = 16,
      face = "bold"
    ),
    
    legend.position = "right",
    
    plot.margin = margin(8, 8, 8, 55)
    
  )
top <- plot_grid(panel_a, panel_b, nrow = 1, rel_widths = c(1, 1))
bottom <- plot_grid(panel_c, panel_d, nrow = 1, rel_widths = c(1.15, 1.0))

fig4 <- plot_grid(
  top,
  bottom,
  ncol = 1,
  rel_heights = c(0.85, 1.35)
)

ggsave(
  file.path(fig4_dir, "Figure4_NatCom_function.pdf"),
  fig4,
  width = 24,
  height = 18,
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  file.path(fig4_dir, "Figure4_NatCom_function.png"),
  fig4,
  width = 24,
  height = 18,
  dpi = 600,
  limitsize = FALSE
)

message("Figure 4 saved to: ", fig4_dir)

