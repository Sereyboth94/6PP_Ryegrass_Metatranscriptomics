###############################################################################
# 🧬 Publication-ready RNA-seq Transcriptomics Pipeline
# Ryegrass × Water × Treatment × Harvest
# Author: Sereyboth Soth
#
# Main features:
#   - DESeq2 differential expression analysis
#   - stricter low-count filtering
#   - VST PCA + sample distance heatmap
#   - gene annotation using AnnotationHub
#   - volcano plots with gene names
#   - top variable gene heatmaps using actual gene names
#   - focused heatmaps of significant genes
#   - GO enrichment analysis with clusterProfiler
#   - stratified subset analyses for biologically interpretable contrasts
#   - true global 6-PP vs Control comparison
#   - reproducible Figure S14 (During vs Post drought treatment volcanos)
###############################################################################

# ─────────────────────────────────────────────────────────────────────────────
# 1️⃣ Load libraries
# ─────────────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(matrixStats)
  library(ggplot2)
  library(ggrepel)
  library(AnnotationHub)
  library(AnnotationDbi)
  library(tibble)
  library(stringr)
  library(here)
  library(patchwork)
  library(vegan)
})

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

bioc_pkgs <- c("clusterProfiler")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}
library(clusterProfiler)

# ─────────────────────────────────────────────────────────────────────────────
# 2️⃣ Define project paths and input files
# ─────────────────────────────────────────────────────────────────────────────
project_dir <- "C:/Users/soths/OneDrive - Lincoln University/Writing/NatMicro/github"

data_dir <- file.path(project_dir, "4.plant_transcriptome", "data")
out_dir  <- file.path(project_dir, "4.plant_transcriptome", "results")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

qc_dir   <- file.path(out_dir, "01_QC")
pca_dir  <- file.path(out_dir, "02_PCA_and_SampleDistance")
heat_dir <- file.path(out_dir, "03_Heatmaps")
de_dir   <- file.path(out_dir, "04_DESeq2")
go_dir   <- file.path(out_dir, "05_GO")
fig_dir  <- file.path(out_dir, "06_Publication_Figures")
rep_dir  <- file.path(out_dir, "07_Reproducibility")

dir.create(qc_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(pca_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(heat_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(de_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(go_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(rep_dir,  recursive = TRUE, showWarnings = FALSE)

counts_fp <- file.path(data_dir, "featureCounts_ryegrass_clean.txt")
meta_fp   <- file.path(data_dir, "Metadata.csv")

stopifnot(
  file.exists(counts_fp),
  file.exists(meta_fp)
)

# ─────────────────────────────────────────────────────────────────────────────
# 3️⃣ Helper functions
# ─────────────────────────────────────────────────────────────────────────────
save_plot <- function(plot_obj, file_base, out_path = out_dir, w = 8, h = 6, dpi = 300) {
  ggsave(
    filename = file.path(out_path, paste0(file_base, ".pdf")),
    plot = plot_obj,
    width = w, height = h,
    useDingbats = FALSE,
    limitsize = FALSE
  )
  ggsave(
    filename = file.path(out_path, paste0(file_base, ".png")),
    plot = plot_obj,
    width = w, height = h,
    dpi = dpi,
    limitsize = FALSE
  )
}

theme_pub <- function(base_size = 13) {
  theme_minimal(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold")
    )
}

make_unique_labels <- function(ids, labels) {
  out <- ifelse(is.na(labels) | labels == "", ids, labels)
  make.unique(out)
}

make_simple_volcano <- function(res_df, panel_title) {
  plot_df <- res_df %>%
    dplyr::mutate(
      padj_plot = ifelse(is.na(padj), 1, padj),
      neglog10_padj = -log10(padj_plot),
      significance = dplyr::case_when(
        !is.na(padj) & padj < 0.05 & log2FoldChange > 1  ~ "Upregulated",
        !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
        TRUE ~ "Not Significant"
      )
    )
  
  ggplot(plot_df, aes(x = log2FoldChange, y = neglog10_padj, color = significance)) +
    geom_point(alpha = 0.7, size = 1.8) +
    scale_color_manual(values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "grey70"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(
      title = panel_title,
      x = "log2 fold change",
      y = expression(-log[10]("adjusted p-value")),
      color = "significance"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
}

# ─────────────────────────────────────────────────────────────────────────────
# 4️⃣ Load input data
# ─────────────────────────────────────────────────────────────────────────────
counts_raw <- read.delim(counts_fp, comment.char = "#", check.names = FALSE)
metadata   <- read.csv(meta_fp, header = TRUE, stringsAsFactors = FALSE)

message("✅ Input files loaded successfully")
message("Counts table rows: ", nrow(counts_raw))
message("Metadata rows: ", nrow(metadata))

# ─────────────────────────────────────────────────────────────────────────────
# 5️⃣ Load and clean featureCounts table
# ─────────────────────────────────────────────────────────────────────────────
count_data <- counts_raw %>%
  dplyr::select(-Chr, -Start, -End, -Strand, -Length)

rownames(count_data) <- count_data$Geneid
count_data <- count_data %>% dplyr::select(-Geneid)

count_mat <- as.matrix(count_data)
storage.mode(count_mat) <- "integer"

message("✅ Count matrix dimensions: ", nrow(count_mat), " genes × ", ncol(count_mat), " samples")

# ─────────────────────────────────────────────────────────────────────────────
# 6️⃣ Load and validate metadata
# ─────────────────────────────────────────────────────────────────────────────
required_cols <- c("Sample", "Ryegrass", "Water", "Treatment", "Harvest")
missing_cols <- setdiff(required_cols, names(metadata))
if (length(missing_cols) > 0) {
  stop("❌ Missing columns in metadata: ", paste(missing_cols, collapse = ", "))
}

if (!all(colnames(count_mat) %in% metadata$Sample)) {
  missing_in_meta <- setdiff(colnames(count_mat), metadata$Sample)
  stop("❌ These count columns are not found in metadata: ", paste(missing_in_meta, collapse = ", "))
}

metadata <- metadata[match(colnames(count_mat), metadata$Sample), , drop = FALSE]
stopifnot(all(colnames(count_mat) == metadata$Sample))

metadata <- metadata %>%
  dplyr::mutate(
    Sample = as.character(Sample),
    Treatment = dplyr::recode(as.character(Treatment), "VOC" = "6-PP"),
    Harvest = dplyr::recode(
      as.character(Harvest),
      "During_Drought" = "During Drought",
      "Post_Drought"   = "Post Drought"
    ),
    Ryegrass  = factor(Ryegrass, levels = c("Nil", "AR37")),
    Water     = factor(Water, levels = c("Irrigated", "Drought")),
    Treatment = factor(Treatment, levels = c("Control", "Polymer", "6-PP")),
    Harvest   = factor(Harvest, levels = c("During Drought", "Post Drought"))
  )

message("✅ Sample names match between count matrix and metadata.")
message("Treatment counts:")
print(table(metadata$Treatment, useNA = "ifany"))
message("Harvest counts:")
print(table(metadata$Harvest, useNA = "ifany"))

write.csv(metadata, file.path(qc_dir, "Metadata_checked.csv"), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 7️⃣ Load Lolium perenne annotation database
# ─────────────────────────────────────────────────────────────────────────────
ah <- AnnotationHub()
lper_db <- ah[["AH117381"]]

gene_info <- AnnotationDbi::select(
  lper_db,
  keys = rownames(count_mat),
  keytype = "SYMBOL",
  columns = c("GENENAME", "ENTREZID")
) %>%
  dplyr::distinct(SYMBOL, .keep_all = TRUE)

gene_annot <- tibble::tibble(SYMBOL = rownames(count_mat)) %>%
  dplyr::left_join(gene_info, by = "SYMBOL") %>%
  dplyr::mutate(
    Display = dplyr::case_when(
      !is.na(GENENAME) & GENENAME != "" ~ GENENAME,
      TRUE ~ SYMBOL
    )
  )

write.csv(gene_annot, file.path(qc_dir, "Gene_annotation_table.csv"), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 8️⃣ Create DESeq2 object and run model
# ─────────────────────────────────────────────────────────────────────────────
dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = metadata,
  design    = ~ Ryegrass + Harvest + Water * Treatment + Ryegrass:Water
)

message("📌 Genes before filtering: ", nrow(dds))

keep <- rowSums(counts(dds) >= 10) >= 8
dds <- dds[keep, ]

message("🧹 Retained genes after filtering: ", nrow(dds))

dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file.path(qc_dir, "Normalized_Counts_Ryegrass.csv"))

vsd <- vst(dds, blind = FALSE)

writeLines(resultsNames(dds), file.path(rep_dir, "DESeq2_resultsNames.txt"))

# ─────────────────────────────────────────────────────────────────────────────
# 9️⃣ PCA
# ─────────────────────────────────────────────────────────────────────────────
pcaData <- plotPCA(
  vsd,
  intgroup = c("Ryegrass", "Water", "Treatment", "Harvest"),
  ntop = 500,
  returnData = TRUE
)

percentVar <- round(100 * attr(pcaData, "percentVar"))

gg_pca <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Treatment, shape = Water)) +
  geom_point(size = 3, alpha = 0.85) +
  facet_wrap(~ Ryegrass + Harvest, ncol = 2) +
  labs(
    title = "PCA of ryegrass transcriptomes",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance")
  ) +
  theme_pub(13)

save_plot(gg_pca, "PCA_Ryegrass_VST", out_path = pca_dir, w = 8, h = 6)

# ─────────────────────────────────────────────────────────────────────────────
# 9️⃣B PCA highlighting endophyte effects (Nil vs AR37)
# ─────────────────────────────────────────────────────────────────────────────

pca_endo <- ggplot(
  pcaData,
  aes(
    x = PC1,
    y = PC2,
    color = Ryegrass,
    shape = Water
  )
) +
  geom_point(size = 3.5, alpha = 0.9) +
  facet_grid(Harvest ~ Treatment) +
  scale_color_manual(values = c(
    "Nil"  = "#0072B2",
    "AR37" = "#D55E00"
  )) +
  labs(
    title = "PCA of ryegrass transcriptomes by endophyte status",
    x = paste0("PC1: ", percentVar[1], "% variance"),
    y = paste0("PC2: ", percentVar[2], "% variance"),
    color = "Endophyte",
    shape = "Water regime"
  ) +
  theme_pub(13)

save_plot(
  pca_endo,
  "PCA_Endophyte_Nil_vs_AR37",
  out_path = pca_dir,
  w = 10,
  h = 7
)

# ─────────────────────────────────────────────────────────────────────────────
# 9️⃣C PERMANOVA of transcriptome variation
# ─────────────────────────────────────────────────────────────────────────────

transcriptome_dist <- dist(t(assay(vsd)))

# Global model: all experimental factors together
transcriptome_permanova_global <- adonis2(
  transcriptome_dist ~ Ryegrass + Water + Treatment + Harvest,
  data = metadata,
  permutations = 999
)

print(transcriptome_permanova_global)

write.csv(
  as.data.frame(transcriptome_permanova_global),
  file.path(pca_dir, "Transcriptome_PERMANOVA_global_model.csv")
)

# Endophyte-only model: direct Nil vs AR37 comparison
transcriptome_permanova_endophyte <- adonis2(
  transcriptome_dist ~ Ryegrass,
  data = metadata,
  permutations = 999
)

print(transcriptome_permanova_endophyte)

write.csv(
  as.data.frame(transcriptome_permanova_endophyte),
  file.path(pca_dir, "Transcriptome_PERMANOVA_Endophyte_only.csv")
)

# ─────────────────────────────────────────────────────────────────────────────
# 🔟 Sample distance heatmap
# ─────────────────────────────────────────────────────────────────────────────
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

ann_cols <- metadata[, c("Ryegrass", "Water", "Treatment", "Harvest"), drop = FALSE]
rownames(ann_cols) <- metadata$Sample

pdf(file.path(pca_dir, "SampleDistanceHeatmap.pdf"), width = 8, height = 7)
pheatmap(
  sampleDistMatrix,
  annotation_col = ann_cols,
  annotation_row = ann_cols,
  main = "Sample-to-sample distances (VST)"
)
dev.off()

png(file.path(pca_dir, "SampleDistanceHeatmap.png"), width = 2400, height = 2100, res = 300)
pheatmap(
  sampleDistMatrix,
  annotation_col = ann_cols,
  annotation_row = ann_cols,
  main = "Sample-to-sample distances (VST)"
)
dev.off()

# ─────────────────────────────────────────────────────────────────────────────
# 1️⃣1️⃣ HEATMAPS (publication-ready)
# ─────────────────────────────────────────────────────────────────────────────
expr_cols <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "RdBu")))(100)

ann_colors_full <- list(
  Treatment = c("Control" = "#27b7e5", "Polymer" = "#24d0b5", "6-PP" = "#7aa6f0"),
  Harvest   = c("During Drought" = "#a8bf00", "Post Drought" = "#e58df2"),
  Water     = c("Irrigated" = "#19d34a", "Drought" = "#f2b000"),
  Ryegrass  = c("Nil" = "#f3a59a", "AR37" = "#e778d9")
)

# A. Top 100 variable genes across all samples
top_n_var <- 100
topVar_idx <- head(order(matrixStats::rowVars(assay(vsd)), decreasing = TRUE), top_n_var)

mat_all <- assay(vsd)[topVar_idx, , drop = FALSE]
mat_all <- mat_all[complete.cases(mat_all), , drop = FALSE]

all_gene_ids <- rownames(mat_all)
all_gene_labels <- gene_annot$Display[match(all_gene_ids, gene_annot$SYMBOL)]
rownames(mat_all) <- make_unique_labels(all_gene_ids, all_gene_labels)

heat_annot_all <- metadata[, c("Treatment", "Harvest", "Water", "Ryegrass"), drop = FALSE]
rownames(heat_annot_all) <- metadata$Sample

heat_annot_all <- heat_annot_all %>%
  dplyr::mutate(
    Treatment = factor(Treatment, levels = c("Control", "Polymer", "6-PP")),
    Harvest   = factor(Harvest, levels = c("During Drought", "Post Drought")),
    Water     = factor(Water, levels = c("Irrigated", "Drought")),
    Ryegrass  = factor(Ryegrass, levels = c("Nil", "AR37"))
  ) %>%
  dplyr::arrange(Treatment, Harvest, Water, Ryegrass)

mat_all <- mat_all[, rownames(heat_annot_all), drop = FALSE]

treat_counts_all <- table(heat_annot_all$Treatment)
gaps_col_all <- cumsum(as.numeric(treat_counts_all))
gaps_col_all <- gaps_col_all[-length(gaps_col_all)]

pheatmap(
  mat_all,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = heat_annot_all,
  annotation_colors = ann_colors_full,
  gaps_col = gaps_col_all,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 5,
  border_color = NA,
  color = expr_cols,
  main = "Top 100 variable genes",
  filename = file.path(heat_dir, "Heatmap_Top100_GroupedByTreatment.pdf"),
  width = 12,
  height = 12
)

pheatmap(
  mat_all,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = heat_annot_all,
  annotation_colors = ann_colors_full,
  gaps_col = gaps_col_all,
  show_rownames = TRUE,
  show_colnames = FALSE,
  fontsize_row = 5,
  border_color = NA,
  color = expr_cols,
  main = "Top 100 variable genes",
  filename = file.path(heat_dir, "Heatmap_Top100_GroupedByTreatment.png"),
  width = 12,
  height = 12
)

# B. Top 50 variable genes averaged by treatment group
meta_vsd <- as.data.frame(colData(vsd))
meta_vsd$Group <- paste(
  meta_vsd$Treatment,
  meta_vsd$Harvest,
  meta_vsd$Water,
  meta_vsd$Ryegrass,
  sep = "__"
)

mat_vst <- assay(vsd)

grouped_means <- sapply(
  split(seq_len(ncol(mat_vst)), meta_vsd$Group),
  function(idx) rowMeans(mat_vst[, idx, drop = FALSE])
)
grouped_means <- as.matrix(grouped_means)

top_n_avg <- 50
topVar_avg <- head(order(matrixStats::rowVars(grouped_means), decreasing = TRUE), top_n_avg)
mat_avg <- grouped_means[topVar_avg, , drop = FALSE]

avg_gene_ids <- rownames(mat_avg)
avg_gene_labels <- gene_annot$Display[match(avg_gene_ids, gene_annot$SYMBOL)]
rownames(mat_avg) <- make_unique_labels(avg_gene_ids, avg_gene_labels)

ann_group <- data.frame(Group = colnames(mat_avg)) %>%
  tidyr::separate(
    Group,
    into = c("Treatment", "Harvest", "Water", "Ryegrass"),
    sep = "__",
    remove = FALSE
  )

ann_group <- ann_group %>%
  dplyr::mutate(
    Treatment = factor(Treatment, levels = c("Control", "Polymer", "6-PP")),
    Harvest   = factor(Harvest, levels = c("During Drought", "Post Drought")),
    Water     = factor(Water, levels = c("Irrigated", "Drought")),
    Ryegrass  = factor(Ryegrass, levels = c("Nil", "AR37"))
  ) %>%
  dplyr::arrange(Treatment, Harvest, Water, Ryegrass)

rownames(ann_group) <- ann_group$Group
ann_group <- ann_group[, c("Treatment", "Harvest", "Water", "Ryegrass"), drop = FALSE]

mat_avg <- mat_avg[, rownames(ann_group), drop = FALSE]

treat_counts_avg <- table(ann_group$Treatment)
gaps_col_avg <- cumsum(as.numeric(treat_counts_avg))
gaps_col_avg <- gaps_col_avg[-length(gaps_col_avg)]

pheatmap(
  mat_avg,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = ann_group,
  annotation_colors = ann_colors_full,
  gaps_col = gaps_col_avg,
  show_rownames = TRUE,
  fontsize_row = 6,
  fontsize_col = 8,
  border_color = NA,
  color = expr_cols,
  main = "Top 50 variable genes (averaged by treatment)",
  filename = file.path(heat_dir, "Heatmap_Top50_Averaged_GroupedByTreatment.pdf"),
  width = 12,
  height = 8
)

pheatmap(
  mat_avg,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = ann_group,
  annotation_colors = ann_colors_full,
  gaps_col = gaps_col_avg,
  show_rownames = TRUE,
  fontsize_row = 6,
  fontsize_col = 8,
  border_color = NA,
  color = expr_cols,
  main = "Top 50 variable genes (averaged by treatment)",
  filename = file.path(heat_dir, "Heatmap_Top50_Averaged_GroupedByTreatment.png"),
  width = 12,
  height = 8
)

# helper for focused DE heatmaps
make_focused_de_heatmap <- function(res_df,
                                    vsd_obj,
                                    metadata_df,
                                    sample_filter_expr,
                                    contrast_title,
                                    out_prefix,
                                    top_n = 30) {
  
  if (is.null(res_df) || nrow(res_df) == 0) {
    message("⚠️ No DE results supplied for heatmap: ", out_prefix)
    return(NULL)
  }
  
  sig_df <- res_df %>%
    dplyr::filter(!is.na(padj), padj < 0.05) %>%
    dplyr::arrange(padj, dplyr::desc(abs(log2FoldChange))) %>%
    dplyr::slice_head(n = top_n)
  
  if (nrow(sig_df) < 2) {
    message("⚠️ Too few significant genes for heatmap: ", out_prefix)
    return(NULL)
  }
  
  meta_sub <- metadata_df %>%
    dplyr::filter(!!rlang::parse_expr(sample_filter_expr)) %>%
    dplyr::mutate(
      Treatment = factor(Treatment, levels = c("Control", "6-PP")),
      Ryegrass  = factor(Ryegrass, levels = c("Nil", "AR37")),
      Water     = factor(Water, levels = c("Irrigated", "Drought")),
      Harvest   = factor(Harvest, levels = c("During Drought", "Post Drought"))
    ) %>%
    dplyr::arrange(Treatment, Ryegrass)
  
  sample_ids <- meta_sub$Sample
  
  mat_sub <- assay(vsd_obj)[sig_df$SYMBOL, sample_ids, drop = FALSE]
  
  gene_ids <- rownames(mat_sub)
  gene_labels <- gene_annot$Display[match(gene_ids, gene_annot$SYMBOL)]
  rownames(mat_sub) <- make_unique_labels(gene_ids, gene_labels)
  
  ann_sub <- meta_sub[, c("Treatment", "Ryegrass", "Water", "Harvest"), drop = FALSE]
  rownames(ann_sub) <- meta_sub$Sample
  
  ann_colors_sub <- list(
    Treatment = c("Control" = "#27b7e5", "6-PP" = "#7aa6f0"),
    Ryegrass  = c("Nil" = "#f3a59a", "AR37" = "#e778d9"),
    Water     = c("Irrigated" = "#19d34a", "Drought" = "#f2b000"),
    Harvest   = c("During Drought" = "#a8bf00", "Post Drought" = "#e58df2")
  )
  
  trt_counts <- table(ann_sub$Treatment)
  gaps_sub <- if (length(trt_counts) > 1) cumsum(as.numeric(trt_counts))[-length(trt_counts)] else NULL
  
  pheatmap(
    mat_sub,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = ann_sub,
    annotation_colors = ann_colors_sub,
    gaps_col = gaps_sub,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 7,
    border_color = NA,
    color = expr_cols,
    main = contrast_title,
    filename = file.path(heat_dir, paste0(out_prefix, ".pdf")),
    width = 9,
    height = 10
  )
  
  pheatmap(
    mat_sub,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = ann_sub,
    annotation_colors = ann_colors_sub,
    gaps_col = gaps_sub,
    show_rownames = TRUE,
    show_colnames = FALSE,
    fontsize_row = 7,
    border_color = NA,
    color = expr_cols,
    main = contrast_title,
    filename = file.path(heat_dir, paste0(out_prefix, ".png")),
    width = 9,
    height = 10
  )
  
  invisible(mat_sub)
}

# ─────────────────────────────────────────────────────────────────────────────
# 1️⃣2️⃣ DE helper: annotate + volcano
# ─────────────────────────────────────────────────────────────────────────────
extract_and_annotate <- function(res_obj, contrast_label, top_label_n = 20) {
  
  res_df <- as.data.frame(res_obj) %>%
    tibble::rownames_to_column("SYMBOL") %>%
    dplyr::left_join(gene_annot, by = "SYMBOL") %>%
    dplyr::mutate(
      Display = dplyr::case_when(
        !is.na(GENENAME) & GENENAME != "" ~ GENENAME,
        TRUE ~ SYMBOL
      ),
      padj_plot = ifelse(is.na(padj), 1, padj),
      neglog10_padj = -log10(padj_plot),
      significance = "Not Significant"
    )
  
  res_df$significance[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
  res_df$significance[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"
  
  write.csv(res_df, file.path(de_dir, paste0("DESeq2_", contrast_label, "_annotated.csv")), row.names = FALSE)
  
  top_labels <- res_df %>%
    dplyr::filter(!is.na(padj), padj < 0.05) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = top_label_n)
  
  volcano_title <- stringr::str_wrap(
    paste("Volcano plot —", contrast_label),
    width = 42
  )
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = neglog10_padj, color = significance)) +
    geom_point(alpha = 0.65, size = 1.5) +
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(label = Display),
      size = 3,
      max.overlaps = 80,
      box.padding = 0.4,
      point.padding = 0.2,
      segment.color = "grey50"
    ) +
    scale_color_manual(values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not Significant" = "grey70"
    )) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(
      title = volcano_title,
      x = "log2 fold change",
      y = expression(-log[10]("adjusted p-value"))
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 11, lineheight = 1.0),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 12, 10, 12)
    )
  
  save_plot(p, paste0("Volcano_", contrast_label), out_path = de_dir, w = 9, h = 6.5)
  
  sig_count <- sum(!is.na(res_df$padj) & res_df$padj < 0.05)
  message("🧬 ", contrast_label, " | Significant genes (FDR < 0.05): ", sig_count)
  
  invisible(res_df)
}

# ─────────────────────────────────────────────────────────────────────────────
# 1️⃣2️⃣B Global treatment model without Water × Treatment interaction
# ─────────────────────────────────────────────────────────────────────────────
dds_global_trt <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData   = metadata,
  design    = ~ Ryegrass + Harvest + Water + Treatment + Ryegrass:Water
)

keep_global <- rowSums(counts(dds_global_trt) >= 10) >= 8
dds_global_trt <- dds_global_trt[keep_global, ]
dds_global_trt <- DESeq(dds_global_trt)

res_6PP_global <- results(
  dds_global_trt,
  contrast = c("Treatment", "6-PP", "Control")
)

res_6PP_global_annot <- extract_and_annotate(
  res_6PP_global,
  "6-PP_vs_Control_GLOBAL"
)

# ─────────────────────────────────────────────────────────────────────────────
# 1️⃣3️⃣ Main-effect contrasts from interaction model
# ─────────────────────────────────────────────────────────────────────────────
res_drought <- results(dds, contrast = c("Water", "Drought", "Irrigated"))
res_polymer <- results(dds, contrast = c("Treatment", "Polymer", "Control"))
res_6PP_ref <- results(dds, contrast = c("Treatment", "6-PP", "Control"))

res_drought_annot <- extract_and_annotate(res_drought, "Drought_vs_Irrigated")
res_polymer_annot <- extract_and_annotate(res_polymer, "Polymer_vs_Control")
res_6PP_ref_annot <- extract_and_annotate(res_6PP_ref, "6-PP_vs_Control_at_IrrigatedReference")

# ─────────────────────────────────────────────────────────────────────────────
# 1️⃣4️⃣ Stratified DE analyses
# ─────────────────────────────────────────────────────────────────────────────
run_subset_deseq <- function(count_matrix, meta_df, subset_expr, design_formula, contrast_vec, out_prefix) {
  
  meta_sub <- meta_df %>% dplyr::filter(!!rlang::parse_expr(subset_expr))
  counts_sub <- count_matrix[, meta_sub$Sample, drop = FALSE]
  
  stopifnot(all(colnames(counts_sub) == meta_sub$Sample))
  
  dds_sub <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData   = meta_sub,
    design    = design_formula
  )
  
  keep_sub <- rowSums(counts(dds_sub) >= 10) >= 4
  dds_sub <- dds_sub[keep_sub, ]
  dds_sub <- DESeq(dds_sub)
  
  res_sub <- results(dds_sub, contrast = contrast_vec)
  res_sub_df <- extract_and_annotate(res_sub, out_prefix)
  
  invisible(list(dds = dds_sub, res = res_sub, res_df = res_sub_df))
}

subset_6PP_drought <- run_subset_deseq(
  count_matrix = count_mat,
  meta_df = metadata,
  subset_expr = 'Water == "Drought" & Harvest == "During Drought"',
  design_formula = ~ Ryegrass + Treatment,
  contrast_vec = c("Treatment", "6-PP", "Control"),
  out_prefix = "6-PP_vs_Control_within_Drought_DuringDrought"
)

subset_polymer_drought <- run_subset_deseq(
  count_matrix = count_mat,
  meta_df = metadata,
  subset_expr = 'Water == "Drought" & Harvest == "During Drought"',
  design_formula = ~ Ryegrass + Treatment,
  contrast_vec = c("Treatment", "Polymer", "Control"),
  out_prefix = "Polymer_vs_Control_within_Drought_DuringDrought"
)

subset_drought_control <- run_subset_deseq(
  count_matrix = count_mat,
  meta_df = metadata,
  subset_expr = 'Treatment == "Control" & Harvest == "During Drought"',
  design_formula = ~ Ryegrass + Water,
  contrast_vec = c("Water", "Drought", "Irrigated"),
  out_prefix = "Drought_vs_Irrigated_within_Control_DuringDrought"
)

subset_6PP_post <- run_subset_deseq(
  count_matrix = count_mat,
  meta_df = metadata,
  subset_expr = 'Water == "Drought" & Harvest == "Post Drought"',
  design_formula = ~ Ryegrass + Treatment,
  contrast_vec = c("Treatment", "6-PP", "Control"),
  out_prefix = "6-PP_vs_Control_within_Drought_PostDrought"
)

make_focused_de_heatmap(
  res_df = subset_6PP_drought$res_df,
  vsd_obj = vsd,
  metadata_df = metadata,
  sample_filter_expr = 'Water == "Drought" & Harvest == "During Drought" & Treatment %in% c("Control", "6-PP")',
  contrast_title = "Top 30 DE genes: 6-PP vs Control within Drought (During Drought)",
  out_prefix = "Heatmap_Top30_6-PP_vs_Control_within_Drought_DuringDrought",
  top_n = 30
)

make_focused_de_heatmap(
  res_df = subset_6PP_post$res_df,
  vsd_obj = vsd,
  metadata_df = metadata,
  sample_filter_expr = 'Water == "Drought" & Harvest == "Post Drought" & Treatment %in% c("Control", "6-PP")',
  contrast_title = "Top 30 DE genes: 6-PP vs Control within Drought (Post Drought)",
  out_prefix = "Heatmap_Top30_6-PP_vs_Control_within_Drought_PostDrought",
  top_n = 30
)

# ─────────────────────────────────────────────────────────────────────────────
# 1️⃣4️⃣B Reproducible combined supplementary volcano (Figure S14)
# ─────────────────────────────────────────────────────────────────────────────
p_s14_during <- make_simple_volcano(
  subset_6PP_drought$res_df,
  "During Drought"
)

p_s14_post <- make_simple_volcano(
  subset_6PP_post$res_df,
  "Post Drought"
)

fig_s14 <- (p_s14_during | p_s14_post) +
  plot_annotation(
    title = "Differential gene expression between 6-PP and control treatments",
  ) &
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )

ggsave(
  filename = file.path(fig_dir, "Figure_S14_6-PP_vs_Control_During_and_PostDrought.pdf"),
  plot = fig_s14,
  width = 13,
  height = 6,
  useDingbats = FALSE
)

ggsave(
  filename = file.path(fig_dir, "Figure_S14_6-PP_vs_Control_During_and_PostDrought.png"),
  plot = fig_s14,
  width = 13,
  height = 6,
  dpi = 300
)

# ─────────────────────────────────────────────────────────────────────────────
# 1️⃣5️⃣ MA plots
# ─────────────────────────────────────────────────────────────────────────────
pdf(file.path(de_dir, "MAplots_mainContrasts.pdf"), width = 8, height = 6)
plotMA(res_drought, main = "Drought vs Irrigated", ylim = c(-5, 5))
plotMA(res_polymer, main = "Polymer vs Control", ylim = c(-5, 5))
plotMA(res_6PP_ref, main = "6-PP vs Control (Irrigated reference)", ylim = c(-5, 5))
plotMA(res_6PP_global, main = "6-PP vs Control (Global)", ylim = c(-5, 5))
dev.off()

# ─────────────────────────────────────────────────────────────────────────────
# 1️⃣6️⃣ Focused heatmap of significant genes from global 6-PP vs Control
# ─────────────────────────────────────────────────────────────────────────────
sig_6PP_global <- res_6PP_global_annot %>%
  dplyr::filter(!is.na(padj), padj < 0.05) %>%
  dplyr::arrange(padj) %>%
  dplyr::slice_head(n = 50)

if (nrow(sig_6PP_global) > 1) {
  mat_sig <- assay(vsd)[sig_6PP_global$SYMBOL, , drop = FALSE]
  rownames(mat_sig) <- make_unique_labels(sig_6PP_global$SYMBOL, sig_6PP_global$Display)
  
  ann_sig <- metadata[, c("Ryegrass", "Water", "Treatment", "Harvest"), drop = FALSE]
  rownames(ann_sig) <- metadata$Sample
  
  pheatmap(
    mat_sig,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = ann_sig,
    show_rownames = TRUE,
    fontsize_row = 6,
    main = "Top significant genes: 6-PP vs Control (Global)",
    color = expr_cols,
    filename = file.path(heat_dir, "Heatmap_TopSig_6-PP_vs_Control_GLOBAL.pdf"),
    width = 10,
    height = 10
  )
  
  pheatmap(
    mat_sig,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    annotation_col = ann_sig,
    show_rownames = TRUE,
    fontsize_row = 6,
    main = "Top significant genes: 6-PP vs Control (Global)",
    color = expr_cols,
    filename = file.path(heat_dir, "Heatmap_TopSig_6-PP_vs_Control_GLOBAL.png"),
    width = 10,
    height = 10
  )
}

# ─────────────────────────────────────────────────────────────────────────────
# 1️⃣7️⃣ GO enrichment helper
# ─────────────────────────────────────────────────────────────────────────────
run_go_enrichment <- function(res_df, contrast_label, ont = "BP", q_cut = 0.05, top_n = 20) {
  
  sig_genes <- res_df %>%
    dplyr::filter(!is.na(padj), padj < q_cut, !is.na(ENTREZID)) %>%
    dplyr::pull(ENTREZID) %>%
    unique()
  
  universe_genes <- res_df %>%
    dplyr::filter(!is.na(ENTREZID)) %>%
    dplyr::pull(ENTREZID) %>%
    unique()
  
  if (length(sig_genes) < 5) {
    message("⚠️ Too few significant annotated genes for GO enrichment: ", contrast_label)
    return(NULL)
  }
  
  ego <- enrichGO(
    gene          = sig_genes,
    universe      = universe_genes,
    OrgDb         = lper_db,
    keyType       = "ENTREZID",
    ont           = ont,
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.20,
    readable      = TRUE
  )
  
  ego_df <- as.data.frame(ego)
  
  if (nrow(ego_df) == 0) {
    message("⚠️ No GO terms enriched for: ", contrast_label)
    return(NULL)
  }
  
  write.csv(
    ego_df,
    file.path(go_dir, paste0("GO_", ont, "_", contrast_label, ".csv")),
    row.names = FALSE
  )
  
  try({
    p_dot <- suppressWarnings(
      dotplot(ego, showCategory = top_n) +
        ggtitle(paste0("GO-", ont, " enrichment — ", contrast_label)) +
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
    )
    save_plot(p_dot, paste0("GO_", ont, "_dotplot_", contrast_label), out_path = go_dir, w = 9, h = 7)
  }, silent = TRUE)
  
  try({
    p_bar <- suppressWarnings(
      barplot(ego, showCategory = top_n) +
        ggtitle(paste0("GO-", ont, " enrichment — ", contrast_label)) +
        theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
    )
    save_plot(p_bar, paste0("GO_", ont, "_barplot_", contrast_label), out_path = go_dir, w = 10, h = 7)
  }, silent = TRUE)
  
  invisible(ego_df)
}

# ─────────────────────────────────────────────────────────────────────────────
# 1️⃣8️⃣ Run GO enrichment for key contrasts
# ─────────────────────────────────────────────────────────────────────────────
go_drought <- run_go_enrichment(res_drought_annot, "Drought_vs_Irrigated", ont = "BP")
go_polymer <- run_go_enrichment(res_polymer_annot, "Polymer_vs_Control", ont = "BP")
go_6PP_global <- run_go_enrichment(res_6PP_global_annot, "6-PP_vs_Control_GLOBAL", ont = "BP")

go_6PP_drought <- run_go_enrichment(
  subset_6PP_drought$res_df,
  "6-PP_vs_Control_within_Drought_DuringDrought",
  ont = "BP"
)

go_6PP_post <- run_go_enrichment(
  subset_6PP_post$res_df,
  "6-PP_vs_Control_within_Drought_PostDrought",
  ont = "BP"
)

# ─────────────────────────────────────────────────────────────────────────────
# 1️⃣9️⃣ Summary table of DE gene counts
# ─────────────────────────────────────────────────────────────────────────────
de_summary <- tibble::tibble(
  Contrast = c(
    "Drought_vs_Irrigated",
    "Polymer_vs_Control",
    "6-PP_vs_Control_at_IrrigatedReference",
    "6-PP_vs_Control_GLOBAL",
    "6-PP_vs_Control_within_Drought_DuringDrought",
    "6-PP_vs_Control_within_Drought_PostDrought",
    "Polymer_vs_Control_within_Drought_DuringDrought",
    "Drought_vs_Irrigated_within_Control_DuringDrought"
  ),
  Significant_FDR_0_05 = c(
    sum(!is.na(res_drought_annot$padj) & res_drought_annot$padj < 0.05),
    sum(!is.na(res_polymer_annot$padj) & res_polymer_annot$padj < 0.05),
    sum(!is.na(res_6PP_ref_annot$padj) & res_6PP_ref_annot$padj < 0.05),
    sum(!is.na(res_6PP_global_annot$padj) & res_6PP_global_annot$padj < 0.05),
    sum(!is.na(subset_6PP_drought$res_df$padj) & subset_6PP_drought$res_df$padj < 0.05),
    sum(!is.na(subset_6PP_post$res_df$padj) & subset_6PP_post$res_df$padj < 0.05),
    sum(!is.na(subset_polymer_drought$res_df$padj) & subset_polymer_drought$res_df$padj < 0.05),
    sum(!is.na(subset_drought_control$res_df$padj) & subset_drought_control$res_df$padj < 0.05)
  )
)

write.csv(de_summary, file.path(de_dir, "DESeq2_summary_counts.csv"), row.names = FALSE)

# ─────────────────────────────────────────────────────────────────────────────
# 2️⃣0️⃣ Save session info
# ─────────────────────────────────────────────────────────────────────────────
writeLines(capture.output(sessionInfo()), file.path(rep_dir, "R_sessionInfo_Ryegrass.txt"))

message("✅ Publication-ready RNA-seq analysis complete.")

# ─────────────────────────────────────────────────────────────────────────────
# 2️⃣1️⃣ Figure 2 — Nature Communications host transcriptome figure
# A = PCA during drought
# B = PCA post drought
# C = Volcano: 6-PP vs Control during drought
# D = Volcano: 6-PP vs Control post drought
# E = GO enrichment during drought
# F = GO enrichment post drought
# ─────────────────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  
  library(patchwork)
  library(cowplot)
  library(ggrepel)
  library(enrichplot)
  
})

fig2_dir <- file.path(project_dir, "7.manuscript_figure", "3.figure2")
dir.create(fig2_dir, recursive = TRUE, showWarnings = FALSE)

theme_natcom <- function(base_size = 15) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold", colour = "black", size = base_size),
      axis.text = element_text(colour = "black", size = base_size - 2),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold", colour = "black", size = base_size - 1),
      legend.title = element_text(face = "bold", size = base_size - 1),
      legend.text = element_text(size = base_size - 2),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 3),
      plot.subtitle = element_text(face = "bold", hjust = 0, size = base_size - 1),
      plot.margin = margin(8, 8, 8, 8)
    )
}

# PCA helper
make_pca_panel <- function(harvest_level, panel_letter) {
  
  ggplot(
    pcaData %>% dplyr::filter(Harvest == harvest_level),
    aes(x = PC1, y = PC2, colour = Treatment, shape = Water)
  ) +
    geom_point(size = 4, alpha = 0.9) +
    facet_wrap(~ Ryegrass, nrow = 1) +
    scale_colour_manual(values = c(
      "Control" = "#F8766D",
      "Polymer" = "#00BA38",
      "6-PP" = "#619CFF"
    )) +
    labs(
      title = panel_letter,
      subtitle = harvest_level,
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      colour = "Seed treatment",
      shape = "Water"
    ) +
    theme_natcom(15) +
    theme(
      plot.subtitle = element_text(face = "bold", size = 14, hjust = 0),
      legend.position = "right"
    )
}

panel_a <- make_pca_panel("During Drought", "A")
panel_b <- make_pca_panel("Post Drought", "B")

# Volcano helper
make_natcom_volcano <- function(res_df, panel_letter, subtitle_text, label_n = 8) {
  
  plot_df <- res_df %>%
    dplyr::mutate(
      padj_plot = ifelse(is.na(padj), 1, padj),
      neglog10_padj = -log10(padj_plot),
      significance = dplyr::case_when(
        !is.na(padj) & padj < 0.05 & log2FoldChange > 1 ~ "Upregulated",
        !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  top_labels <- plot_df %>%
    dplyr::filter(!is.na(padj), padj < 0.05) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = label_n)
  
  ggplot(plot_df, aes(x = log2FoldChange, y = neglog10_padj, colour = significance)) +
    geom_point(size = 1.8, alpha = 0.7) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", linewidth = 0.45) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", linewidth = 0.45) +
    ggrepel::geom_text_repel(
      data = top_labels,
      aes(label = Display),
      size = 3.4,
      max.overlaps = 50,
      box.padding = 0.35,
      point.padding = 0.25,
      segment.color = "grey50"
    ) +
    scale_colour_manual(values = c(
      "Upregulated" = "red",
      "Downregulated" = "blue",
      "Not significant" = "grey70"
    )) +
    labs(
      title = panel_letter,
      subtitle = subtitle_text,
      x = expression(log[2]~fold~change),
      y = expression(-log[10]~adjusted~italic(P)),
      colour = "Significance"
    ) +
    theme_natcom(15) +
    theme(
      legend.position = "right"
    )
}

panel_c <- make_natcom_volcano(
  subset_6PP_drought$res_df,
  panel_letter = "C",
  subtitle_text = "6-PP vs Control during drought"
)

panel_d <- make_natcom_volcano(
  subset_6PP_post$res_df,
  panel_letter = "D",
  subtitle_text = "6-PP vs Control post drought"
)

# GO plotting helper
make_go_dotplot <- function(go_df, panel_letter, subtitle_text, show_n = 10) {
  
  if (is.null(go_df) || nrow(go_df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = "No enriched GO terms", size = 6) +
        labs(title = panel_letter, subtitle = subtitle_text) +
        theme_void() +
        theme(
          plot.title = element_text(face = "bold", size = 18, hjust = 0),
          plot.subtitle = element_text(face = "bold", size = 14, hjust = 0)
        )
    )
  }
  
  plot_df <- go_df %>%
    dplyr::arrange(p.adjust) %>%
    dplyr::slice_head(n = show_n) %>%
    dplyr::mutate(
      Description = stringr::str_wrap(Description, width = 20),
      Description = factor(Description, levels = rev(Description)),
      GeneRatio_num = sapply(
        strsplit(GeneRatio, "/"),
        function(x) as.numeric(x[1]) / as.numeric(x[2])
      )
    )
  
  ggplot(
    plot_df,
    aes(
      x = GeneRatio_num,
      y = Description,
      size = Count,
      colour = p.adjust
    )
  ) +
    geom_point(alpha = 0.9, stroke = 0.3) +
    scale_colour_gradient(
      low = "red",
      high = "blue",
      name = "Adjusted P"
    ) +
    scale_size_continuous(
      name = "Gene count",
      range = c(4, 9)
    ) +
    labs(
      title = panel_letter,
      subtitle = subtitle_text,
      x = "Gene ratio",
      y = NULL
    ) +
    theme_natcom(16) +
    theme(
      axis.text.y = element_text(
        size = 13,
        colour = "black",
        hjust = 1
      ),
      axis.text.x = element_text(
        size = 12,
        colour = "black"
      ),
      axis.title.x = element_text(
        size = 14,
        face = "bold"
      ),
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      legend.text = element_text(
        size = 11
      ),
      plot.subtitle = element_text(
        size = 14,
        face = "bold"
      ),
      legend.position = "right",
      plot.margin = margin(6, 8, 6, 4)
    )
}

panel_e <- make_go_dotplot(
  go_6PP_drought,
  panel_letter = "E",
  subtitle_text = "GO enrichment during drought",
  show_n = 10
)

panel_f <- make_go_dotplot(
  go_6PP_post,
  panel_letter = "F",
  subtitle_text = "GO enrichment post drought",
  show_n = 10
)

#################################################
# Assemble Figure 2
#################################################

top <- plot_grid(
  
  panel_a,
  panel_b,
  
  nrow = 1,
  
  rel_widths = c(1,1)
  
)


middle <- plot_grid(
  
  panel_c,
  panel_d,
  
  nrow = 1,
  
  rel_widths = c(1,1)
  
)


bottom <- plot_grid(
  
  panel_e,
  panel_f,
  
  nrow = 1,
  
  rel_widths = c(0.75,0.75)
  
)


fig2 <- plot_grid(
  top,
  middle,
  bottom,
  ncol = 1,
  rel_heights = c(0.95, 0.95, 1.15)
)

ggsave(
  filename = file.path(fig2_dir, "Figure2_NatCom.pdf"),
  plot = fig2,
  width = 20,
  height = 16,
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  filename = file.path(fig2_dir, "Figure2_NatCom.png"),
  plot = fig2,
  width = 20,
  height = 16,
  dpi = 600,
  limitsize = FALSE
)

write.csv(
  de_summary,
  file.path(fig2_dir, "Figure2_DESeq2_summary_counts.csv"),
  row.names = FALSE
)

message("Done. Figure 2 saved to: ", fig2_dir)

