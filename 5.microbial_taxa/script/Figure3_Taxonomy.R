## =========================================================
## Metatranscriptomics genus-level analysis pipeline
## Publication-ready community structure + differential abundance
##
## INPUTS
## 1) Genus count table (wide): rows = genera, cols = samples
## 2) Metadata table with at least:
##    SampleID, Endophyte, Treatment, Water, Harvest
##    Optional: Flowcell, Batch, Replicate
##
## OUTPUTS
## - QC summaries
## - Alpha diversity tables + figures
## - Relative abundance tables
## - Bray-Curtis and Aitchison ordinations
## - PERMANOVA + pairwise PERMANOVA + dispersion tests
## - ALDEx2 differential abundance outputs
## - Volcano plots
## - Stacked bar plots
## - Top-genus heatmap
## - Reproducibility log
##
## NOTES
## - Suitable for genus-level metatranscriptomics compositional analysis
## - Interpret taxa as transcript-associated representation, not absolute organism abundance
## =========================================================

## -----------------------------
## 0) Setup
## -----------------------------
set.seed(1)

required_pkgs <- c(
  "dplyr","tidyr","readr","stringr","tibble","forcats",
  "ggplot2","ggrepel","patchwork","scales",
  "vegan","ape","ALDEx2","pheatmap","RColorBrewer","ggsci","here"
)

to_install <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

invisible(lapply(required_pkgs, library, character.only = TRUE))

## -----------------------------
## 1) File paths
## -----------------------------
project_dir <- "C:/Users/soths/OneDrive - Lincoln University/Writing/NatMicro/github"

data_dir <- file.path(project_dir, "5.microbial_taxa", "data")
out_dir  <- file.path(project_dir, "5.microbial_taxa", "results")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

qc_dir      <- file.path(out_dir, "01_QC")
alpha_dir   <- file.path(out_dir, "02_AlphaDiversity")
ord_dir     <- file.path(out_dir, "04_Ordination")
perm_dir    <- file.path(out_dir, "05_PERMANOVA")
aldex_dir   <- file.path(out_dir, "06_ALDEx2")
heat_dir    <- file.path(out_dir, "07_Heatmaps")
repr_dir    <- file.path(out_dir, "08_Reproducibility")
fig_dir     <- file.path(out_dir, "09_Publication_Figures")

dir.create(qc_dir,    recursive = TRUE, showWarnings = FALSE)
dir.create(alpha_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(ord_dir,   recursive = TRUE, showWarnings = FALSE)
dir.create(perm_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(aldex_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(heat_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(repr_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,   recursive = TRUE, showWarnings = FALSE)

genus_fp  <- file.path(data_dir, "Genus.tsv")
phylum_fp <- file.path(data_dir, "Kraken2_Phylum_long_with_metadata.csv")
meta_fp   <- file.path(data_dir, "Metadata_tidy.csv")

## -----------------------------
## 2) Analysis parameters
## -----------------------------
min_total_counts_per_sample <- 1000
min_prevalence              <- 0.10
min_total_counts_per_taxon  <- 20
pseudocount                 <- 0.5
aldex_mc_samples            <- 128
fdr_cutoff                  <- 0.05
top_n_barplot               <- 15
top_n_heatmap               <- 40
nperm                       <- 999

## -----------------------------
## 3) Plot theme
## -----------------------------
theme_pub <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(size = base_size - 1, hjust = 0),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", color = "black"),
      legend.title = element_text(face = "bold"),
      legend.key = element_blank(),
      panel.spacing = unit(5, "mm")
    )
}

save_pub <- function(plot_obj, filename, width, height, units = "mm", dpi = 600) {
  ggsave(
    filename, plot_obj,
    width = width, height = height, units = units,
    dpi = dpi, limitsize = FALSE
  )
  if (grepl("\\.png$", filename)) {
    ggsave(
      sub("\\.png$", ".pdf", filename), plot_obj,
      width = width, height = height, units = units,
      dpi = dpi, useDingbats = FALSE, limitsize = FALSE
    )
  }
}

## -----------------------------
## 4) Read metadata
## -----------------------------
meta <- read_csv(meta_fp, show_col_types = FALSE) %>%
  dplyr::mutate(
    SampleID  = stringr::str_trim(SampleID),
    Treatment = dplyr::recode(as.character(Treatment), "VOC" = "6-PP"),
    Endophyte = factor(Endophyte, levels = c("Nil","AR37")),
    Treatment = factor(Treatment, levels = c("Control","Polymer","6-PP")),
    Water     = factor(Water, levels = c("Irrigated","Drought")),
    Harvest   = factor(Harvest, levels = c("During Drought","Post Drought"))
  )

stopifnot(all(c("SampleID","Endophyte","Treatment","Water","Harvest") %in% names(meta)))
stopifnot(!anyDuplicated(meta$SampleID))

has_flowcell <- "Flowcell" %in% names(meta)
has_batch    <- "Batch" %in% names(meta)

## -----------------------------
## 5) Read genus table
## -----------------------------
wide <- read_tsv(genus_fp, show_col_types = FALSE)
names(wide)[1] <- "Genus"

genus_mat <- wide %>%
  dplyr::mutate(Genus = as.character(Genus)) %>%
  column_to_rownames("Genus") %>%
  as.matrix()

colnames(genus_mat) <- str_trim(colnames(genus_mat))
colnames(genus_mat) <- str_replace(colnames(genus_mat), "^([A-H])(\\d)$", "\\10\\2")
mode(genus_mat) <- "numeric"

common_samples <- intersect(colnames(genus_mat), meta$SampleID)
if (length(common_samples) < 2) stop("No matching samples between genus table and metadata.")

genus_mat <- genus_mat[, common_samples, drop = FALSE]
meta_use <- meta %>%
  dplyr::filter(SampleID %in% common_samples) %>%
  dplyr::arrange(match(SampleID, common_samples))

stopifnot(identical(meta_use$SampleID, colnames(genus_mat)))

## -----------------------------
## 6) Initial QC
## -----------------------------
sample_depth <- colSums(genus_mat, na.rm = TRUE)
tax_total    <- rowSums(genus_mat, na.rm = TRUE)
tax_prev     <- rowMeans(genus_mat > 0, na.rm = TRUE)

qc_samples <- tibble(
  SampleID = names(sample_depth),
  Depth = as.numeric(sample_depth)
) %>%
  left_join(meta_use, by = "SampleID") %>%
  dplyr::mutate(Keep = Depth >= min_total_counts_per_sample)

write_csv(qc_samples, file.path(qc_dir, "QC_sample_depths.csv"))

qc_taxa <- tibble(
  Genus = rownames(genus_mat),
  TotalCounts = as.numeric(tax_total),
  Prevalence = as.numeric(tax_prev),
  Keep = (tax_total >= min_total_counts_per_taxon) & (tax_prev >= min_prevalence)
)

write_csv(qc_taxa, file.path(qc_dir, "QC_taxa_filtering.csv"))

p_depth <- ggplot(qc_samples, aes(x = reorder(SampleID, Depth), y = Depth, fill = Keep)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = min_total_counts_per_sample, linetype = "dashed") +
  coord_flip() +
  labs(
    title = "Sequencing depth per sample",
    x = "Sample",
    y = "Total genus-level counts"
  ) +
  scale_y_continuous(labels = comma) +
  scale_fill_manual(values = c("TRUE" = "black", "FALSE" = "grey70")) +
  theme_pub(9) +
  theme(legend.position = "none")

save_pub(p_depth, file.path(qc_dir, "QC_sample_depth.png"), 140, 180)

keep_samples <- sample_depth >= min_total_counts_per_sample
genus_mat <- genus_mat[, keep_samples, drop = FALSE]
meta_use  <- meta_use[keep_samples, , drop = FALSE]

tax_total_f <- rowSums(genus_mat, na.rm = TRUE)
tax_prev_f  <- rowMeans(genus_mat > 0, na.rm = TRUE)
keep_taxa <- (tax_total_f >= min_total_counts_per_taxon) & (tax_prev_f >= min_prevalence)

genus_mat_filt <- genus_mat[keep_taxa, , drop = FALSE]

write_tsv(
  as_tibble(genus_mat_filt, rownames = "Genus"),
  file.path(qc_dir, "Genus_counts_filtered.tsv")
)

## -----------------------------
## 7) Relative abundance + CLR
## -----------------------------
rel_abund <- sweep(genus_mat_filt, 2, colSums(genus_mat_filt), "/")
rel_abund[is.na(rel_abund)] <- 0

clr_mat <- apply(genus_mat_filt + pseudocount, 2, function(x) {
  log(x) - mean(log(x))
})

saveRDS(rel_abund, file.path(repr_dir, "rel_abund_matrix.rds"))
saveRDS(clr_mat, file.path(repr_dir, "clr_matrix.rds"))

## -----------------------------
## 8) Alpha diversity
## -----------------------------
alpha_df <- tibble(
  SampleID = colnames(genus_mat_filt),
  Observed = colSums(genus_mat_filt > 0),
  Shannon  = vegan::diversity(t(genus_mat_filt), index = "shannon"),
  Simpson  = vegan::diversity(t(genus_mat_filt), index = "simpson"),
  Evenness = Shannon / log(pmax(Observed, 2))
) %>%
  left_join(meta_use, by = "SampleID")

write_csv(alpha_df, file.path(alpha_dir, "Alpha_diversity_metrics.csv"))

plot_alpha <- function(df, metric, ylab_text, fname) {
  g <- ggplot(df, aes(x = Water, y = .data[[metric]], fill = Water)) +
    geom_boxplot(
      outlier.shape = NA,
      width = 0.65,
      alpha = 0.85,
      color = "black",
      linewidth = 0.5
    ) +
    geom_jitter(
      aes(shape = Endophyte),
      width = 0.12,
      size = 2.2,
      alpha = 0.9,
      color = "black"
    ) +
    facet_grid(Harvest ~ Treatment) +
    labs(
      title = paste0(metric, " diversity"),
      x = NULL,
      y = ylab_text
    ) +
    scale_fill_manual(values = c("Irrigated" = "#F8766D", "Drought" = "#00BFC4")) +
    scale_shape_manual(values = c("Nil" = 16, "AR37" = 17)) +
    guides(
      fill  = guide_legend(order = 1),
      shape = guide_legend(order = 2)
    ) +
    theme_pub(9)
  
  save_pub(g, file.path(alpha_dir, fname), 240, 140)
  g
}

p_alpha_obs  <- plot_alpha(alpha_df, "Observed", "Observed genera", "Alpha_Observed.png")
p_alpha_shan <- plot_alpha(alpha_df, "Shannon", "Shannon diversity", "Alpha_Shannon.png")
p_alpha_simp <- plot_alpha(alpha_df, "Simpson", "Simpson diversity", "Alpha_Simpson.png")
p_alpha_even <- plot_alpha(alpha_df, "Evenness", "Pielou evenness", "Alpha_Evenness.png")

alpha_lm_results <- list()
for (metric in c("Observed", "Shannon", "Simpson", "Evenness")) {
  fit <- lm(as.formula(paste(metric, "~ Endophyte + Treatment + Water + Harvest")), data = alpha_df)
  alpha_lm_results[[metric]] <- as.data.frame(anova(fit)) %>%
    rownames_to_column("Term")
  
  write_csv(
    alpha_lm_results[[metric]],
    file.path(alpha_dir, paste0("ANOVA_", metric, ".csv"))
  )
}

## -----------------------------
## 9) Publication-style composition figures
## -----------------------------
theme_natureish_comp <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      axis.title.x = element_text(size = base_size, margin = margin(t = 5)),
      axis.title.y = element_text(size = base_size, margin = margin(r = 5)),
      axis.text = element_text(size = base_size, color = "black"),
      axis.line = element_line(linewidth = 0.25),
      axis.ticks = element_line(linewidth = 0.25),
      axis.ticks.length = unit(2, "pt"),
      strip.background = element_blank(),
      strip.text = element_text(size = base_size, face = "bold"),
      strip.clip = "off",
      strip.text.x = element_text(margin = margin(t = 1, r = 2, b = 1, l = 2)),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text  = element_text(size = base_size - 0.5),
      legend.key.size = unit(3, "mm"),
      panel.spacing = unit(4, "mm"),
      plot.margin = margin(6, 6, 6, 6)
    )
}

## =========================================================
## 9A) Genus composition figure
## =========================================================
genus_coverage_target <- 0.90
genus_max_keep        <- 50
genus_min_keep        <- 20
genus_legend_ncol_max <- 4
genus_suffix          <- "COMPREHENSIVE_TOP50"

genus_long <- as_tibble(rel_abund, rownames = "Genus") %>%
  pivot_longer(cols = -Genus, names_to = "SampleID", values_to = "rel_abund") %>%
  left_join(meta_use, by = "SampleID") %>%
  dplyr::mutate(pct_rel = rel_abund * 100)

gen_mean0 <- genus_long %>%
  group_by(Endophyte, Treatment, Water, Harvest, Genus) %>%
  dplyr::summarise(mean_pct = mean(pct_rel, na.rm = TRUE), .groups = "drop")

gen_ranked <- gen_mean0 %>%
  group_by(Genus) %>%
  dplyr::summarise(overall = mean(mean_pct, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(overall)) %>%
  dplyr::mutate(cumprop = cumsum(overall) / sum(overall))

keep_genera <- gen_ranked %>%
  dplyr::filter(cumprop <= genus_coverage_target) %>%
  pull(Genus)

if (length(keep_genera) < genus_min_keep) {
  keep_genera <- gen_ranked %>%
    slice_head(n = genus_min_keep) %>%
    pull(Genus)
}

keep_genera <- head(keep_genera, genus_max_keep)

gen_mean <- gen_mean0 %>%
  dplyr::mutate(Genus2 = ifelse(Genus %in% keep_genera, Genus, "Other")) %>%
  group_by(Endophyte, Treatment, Water, Harvest, Genus2) %>%
  dplyr::summarise(mean_pct = sum(mean_pct), .groups = "drop")

gen_overall_tbl <- gen_mean %>%
  group_by(Genus2) %>%
  dplyr::summarise(overall = mean(mean_pct, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(overall))

gen_levels_big_to_small <- gen_overall_tbl$Genus2
gen_levels_small_to_big <- rev(gen_levels_big_to_small)

gen_mean <- gen_mean %>%
  dplyr::mutate(Genus2 = factor(Genus2, levels = gen_levels_small_to_big))

if ("Other" %in% levels(gen_mean$Genus2)) {
  gen_mean$Genus2 <- forcats::fct_relevel(gen_mean$Genus2, "Other", after = 0)
}

gen_mean <- gen_mean %>%
  dplyr::mutate(Endophyte_plot = as.character(Endophyte))

gen_spacer_rows <- gen_mean %>%
  distinct(Treatment, Water, Harvest, Genus2) %>%
  dplyr::mutate(Endophyte_plot = " ", mean_pct = 0)

gen_mean <- bind_rows(gen_mean, gen_spacer_rows) %>%
  dplyr::mutate(Endophyte_plot = factor(Endophyte_plot, levels = c("Nil", " ", "AR37")))

gen_levels_fill <- levels(gen_mean$Genus2)
gen_legend_order <- gen_overall_tbl %>%
  dplyr::filter(Genus2 != "Other") %>%
  pull(Genus2)

gen_pal <- setNames(rep(NA_character_, length(gen_levels_fill)), gen_levels_fill)

gen_npg10 <- ggsci::pal_npg("nrc")(10)
gen_npg10_rot <- gen_npg10[c(3:10, 1:2)]

top10_genera <- head(gen_legend_order, 10)
gen_pal[top10_genera] <- gen_npg10_rot[seq_len(length(top10_genera))]

remaining_genera <- setdiff(gen_legend_order, top10_genera)

if (length(remaining_genera) > 0) {
  okabe_ito <- c(
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#000000"
  )
  
  extra_pool <- c(
    okabe_ito,
    RColorBrewer::brewer.pal(8, "Dark2"),
    RColorBrewer::brewer.pal(12, "Set3"),
    RColorBrewer::brewer.pal(8, "Accent")
  )
  
  extra_pool <- setdiff(extra_pool, gen_npg10_rot)
  
  too_pale <- c("#FFFFB3", "#FFED6F", "#CCEBC5", "#E5D8BD", "#F2F2F2")
  extra_pool <- setdiff(extra_pool, too_pale)
  
  if (length(extra_pool) < length(remaining_genera)) {
    extra_pool <- c(
      extra_pool,
      grDevices::hcl.colors(
        length(remaining_genera) - length(extra_pool),
        palette = "Dark 3"
      )
    )
  }
  
  gen_pal[remaining_genera] <- extra_pool[seq_len(length(remaining_genera))]
}

gen_pal["Other"] <- "grey75"
gen_pal <- gen_pal[gen_levels_fill]

plot_one_treatment_genus <- function(df, treatment_label, legend_ncol = 2, show_y = TRUE) {
  df1 <- df %>% dplyr::filter(Treatment == treatment_label)
  
  ggplot(df1, aes(x = Endophyte_plot, y = mean_pct, fill = Genus2)) +
    geom_col(width = 1.5, color = NA, position = position_stack(reverse = TRUE)) +
    facet_grid(Harvest ~ Water) +
    coord_cartesian(clip = "off") +
    scale_y_continuous(
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_x_discrete(
      labels = c(Nil = "Nil", " " = "", AR37 = "AR37"),
      expand = expansion(add = 0.30)
    ) +
    scale_fill_manual(values = gen_pal, drop = FALSE) +
    guides(fill = guide_legend(ncol = legend_ncol, reverse = TRUE)) +
    labs(
      title = treatment_label,
      x = NULL,
      y = if (show_y) "Mean relative abundance (%)" else NULL,
      fill = NULL
    ) +
    theme_natureish_comp(base_size = 8) +
    theme(
      axis.text.x = element_text(size = 7),
      plot.title = element_text(face = "bold", size = 9, hjust = 0.5),
      legend.position = "none"
    )
}

gen_legend_cols <- if (length(levels(gen_mean$Genus2)) <= 18) 2 else genus_legend_ncol_max

p_gen_control <- plot_one_treatment_genus(gen_mean, "Control", gen_legend_cols, show_y = TRUE)
p_gen_polymer <- plot_one_treatment_genus(gen_mean, "Polymer", gen_legend_cols, show_y = FALSE)
p_gen_6pp     <- plot_one_treatment_genus(gen_mean, "6-PP", gen_legend_cols, show_y = FALSE)

final_genus_fig <- (p_gen_control | p_gen_polymer | p_gen_6pp) +
  patchwork::plot_annotation(tag_levels = "A") +
  patchwork::plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    plot.margin = margin(12, 16, 10, 10)
  )

ggsave(
  filename = file.path(
    fig_dir,
    paste0("Fig_Genus_ABC_", genus_suffix, "_Coverage", genus_coverage_target, "_Cap", genus_max_keep, ".pdf")
  ),
  plot = final_genus_fig,
  width = 260, height = 110, units = "mm",
  dpi = 600, useDingbats = FALSE, limitsize = FALSE
)

ggsave(
  filename = file.path(
    fig_dir,
    paste0("Fig_Genus_ABC_", genus_suffix, "_Coverage", genus_coverage_target, "_Cap", genus_max_keep, ".png")
  ),
  plot = final_genus_fig,
  width = 260, height = 110, units = "mm",
  dpi = 600, limitsize = FALSE
)

write_csv(
  gen_ranked %>% dplyr::mutate(kept = Genus %in% keep_genera),
  file.path(fig_dir, paste0("Fig_Genus_ABC_", genus_suffix, "_kept_genera.csv"))
)

## =========================================================
## 9B) Phylum composition figure
## =========================================================
if (file.exists(phylum_fp)) {
  
  phy <- read_csv(phylum_fp, show_col_types = FALSE)
  
  req_phy_base <- c("SampleID", "Phylum", "pct", "Endophyte", "Water", "Harvest")
  missing_phy_base <- setdiff(req_phy_base, names(phy))
  if (length(missing_phy_base) > 0) {
    stop("Missing columns in phylum input: ", paste(missing_phy_base, collapse = ", "))
  }
  
  if (!("Treatment" %in% names(phy)) && !("Amendment" %in% names(phy))) {
    stop("Phylum input must contain either 'Treatment' or 'Amendment'.")
  }
  
  phy <- phy %>%
    dplyr::mutate(
      Treatment_clean = dplyr::case_when(
        "Amendment" %in% names(.) ~ as.character(.data[["Amendment"]]),
        TRUE ~ as.character(.data[["Treatment"]])
      ),
      Treatment_clean = dplyr::recode(Treatment_clean, "VOC" = "6-PP"),
      Endophyte = factor(Endophyte, levels = c("Nil", "AR37")),
      Water     = factor(Water, levels = c("Irrigated", "Drought")),
      Harvest   = factor(Harvest, levels = c("During Drought", "Post Drought")),
      Treatment = factor(Treatment_clean, levels = c("Control", "Polymer", "6-PP")),
      Phylum    = as.character(Phylum)
    ) %>%
    dplyr::filter(!Phylum %in% c("Chordata")) %>%
    dplyr::select(-Treatment_clean)
  
  phy <- phy %>%
    group_by(SampleID) %>%
    dplyr::mutate(pct_rel = pct / sum(pct, na.rm = TRUE) * 100) %>%
    ungroup()
  
  phy_mean0 <- phy %>%
    group_by(Endophyte, Treatment, Water, Harvest, Phylum) %>%
    dplyr::summarise(mean_pct = mean(pct_rel, na.rm = TRUE), .groups = "drop")
  
  top_n_phylum <- 10
  
  top_phyla <- phy_mean0 %>%
    group_by(Phylum) %>%
    dplyr::summarise(overall = mean(mean_pct, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(overall)) %>%
    slice_head(n = top_n_phylum) %>%
    pull(Phylum)
  
  phy_mean <- phy_mean0 %>%
    dplyr::mutate(Phylum = ifelse(Phylum %in% top_phyla, Phylum, "Other")) %>%
    group_by(Endophyte, Treatment, Water, Harvest, Phylum) %>%
    dplyr::summarise(mean_pct = sum(mean_pct), .groups = "drop")
  
  phy_overall_tbl <- phy_mean %>%
    group_by(Phylum) %>%
    dplyr::summarise(overall = mean(mean_pct, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(desc(overall))
  
  phy_named_big_to_small <- phy_overall_tbl %>%
    dplyr::filter(Phylum != "Other") %>%
    pull(Phylum)
  
  phy_levels_small_to_big <- c("Other", rev(phy_named_big_to_small))
  phy_levels_small_to_big <- phy_levels_small_to_big[phy_levels_small_to_big %in% unique(phy_mean$Phylum)]
  
  phy_mean <- phy_mean %>%
    dplyr::mutate(Phylum = factor(Phylum, levels = phy_levels_small_to_big))
  
  phy_mean <- phy_mean %>%
    dplyr::mutate(Endophyte_plot = as.character(Endophyte))
  
  phy_spacer_rows <- phy_mean %>%
    distinct(Treatment, Water, Harvest, Phylum) %>%
    dplyr::mutate(Endophyte_plot = " ", mean_pct = 0)
  
  phy_mean <- bind_rows(phy_mean, phy_spacer_rows) %>%
    dplyr::mutate(Endophyte_plot = factor(Endophyte_plot, levels = c("Nil", " ", "AR37")))
  
  phy_levels_fill <- levels(phy_mean$Phylum)
  
  phy_non_other <- phy_overall_tbl %>%
    dplyr::filter(Phylum != "Other") %>%
    pull(Phylum)
  
  phy_non_other <- phy_non_other[phy_non_other %in% phy_levels_fill]
  phy_n_non_other <- length(phy_non_other)
  
  if (phy_n_non_other > 10) {
    stop("More than 10 non-Other phyla. Reduce top_n_phylum.")
  }
  
  phy_npg10 <- ggsci::pal_npg("nrc")(10)
  phy_npg10_rot <- phy_npg10[c(3:10, 1:2)]
  
  phy_pal_non_other <- phy_npg10_rot[seq_len(phy_n_non_other)]
  names(phy_pal_non_other) <- phy_non_other
  
  phy_pal <- phy_pal_non_other
  if ("Other" %in% phy_levels_fill) phy_pal <- c(phy_pal, Other = "grey75")
  
  phy_pal <- phy_pal[phy_levels_fill]
  
  plot_one_treatment_phylum <- function(df, treatment_label, show_y = TRUE) {
    df1 <- df %>% dplyr::filter(Treatment == treatment_label)
    
    ggplot(df1, aes(x = Endophyte_plot, y = mean_pct, fill = Phylum)) +
      geom_col(width = 1.5, color = NA, position = position_stack(reverse = TRUE)) +
      facet_grid(Harvest ~ Water) +
      coord_cartesian(clip = "off") +
      scale_y_continuous(
        limits = c(0, 100),
        breaks = c(0, 25, 50, 75, 100),
        expand = expansion(mult = c(0, 0.02))
      ) +
      scale_x_discrete(
        labels = c(Nil = "Nil", " " = "", AR37 = "AR37"),
        expand = expansion(add = 0.30)
      ) +
      scale_fill_manual(values = phy_pal, drop = FALSE) +
      guides(fill = guide_legend(reverse = TRUE)) +
      labs(
        title = treatment_label,
        x = NULL,
        y = if (show_y) "Mean relative abundance (%)" else NULL,
        fill = NULL
      ) +
      theme_natureish_comp(base_size = 8) +
      theme(
        axis.text.x = element_text(size = 7),
        plot.title = element_text(face = "bold", size = 8, hjust = 0.5),
        legend.position = "none"
      )
  }
  
  p_phy_control <- plot_one_treatment_phylum(phy_mean, "Control", show_y = TRUE)
  p_phy_polymer <- plot_one_treatment_phylum(phy_mean, "Polymer", show_y = FALSE)
  p_phy_6pp     <- plot_one_treatment_phylum(phy_mean, "6-PP", show_y = FALSE)
  
  final_phylum_fig <- (p_phy_control | p_phy_polymer | p_phy_6pp) +
    patchwork::plot_annotation(tag_levels = "A") +
    patchwork::plot_layout(guides = "collect") &
    theme(
      legend.position = "right",
      plot.margin = margin(12, 16, 10, 10)
    )
  
  ggsave(
    filename = file.path(fig_dir, paste0("Fig_Phylum_ABC_Top", top_n_phylum, "_BigOnTop_Spacer_NPG_OtherGrey.pdf")),
    plot = final_phylum_fig,
    width = 220, height = 105, units = "mm",
    dpi = 600, useDingbats = FALSE, limitsize = FALSE
  )
  
  ggsave(
    filename = file.path(fig_dir, paste0("Fig_Phylum_ABC_Top", top_n_phylum, "_BigOnTop_Spacer_NPG_OtherGrey.png")),
    plot = final_phylum_fig,
    width = 220, height = 105, units = "mm",
    dpi = 600, limitsize = FALSE
  )
  
} else {
  message("Phylum input file not found: ", phylum_fp)
}

## -----------------------------
## 10) Distance matrices
## -----------------------------
dist_bray  <- vegdist(t(rel_abund), method = "bray")
dist_aitch <- dist(t(clr_mat), method = "euclidean")

saveRDS(dist_bray,  file.path(repr_dir, "dist_BrayCurtis.rds"))
saveRDS(dist_aitch, file.path(repr_dir, "dist_Aitchison_CLR_Euclidean.rds"))

## -----------------------------
## 11) Ordination functions
## -----------------------------
ordination_plot <- function(dist_obj, meta_df, method_name, out_prefix) {
  ord <- ape::pcoa(dist_obj)
  scores <- as.data.frame(ord$vectors[, 1:2, drop = FALSE])
  scores$SampleID <- rownames(scores)
  
  df <- scores %>% left_join(meta_df, by = "SampleID")
  var1 <- round(100 * ord$values$Relative_eig[1], 1)
  var2 <- round(100 * ord$values$Relative_eig[2], 1)
  
  g <- ggplot(df, aes(x = Axis.1, y = Axis.2)) +
    geom_point(aes(shape = Endophyte, colour = Water), size = 3) +
    facet_grid(Harvest ~ Treatment) +
    labs(
      title = paste0("PCoA (", method_name, ")"),
      x = paste0("PCoA1 (", var1, "%)"),
      y = paste0("PCoA2 (", var2, "%)")
    ) +
    scale_shape_manual(values = c("Nil" = 16, "AR37" = 17)) +
    scale_colour_manual(values = c("Irrigated" = "#F8766D", "Drought" = "#00BFC4")) +
    guides(
      colour = guide_legend(order = 1),
      shape  = guide_legend(order = 2)
    ) +
    theme_pub(9)
  
  write_csv(df, paste0(out_prefix, "_scores.csv"))
  save_pub(g, paste0(out_prefix, ".png"), 240, 140)
  list(plot = g, scores = df, ord = ord)
}

nmds_plot <- function(dist_obj, meta_df, method_name, out_prefix, trymax = 100) {
  set.seed(1)
  nmds <- vegan::metaMDS(dist_obj, k = 2, trymax = trymax, autotransform = FALSE, trace = FALSE)
  scr <- as.data.frame(scores(nmds, display = "sites"))
  scr$SampleID <- rownames(scr)
  
  df <- scr %>% left_join(meta_df, by = "SampleID")
  
  g <- ggplot(df, aes(x = NMDS1, y = NMDS2)) +
    geom_point(aes(shape = Endophyte, colour = Water), size = 3) +
    facet_grid(Harvest ~ Treatment) +
    labs(
      title = paste0("NMDS (", method_name, ", stress = ", round(nmds$stress, 3), ")"),
      x = "NMDS1",
      y = "NMDS2"
    ) +
    scale_shape_manual(values = c("Nil" = 16, "AR37" = 17)) +
    scale_colour_manual(values = c("Irrigated" = "#F8766D", "Drought" = "#00BFC4")) +
    guides(
      colour = guide_legend(order = 1),
      shape  = guide_legend(order = 2)
    ) +
    theme_pub(9)
  
  write_csv(df, paste0(out_prefix, "_scores.csv"))
  saveRDS(nmds, paste0(out_prefix, ".rds"))
  save_pub(g, paste0(out_prefix, ".png"), 240, 140)
  list(plot = g, scores = df, nmds = nmds)
}

ord_bray  <- ordination_plot(dist_bray,  meta_use, "Bray–Curtis", file.path(ord_dir, "PCoA_BrayCurtis"))
ord_aitch <- ordination_plot(dist_aitch, meta_use, "Aitchison",   file.path(ord_dir, "PCoA_Aitchison"))

nmds_bray  <- nmds_plot(dist_bray,  meta_use, "Bray–Curtis", file.path(ord_dir, "NMDS_BrayCurtis"))
nmds_aitch <- nmds_plot(dist_aitch, meta_use, "Aitchison",   file.path(ord_dir, "NMDS_Aitchison"))

## -----------------------------
## 12) PERMANOVA + pairwise PERMANOVA
## -----------------------------
perm_rhs <- paste(
  "Endophyte +",
  "Treatment +",
  "Water +",
  "Harvest +",
  "Endophyte:Treatment +",
  "Endophyte:Water +",
  "Endophyte:Harvest +",
  "Treatment:Water +",
  "Treatment:Harvest +",
  "Water:Harvest"
)

run_permanova <- function(dist_obj, meta_df, label, rhs = perm_rhs, permutations = 999) {
  
  set.seed(1)
  
  fml <- as.formula(paste0("dist_obj ~ ", rhs))
  
  ad <- adonis2(
    fml,
    data = meta_df,
    permutations = permutations,
    by = "terms"
  )
  
  out <- as.data.frame(ad) %>%
    rownames_to_column("Term")
  
  write_csv(
    out,
    file.path(perm_dir, paste0("PERMANOVA_", label, ".csv"))
  )
  
  out
}

pairwise_permanova <- function(dist_obj, meta_df, group_var, label, permutations = 999) {
  levs <- levels(droplevels(meta_df[[group_var]]))
  combs <- combn(levs, 2, simplify = FALSE)
  
  res_list <- lapply(combs, function(x) {
    sub_meta <- meta_df %>%
      dplyr::filter(.data[[group_var]] %in% x) %>%
      dplyr::mutate(GroupTmp = droplevels(factor(.data[[group_var]])))
    
    sub_ids <- sub_meta$SampleID
    sub_dist <- as.dist(as.matrix(dist_obj)[sub_ids, sub_ids])
    
    set.seed(1)
    ad <- adonis2(sub_dist ~ GroupTmp, data = sub_meta, permutations = permutations)
    
    tibble(
      Variable = group_var,
      Contrast = paste(x, collapse = " vs "),
      Df = ad$Df[1],
      R2 = ad$R2[1],
      F = ad$F[1],
      P = ad$`Pr(>F)`[1]
    )
  })
  
  out <- bind_rows(res_list) %>%
    dplyr::mutate(FDR_BH = p.adjust(P, method = "BH"))
  
  write_csv(
    out,
    file.path(perm_dir, paste0("Pairwise_PERMANOVA_", label, "_", group_var, ".csv"))
  )
  
  out
}

perm_bray  <- run_permanova(dist_bray,  meta_use, "BrayCurtis")
perm_aitch <- run_permanova(dist_aitch, meta_use, "Aitchison")

pw_bray_water     <- pairwise_permanova(dist_bray,  meta_use, "Water", "BrayCurtis")
pw_bray_treatment <- pairwise_permanova(dist_bray,  meta_use, "Treatment", "BrayCurtis")
pw_bray_endo      <- pairwise_permanova(dist_bray,  meta_use, "Endophyte", "BrayCurtis")
pw_bray_harv      <- pairwise_permanova(dist_bray,  meta_use, "Harvest", "BrayCurtis")

pw_aitch_water     <- pairwise_permanova(dist_aitch, meta_use, "Water", "Aitchison")
pw_aitch_treatment <- pairwise_permanova(dist_aitch, meta_use, "Treatment", "Aitchison")
pw_aitch_endo      <- pairwise_permanova(dist_aitch, meta_use, "Endophyte", "Aitchison")
pw_aitch_harv      <- pairwise_permanova(dist_aitch, meta_use, "Harvest", "Aitchison")

## -----------------------------
## 13) Dispersion checks
## -----------------------------
run_betadisper <- function(dist_obj, meta_df, group_var, label) {
  grp <- droplevels(meta_df[[group_var]])
  bd  <- betadisper(dist_obj, grp)
  an  <- as.data.frame(anova(bd)) %>% rownames_to_column("Term")
  pe  <- as.data.frame(permutest(bd, permutations = 999)$tab) %>% rownames_to_column("Term")
  
  write_csv(an, file.path(perm_dir, paste0("Betadisper_ANOVA_", label, "_", group_var, ".csv")))
  write_csv(pe, file.path(perm_dir, paste0("Betadisper_PERM_", label, "_", group_var, ".csv")))
  invisible(list(anova = an, perm = pe))
}

bd_bray_water     <- run_betadisper(dist_bray,  meta_use, "Water", "BrayCurtis")
bd_bray_treatment <- run_betadisper(dist_bray,  meta_use, "Treatment", "BrayCurtis")
bd_bray_endo      <- run_betadisper(dist_bray,  meta_use, "Endophyte", "BrayCurtis")
bd_bray_harv      <- run_betadisper(dist_bray,  meta_use, "Harvest", "BrayCurtis")

bd_aitch_water     <- run_betadisper(dist_aitch, meta_use, "Water", "Aitchison")
bd_aitch_treatment <- run_betadisper(dist_aitch, meta_use, "Treatment", "Aitchison")
bd_aitch_endo      <- run_betadisper(dist_aitch, meta_use, "Endophyte", "Aitchison")
bd_aitch_harv      <- run_betadisper(dist_aitch, meta_use, "Harvest", "Aitchison")

## -----------------------------
## 14) Technical-effect check
## -----------------------------
if (has_flowcell) {
  tech_perm_bray  <- adonis2(dist_bray ~ Flowcell, data = meta_use, permutations = nperm)
  tech_perm_aitch <- adonis2(dist_aitch ~ Flowcell, data = meta_use, permutations = nperm)
  
  write_csv(
    as.data.frame(tech_perm_bray) %>% rownames_to_column("Term"),
    file.path(perm_dir, "TechnicalCheck_Flowcell_BrayCurtis.csv")
  )
  write_csv(
    as.data.frame(tech_perm_aitch) %>% rownames_to_column("Term"),
    file.path(perm_dir, "TechnicalCheck_Flowcell_Aitchison.csv")
  )
}

if (has_batch) {
  tech_perm_bray_batch  <- adonis2(dist_bray ~ Batch, data = meta_use, permutations = nperm)
  tech_perm_aitch_batch <- adonis2(dist_aitch ~ Batch, data = meta_use, permutations = nperm)
  
  write_csv(
    as.data.frame(tech_perm_bray_batch) %>% rownames_to_column("Term"),
    file.path(perm_dir, "TechnicalCheck_Batch_BrayCurtis.csv")
  )
  write_csv(
    as.data.frame(tech_perm_aitch_batch) %>% rownames_to_column("Term"),
    file.path(perm_dir, "TechnicalCheck_Batch_Aitchison.csv")
  )
}

## -----------------------------
## 15) ALDEx2 differential abundance
## -----------------------------
counts <- as.matrix(genus_mat_filt)
storage.mode(counts) <- "numeric"

run_aldex2_contrast <- function(count_mat, meta_df, group_col, level_a, level_b, label,
                                mc = aldex_mc_samples, denom = "all") {
  
  keep <- meta_df[[group_col]] %in% c(level_a, level_b)
  meta2 <- meta_df[keep, , drop = FALSE]
  mat2 <- count_mat[, meta2$SampleID, drop = FALSE]
  storage.mode(mat2) <- "numeric"
  
  conds <- as.character(meta2[[group_col]])
  
  set.seed(1)
  clr <- aldex.clr(mat2, conds = conds, mc.samples = mc, denom = denom, verbose = FALSE)
  tt  <- aldex.ttest(clr, verbose = FALSE)
  ef  <- aldex.effect(clr, verbose = FALSE)
  
  res <- data.frame(
    Genus = rownames(tt),
    we.ep = tt$we.ep,
    we.eBH = p.adjust(tt$we.ep, method = "BH"),
    effect = ef$effect,
    overlap = ef$overlap,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::arrange(we.eBH, desc(abs(effect)))
  
  safe_a <- str_replace_all(level_a, "[^A-Za-z0-9]+", "_")
  safe_b <- str_replace_all(level_b, "[^A-Za-z0-9]+", "_")
  
  out_fp <- file.path(
    aldex_dir,
    paste0("ALDEx2_", label, "_", group_col, "_", safe_a, "_vs_", safe_b, ".csv")
  )
  write_csv(res, out_fp)
  res
}

aldex_water   <- run_aldex2_contrast(counts, meta_use, "Water", "Irrigated", "Drought", "Genus")
aldex_endo    <- run_aldex2_contrast(counts, meta_use, "Endophyte", "Nil", "AR37", "Genus")
aldex_harvest <- run_aldex2_contrast(counts, meta_use, "Harvest", "During Drought", "Post Drought", "Genus")
aldex_6pp     <- run_aldex2_contrast(counts, meta_use, "Treatment", "Control", "6-PP", "Genus")
aldex_poly    <- run_aldex2_contrast(counts, meta_use, "Treatment", "Control", "Polymer", "Genus")
aldex_6pp_poly <- run_aldex2_contrast(counts, meta_use, "Treatment", "Polymer", "6-PP", "Genus")

run_aldex2_stratified <- function(count_mat, meta_df, stratum_col, stratum_level, group_col, a, b, label) {
  sub_meta <- meta_df %>% dplyr::filter(.data[[stratum_col]] == stratum_level)
  run_aldex2_contrast(count_mat, sub_meta, group_col, a, b,
                      paste0(label, "_within_", stratum_col, "_", str_replace_all(stratum_level, "[^A-Za-z0-9]+", "_")))
}

aldex_water_during <- run_aldex2_stratified(
  counts, meta_use, "Harvest", "During Drought",
  "Water", "Irrigated", "Drought", "Genus"
)

aldex_water_post <- run_aldex2_stratified(
  counts, meta_use, "Harvest", "Post Drought",
  "Water", "Irrigated", "Drought", "Genus"
)

sig_summary <- function(res, fdr = fdr_cutoff) {
  res %>%
    dplyr::filter(we.eBH <= fdr) %>%
    dplyr::select(Genus, we.eBH, effect, overlap)
}

write_csv(sig_summary(aldex_water),   file.path(aldex_dir, "ALDEx2_SIGNIFICANT_Water.csv"))
write_csv(sig_summary(aldex_endo),    file.path(aldex_dir, "ALDEx2_SIGNIFICANT_Endophyte.csv"))
write_csv(sig_summary(aldex_harvest), file.path(aldex_dir, "ALDEx2_SIGNIFICANT_Harvest.csv"))
write_csv(sig_summary(aldex_6pp),     file.path(aldex_dir, "ALDEx2_SIGNIFICANT_Control_vs_6PP_plus_polymer.csv"))
write_csv(sig_summary(aldex_poly),    file.path(aldex_dir, "ALDEx2_SIGNIFICANT_Control_vs_Polymer.csv"))
write_csv(sig_summary(aldex_6pp_poly), file.path(aldex_dir, "ALDEx2_SIGNIFICANT_Polymer_vs_6PP_plus_polymer.csv"))

sig_all <- bind_rows(
  sig_summary(aldex_water)        %>% dplyr::mutate(Contrast = "Water: Drought vs Irrigated"),
  sig_summary(aldex_endo)         %>% dplyr::mutate(Contrast = "Endophyte: AR37 vs Nil"),
  sig_summary(aldex_harvest)      %>% dplyr::mutate(Contrast = "Harvest: Post vs During"),
  sig_summary(aldex_6pp)          %>% dplyr::mutate(Contrast = "Treatment: 6-PP vs Control"),
  sig_summary(aldex_poly)         %>% dplyr::mutate(Contrast = "Treatment: Polymer vs Control"),
  sig_summary(aldex_6pp_poly)     %>% dplyr::mutate(Contrast = "Treatment: 6-PP vs Polymer"),
  sig_summary(aldex_water_during) %>% dplyr::mutate(Contrast = "Water within During Drought"),
  sig_summary(aldex_water_post)   %>% dplyr::mutate(Contrast = "Water within Post Drought")
) %>%
  dplyr::select(Contrast, everything())

write_csv(sig_all, file.path(aldex_dir, "ALDEx2_SIGNIFICANT_ALL_CONTRASTS.csv"))

## -----------------------------
## 16) Volcano plots
## -----------------------------
volcano_plot <- function(res, title, out_prefix, fdr = fdr_cutoff, label_top = 10) {
  df <- res %>%
    dplyr::mutate(
      neglog10FDR = -log10(pmax(we.eBH, 1e-300)),
      Significant = we.eBH <= fdr
    ) %>%
    dplyr::arrange(we.eBH)
  
  lab_df <- df %>% dplyr::filter(Significant) %>% slice_head(n = label_top)
  
  g <- ggplot(df, aes(x = effect, y = neglog10FDR)) +
    geom_hline(yintercept = -log10(fdr), linetype = "dashed", linewidth = 0.35) +
    geom_point(aes(fill = Significant), shape = 21, color = "black", size = 2, stroke = 0.2) +
    geom_text_repel(
      data = lab_df,
      aes(label = Genus),
      size = 2.7,
      max.overlaps = 50,
      box.padding = 0.35,
      point.padding = 0.2,
      segment.color = "grey50"
    ) +
    scale_fill_manual(values = c("FALSE" = "grey80", "TRUE" = "black")) +
    labs(
      title = title,
      x = "ALDEx2 effect size",
      y = expression(-log[10]("FDR"))
    ) +
    theme_pub(9) +
    theme(legend.position = "none")
  
  save_pub(g, paste0(out_prefix, ".png"), 95, 85)
  g
}

volcano_plot(aldex_water,   "Drought vs Irrigated",      file.path(aldex_dir, "Volcano_Water"))
volcano_plot(aldex_endo,    "AR37 vs Nil",               file.path(aldex_dir, "Volcano_Endophyte"))
volcano_plot(aldex_harvest, "Post vs During drought",    file.path(aldex_dir, "Volcano_Harvest"))
volcano_plot(aldex_6pp,     "6-PP vs Control", file.path(aldex_dir, "Volcano_6PP_plus_polymer_vs_Control"))
volcano_plot(aldex_poly,    "Polymer vs Control",        file.path(aldex_dir, "Volcano_Polymer_vs_Control"))
volcano_plot(aldex_6pp_poly,"6-PP vs Polymer", file.path(aldex_dir, "Volcano_6PP_plus_polymer_vs_Polymer"))

## -----------------------------
## 17) Heatmap of top variable genera
## -----------------------------
genus_var <- apply(clr_mat, 1, var)
top_heat_genera <- names(sort(genus_var, decreasing = TRUE))[seq_len(min(top_n_heatmap, length(genus_var)))]
heat_mat <- clr_mat[top_heat_genera, , drop = FALSE]

annotation_col <- meta_use %>%
  dplyr::select(SampleID, Endophyte, Treatment, Water, Harvest) %>%
  as.data.frame()

rownames(annotation_col) <- annotation_col$SampleID
annotation_col$SampleID <- NULL

heat_mat_scaled <- t(scale(t(heat_mat)))
heat_mat_scaled[is.na(heat_mat_scaled)] <- 0

pheatmap(
  mat = heat_mat_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  fontsize_row = 7,
  border_color = NA,
  filename = file.path(heat_dir, "Heatmap_top_variable_genera.pdf"),
  width = 10,
  height = 9
)

png(file.path(heat_dir, "Heatmap_top_variable_genera.png"),
    width = 2600, height = 2200, res = 300)
pheatmap(
  mat = heat_mat_scaled,
  annotation_col = annotation_col,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_colnames = FALSE,
  fontsize_row = 7,
  border_color = NA
)
dev.off()

write_csv(
  tibble(Genus = top_heat_genera, Variance = genus_var[top_heat_genera]),
  file.path(heat_dir, "Heatmap_top_variable_genera_list.csv")
)

## -----------------------------
## 18) Save key summary tables
## -----------------------------
summary_counts <- tibble(
  Metric = c(
    "Samples before filtering",
    "Samples after filtering",
    "Genera before filtering",
    "Genera after filtering"
  ),
  Value = c(
    nrow(qc_samples),
    ncol(genus_mat_filt),
    nrow(genus_mat),
    nrow(genus_mat_filt)
  )
)

write_csv(summary_counts, file.path(repr_dir, "Summary_counts.csv"))

## -----------------------------
## 19) Reproducibility log
## -----------------------------
log_fp <- file.path(repr_dir, "Reproducibility_log.txt")
sink(log_fp)

cat("Metatranscriptomics genus-level analysis pipeline\n")
cat("Run date/time: ", as.character(Sys.time()), "\n\n")

cat("Input files:\n")
cat("  Genus table: ", genus_fp, "\n")
cat("  Metadata:    ", meta_fp, "\n\n")

cat("Parameters:\n")
cat("  min_total_counts_per_sample: ", min_total_counts_per_sample, "\n")
cat("  min_prevalence:              ", min_prevalence, "\n")
cat("  min_total_counts_per_taxon:  ", min_total_counts_per_taxon, "\n")
cat("  pseudocount:                 ", pseudocount, "\n")
cat("  ALDEx2 mc.samples:           ", aldex_mc_samples, "\n")
cat("  FDR cutoff:                  ", fdr_cutoff, "\n")
cat("  PERMANOVA permutations:      ", nperm, "\n")
cat("  PERMANOVA rhs:               ", perm_rhs, "\n\n")

cat("Retained samples: ", ncol(genus_mat_filt), "\n")
cat("Retained genera:  ", nrow(genus_mat_filt), "\n\n")

cat("Metadata columns:\n")
print(names(meta_use))
cat("\n")

cat("Session info:\n")
print(sessionInfo())

sink()

message("Analysis completed successfully.")
message("Outputs saved to: ", out_dir)

# =========================================================
# Figure 3 — Nature Communications taxonomy figure
# A = Bray-Curtis PCoA
# B = Shannon diversity
# C = Genus-level composition
# D = PERMANOVA variance explained
# =========================================================

suppressPackageStartupMessages({
  library(cowplot)
})

fig3_dir <- file.path(project_dir, "7.manuscript_figure", "4.figure3")
dir.create(fig3_dir, recursive = TRUE, showWarnings = FALSE)

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

# Panel A — Bray-Curtis PCoA
pcoa_scores <- ord_bray$scores
pcoa_var1 <- round(100 * ord_bray$ord$values$Relative_eig[1], 1)
pcoa_var2 <- round(100 * ord_bray$ord$values$Relative_eig[2], 1)

panel_a <- ggplot(
  pcoa_scores,
  aes(x = Axis.1, y = Axis.2, colour = Treatment, shape = Water)
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
    subtitle = "Genus-level Bray–Curtis PCoA",
    x = paste0("PCoA1 (", pcoa_var1, "%)"),
    y = paste0("PCoA2 (", pcoa_var2, "%)"),
    colour = "Seed treatment",
    shape = "Water"
  ) +
  theme_natcom(15)

# Panel B — Shannon diversity
panel_b <- ggplot(
  alpha_df,
  aes(x = Treatment, y = Shannon, fill = Water)
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.55,
    outlier.shape = NA,
    colour = "black",
    linewidth = 0.45,
    alpha = 0.85
  ) +
  geom_point(
    aes(shape = Endophyte),
    position = position_jitterdodge(
      jitter.width = 0.08,
      dodge.width = 0.75
    ),
    colour = "black",
    size = 1.7,
    alpha = 0.65
  ) +
  facet_wrap(~ Harvest, nrow = 1) +
  scale_fill_manual(values = c(
    "Irrigated" = "#F8766D",
    "Drought" = "#00BFC4"
  )) +
  scale_shape_manual(values = c("Nil" = 16, "AR37" = 17)) +
  labs(
    title = "B",
    subtitle = "Genus-level Shannon diversity",
    x = "Seed treatment",
    y = "Shannon diversity",
    fill = "Water",
    shape = "Endophyte"
  ) +
  theme_natcom(15) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

# Panel C — Genus-level composition with top genera only

top_n_fig3_genera <- 20

top_genera_fig3 <- genus_long %>%
  group_by(Genus) %>%
  summarise(overall_pct = mean(pct_rel, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall_pct)) %>%
  slice_head(n = top_n_fig3_genera) %>%
  pull(Genus)

genus_fig3 <- genus_long %>%
  mutate(Genus2 = ifelse(Genus %in% top_genera_fig3, Genus, "Other")) %>%
  group_by(Treatment, Water, Harvest, Endophyte, Genus2) %>%
  summarise(mean_pct = mean(pct_rel, na.rm = TRUE), .groups = "drop") %>%
  group_by(Treatment, Water, Harvest, Endophyte) %>%
  mutate(mean_pct = mean_pct / sum(mean_pct, na.rm = TRUE) * 100) %>%
  ungroup() %>%
  mutate(Endophyte_plot = as.character(Endophyte))

genus_fig3_order <- genus_fig3 %>%
  group_by(Genus2) %>%
  summarise(overall = mean(mean_pct, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(overall)) %>%
  pull(Genus2)

genus_fig3_order <- c(setdiff(genus_fig3_order, "Other"), "Other")

genus_fig3 <- genus_fig3 %>%
  mutate(Genus2 = factor(Genus2, levels = genus_fig3_order))

genus_fig3_pal <- gen_pal[levels(genus_fig3$Genus2)]

if (any(is.na(genus_fig3_pal))) {
  genus_fig3_pal[is.na(genus_fig3_pal)] <- grDevices::hcl.colors(
    sum(is.na(genus_fig3_pal)),
    palette = "Dark 3"
  )
}

genus_fig3_pal["Other"] <- "grey70"

panel_c <- ggplot(
  genus_fig3,
  aes(x = Endophyte_plot, y = mean_pct, fill = Genus2)
) +
  geom_col(width = 0.70, colour = NA) +
  facet_grid(Harvest + Water ~ Treatment) +
  scale_y_continuous(
    limits = c(0, 100.5),
    breaks = c(0, 25, 50, 75, 100),
    expand = expansion(mult = c(0, 0.02))
  ) +
  scale_fill_manual(values = genus_fig3_pal, drop = FALSE) +
  guides(
    fill = guide_legend(
      ncol = 1,
      byrow = TRUE,
      override.aes = list(colour = NA)
    )
  ) +
  labs(
    title = "C",
    subtitle = "Dominant genus-level transcript profiles",
    x = NULL,
    y = "Mean relative abundance (%)",
    fill = "Genus"
  ) +
  theme_natcom(12) +
  theme(
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.title.y = element_text(size = 12, face = "bold"),
    strip.text = element_text(size = 10, face = "bold"),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9, face = "bold"),
    legend.key.height = unit(2.7, "mm"),
    legend.key.width = unit(2.7, "mm"),
    plot.subtitle = element_text(size = 12, face = "bold"),
    panel.spacing = unit(3, "mm")
  )

# Panel D — PERMANOVA variance explained using full model + residual

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
    x == "Endophyte:Treatment" ~ "Endophyte × Seed treatment",
    x == "Endophyte:Water" ~ "Endophyte × Water regime",
    x == "Endophyte:Harvest" ~ "Endophyte × Harvest stage",
    x == "Treatment:Water" ~ "Seed treatment × Water regime",
    x == "Treatment:Harvest" ~ "Seed treatment × Harvest stage",
    x == "Water:Harvest" ~ "Water regime × Harvest stage",
    x == "Residual" ~ "Unexplained variance",
    TRUE ~ x
  )
}

perm_plot_df <- perm_bray %>%
  dplyr::rename(p_value = `Pr(>F)`) %>%
  dplyr::mutate(
    Term = trimws(as.character(Term)),
    Term_clean = clean_perm_term(Term),
    Analysis = dplyr::case_when(
      Term %in% c("Endophyte", "Treatment", "Water", "Harvest") ~ "Main effect",
      grepl(":", Term, fixed = TRUE) ~ "Interaction",
      Term == "Residual" ~ "Residual",
      TRUE ~ "Other"
    ),
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
    )
  ) %>%
  dplyr::filter(Analysis %in% c("Main effect", "Interaction", "Residual")) %>%
  dplyr::arrange(R2) %>%
  dplyr::mutate(
    Term_clean = factor(Term_clean, levels = Term_clean),
    significant = factor(
      significant,
      levels = c("Significant", "Not significant", "Unexplained variance")
    )
  )

panel_d <- ggplot(
  perm_plot_df,
  aes(x = Term_clean, y = R2, fill = significant)
) +
  geom_col(width = 0.85) +
  geom_text(
    aes(label = label),
    hjust = -0.15,
    size = 3.7
  ) +
  coord_flip(clip = "off") +
  scale_fill_manual(
    values = c(
      "Significant" = "#8A3FFC",
      "Not significant" = "grey70",
      "Unexplained variance" = "grey35"
    ),
    drop = FALSE,
    name = NULL
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.18))
  ) +
  labs(
    title = "D",
    subtitle = "Variance explained by PERMANOVA model terms",
    x = NULL,
    y = expression(R^2~"(variance explained)")
  ) +
  theme_natcom(13) +
  theme(
    axis.text.y = element_text(size = 10, colour = "black"),
    axis.text.x = element_text(size = 10, colour = "black"),
    axis.title.x = element_text(size = 12, face = "bold"),
    legend.position = "right",
    legend.text = element_text(size = 9),
    plot.subtitle = element_text(size = 12, face = "bold"),
    plot.margin = margin(8, 35, 8, 8)
  )

# Assemble Figure 3
top <- plot_grid(
  panel_a,
  panel_b,
  nrow = 1,
  rel_widths = c(1.1, 1)
)

bottom <- plot_grid(
  panel_c,
  panel_d,
  nrow = 1,
  rel_widths = c(0.85, 1.15)
)

fig3 <- plot_grid(
  top,
  bottom,
  ncol = 1,
  rel_heights = c(0.95, 1.25)
)

ggsave(
  file.path(fig3_dir, "Figure3_NatCom_taxonomy.pdf"),
  fig3,
  width = 20,
  height = 13,
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  file.path(fig3_dir, "Figure3_NatCom_taxonomy.png"),
  fig3,
  width = 20,
  height = 13,
  dpi = 600,
  limitsize = FALSE
)

message("Figure 3 saved to: ", fig3_dir)

