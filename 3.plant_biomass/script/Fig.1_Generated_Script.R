# ============================================================================
# Larger labels for the final, assembled manuscript figure
# Layout revision: 2026-10-12 (wider exports, explicit ticks and edge padding)
# Text sizes retain the previous export calibration as the canvas widens.
# The PDF, PNG and LZW-compressed TIFF use the same layout; rasters use 600 dpi.
# Change MANUSCRIPT_FONT_SCALE to 1.10 for a further 10% increase.
# MANUSCRIPT_PROJECT_ROOT optionally overrides the original Windows project path.
# ============================================================================
FIGURE_DISPLAY_WIDTH_IN <- as.numeric(Sys.getenv("MANUSCRIPT_FIGURE_WIDTH_IN", "6.27"))
FIGURE_EXPORT_WIDTH_IN <- as.numeric(Sys.getenv("MANUSCRIPT_EXPORT_WIDTH_IN", "8.0"))
# Keep the previous exported text and symbol sizes as the canvas widens.
FIGURE_REFERENCE_WIDTH_IN <- 7.5
FIGURE_EXPORT_HEIGHT_IN <- 5.0
FIGURE_REVISION <- "2026-10-12 wider-export update"
message("Figure 1: ", FIGURE_REVISION)
FIGURE_FONT_SCALE <- as.numeric(Sys.getenv("MANUSCRIPT_FONT_SCALE", "1.0"))
# Reduce points and line keys by 20%; text remains independently sized.
FIGURE_SYMBOL_SCALE <- as.numeric(Sys.getenv("MANUSCRIPT_SYMBOL_SCALE", "0.80"))

if (any(!is.finite(c(FIGURE_EXPORT_WIDTH_IN, FIGURE_DISPLAY_WIDTH_IN, FIGURE_FONT_SCALE, FIGURE_SYMBOL_SCALE))) ||
    any(c(FIGURE_EXPORT_WIDTH_IN, FIGURE_DISPLAY_WIDTH_IN, FIGURE_FONT_SCALE, FIGURE_SYMBOL_SCALE) <= 0)) {
  stop("Export width, manuscript width, font scale and symbol scale must be positive numbers.")
}

FIGURE_TEXT <- list(
  axis = 9, axis_title = 9.5, facet = 8.5,
  legend = 8, legend_title = 8.5,
  subtitle = 9, panel = 12, annotation = 8.5, dense = 8
)

# Themes use points; geom_text/geom_label use millimetres.
figure_pt <- function(pt) {
  pt * FIGURE_FONT_SCALE * FIGURE_REFERENCE_WIDTH_IN / FIGURE_DISPLAY_WIDTH_IN
}
figure_mm <- function(pt) figure_pt(pt) / (72.27 / 25.4)
figure_symbol <- function(mm) {
  mm * FIGURE_SYMBOL_SCALE * FIGURE_REFERENCE_WIDTH_IN / FIGURE_DISPLAY_WIDTH_IN
}
figure_subtitle <- function(x, width = 32) stringr::str_wrap(x, width = width)

# Extract a single legend without reserving space inside the data panel.
# Supports both the older guide-box and the newer position-specific guide boxes.
figure_legend <- function(plot_obj, position = "bottom") {
  table <- ggplot2::ggplotGrob(
    plot_obj + ggplot2::theme(legend.position = position)
  )
  index <- which(table$layout$name == paste0("guide-box-", position))
  if (!length(index)) index <- which(table$layout$name == "guide-box")
  if (length(index)) {
    index <- index[!vapply(table$grobs[index], inherits, logical(1), "zeroGrob")]
  }
  if (!length(index)) return(grid::nullGrob())
  table$grobs[[index[1L]]]
}

save_manuscript_figure <- function(plot_obj, out_path, stem) {
  dir.create(out_path, recursive = TRUE, showWarnings = FALSE)
  stem <- paste0(stem, "_LargeLabels")
  pdf_fp <- file.path(out_path, paste0(stem, ".pdf"))

  # Cairo embeds the system fonts when this device is available.
  if (isTRUE(capabilities("cairo"))) {
    ggplot2::ggsave(
      pdf_fp, plot_obj, device = grDevices::cairo_pdf,
      width = FIGURE_EXPORT_WIDTH_IN, height = FIGURE_EXPORT_HEIGHT_IN,
      units = "in", bg = "white", limitsize = FALSE
    )
  } else {
    ggplot2::ggsave(
      pdf_fp, plot_obj, device = grDevices::pdf, useDingbats = FALSE,
      width = FIGURE_EXPORT_WIDTH_IN, height = FIGURE_EXPORT_HEIGHT_IN,
      units = "in", bg = "white", limitsize = FALSE
    )
    message("Cairo is unavailable; used the standard PDF device for ", pdf_fp)
  }

  ggplot2::ggsave(
    file.path(out_path, paste0(stem, ".png")), plot_obj,
    width = FIGURE_EXPORT_WIDTH_IN, height = FIGURE_EXPORT_HEIGHT_IN,
    units = "in", dpi = 600, bg = "white", limitsize = FALSE
  )
  ggplot2::ggsave(
    file.path(out_path, paste0(stem, ".tiff")), plot_obj,
    device = grDevices::tiff, compression = "lzw",
    width = FIGURE_EXPORT_WIDTH_IN, height = FIGURE_EXPORT_HEIGHT_IN,
    units = "in", dpi = 600, bg = "white", limitsize = FALSE
  )
  message("Saved larger-label PDF, PNG and 600-dpi TIFF to: ", out_path)
}

###############################################################################
# Figure 1 — Nature Communications
# A fungal volatile enhances drought resilience in perennial ryegrass
# through host and microbial modulation
#
# Outputs:
#   Figure1_NatCom_LargeLabels.pdf
#   Figure1_NatCom_LargeLabels.png
#   Figure1_NatCom_LargeLabels.tiff
#   Figure1_PanelB_SPAD_Drought.csv
#   Figure1_SPAD_Letters_by_Water.csv
#   Figure1_SPAD_Tukey_by_Water.csv
#   Figure1_Biomass_Drought.csv
###############################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(stringr)
  library(emmeans)
  library(multcomp)
  library(patchwork)
  library(cowplot)
  library(png)
  library(grid)
})

# ─────────────────────────────────────────────────────────────────────────────
# 1. Paths
# ─────────────────────────────────────────────────────────────────────────────

project_root <- Sys.getenv("MANUSCRIPT_PROJECT_ROOT", "C:/Users/soths/OneDrive - Lincoln University/Writing/NatMicro/github")

data_dir <- file.path(project_root, "3.plant_biomass", "data")

results_dir <- file.path(project_root, "3.plant_biomass", "results")

manuscript_dir <- file.path(project_root,
                            "7.manuscript_figure",
                            "2.figure1")

dir.create(results_dir,
           recursive = TRUE,
           showWarnings = FALSE)

dir.create(manuscript_dir,
           recursive = TRUE,
           showWarnings = FALSE)

spad_dir <- data_dir

biomass_dir <- data_dir

design_png <- file.path(data_dir,
                        "Biotron_Experiment.png")

# ─────────────────────────────────────────────────────────────────────────────
# 2. Theme and colours
# ─────────────────────────────────────────────────────────────────────────────

# Panel A uses Biotron_Experiment.png: its own labels must be enlarged in that source image.
theme_natcom <- function(base_size = FIGURE_TEXT$axis) {
  theme_classic(base_size = figure_pt(base_size), base_family = "sans") +
    theme(
      axis.title = element_text(size = figure_pt(FIGURE_TEXT$axis_title),
                                face = "bold", colour = "black"),
      axis.text = element_text(size = figure_pt(FIGURE_TEXT$axis), colour = "black"),
      axis.text.x = element_text(size = figure_pt(FIGURE_TEXT$axis), angle = 35, hjust = 1),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(size = figure_pt(FIGURE_TEXT$facet),
                                face = "bold", colour = "black"),
      legend.title = element_text(size = figure_pt(FIGURE_TEXT$legend_title),
                                  face = "bold"),
      legend.text = element_text(size = figure_pt(FIGURE_TEXT$legend), colour = "black"),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "vertical",
      legend.key.size = unit(figure_pt(8), "pt"),
      legend.spacing.y = unit(figure_pt(2), "pt"),
      legend.box.spacing = unit(figure_pt(2), "pt"),
      plot.title = element_text(size = figure_pt(FIGURE_TEXT$panel), face = "bold", hjust = 0),
      plot.subtitle = element_text(size = figure_pt(FIGURE_TEXT$subtitle),
                                   face = "bold", hjust = 0, lineheight = 1.05,
                                   margin = margin(b = figure_pt(2))),
      plot.title.position = "plot",
      legend.margin = margin(0, 0, 0, 0),
      panel.spacing = unit(figure_pt(4), "pt"),
      plot.margin = margin(figure_pt(3), figure_pt(3), figure_pt(3), figure_pt(3))
    )
}

endo_cols <- c(
  "Nil"  = "#6baed6",
  "AR37" = "#f0a35e"
)

treat_cols <- c(
  "Control" = "#86c1b5",
  "Polymer" = "#e3e59a",
  "6-PP"    = "#b2b0cf"
)

# ─────────────────────────────────────────────────────────────────────────────
# 3. Panel A — experimental design
# ─────────────────────────────────────────────────────────────────────────────

if (!file.exists(design_png)) {
  stop("Experimental design image not found: ", design_png)
}

panel_a <- ggdraw() +
  draw_image(
    png::readPNG(design_png),
    x = 0,
    y = 0,
    width = 1,
    height = 1
  ) +
  theme_void()

# ─────────────────────────────────────────────────────────────────────────────
# 4. SPAD data
# ─────────────────────────────────────────────────────────────────────────────

get_timepoint <- function(path) {
  fn <- basename(path)
  
  case_when(
    str_detect(fn, regex("before[_ ]?drought", ignore_case = TRUE)) ~ "S1",
    str_detect(fn, regex("during[_ ]?drought", ignore_case = TRUE)) ~ "S2",
    str_detect(fn, regex("post[_ ]?drought", ignore_case = TRUE))   ~ "S3",
    TRUE ~ NA_character_
  )
}

read_spad_file <- function(path) {
  tp <- get_timepoint(path)
  
  if (is.na(tp)) {
    stop("Could not determine SPAD timepoint from filename: ", basename(path))
  }
  
  df <- read_excel(path)
  
  names(df) <- names(df) %>%
    tolower() %>%
    gsub("[^a-z0-9]+", "_", .) %>%
    gsub("^_|_$", "", .)
  
  required <- c("treatment", "endophyte", "water", "replicate", "spad")
  missing <- setdiff(required, names(df))
  
  if (length(missing) > 0) {
    stop("Missing columns in ", basename(path), ": ", paste(missing, collapse = ", "))
  }
  
  df %>%
    transmute(
      Timepoint = tp,
      Treatment = str_squish(as.character(treatment)),
      Endophyte = str_squish(as.character(endophyte)),
      Water     = str_squish(as.character(water)),
      Replicate = as.integer(replicate),
      SPAD      = as.numeric(spad)
    )
}

spad_files <- list.files(
  spad_dir,
  pattern = "\\.xlsx$",
  full.names = TRUE,
  ignore.case = TRUE
)

spad_files <- spad_files[
  grepl("spad|before|during|post", basename(spad_files), ignore.case = TRUE) &
    !grepl("fresh|dry|weight|biomass", basename(spad_files), ignore.case = TRUE)
]

if (length(spad_files) == 0) {
  stop("No SPAD .xlsx files found in: ", spad_dir)
}

spad_raw <- bind_rows(lapply(spad_files, read_spad_file)) %>%
  mutate(
    Timepoint = factor(Timepoint, levels = c("S1", "S2", "S3")),
    Treatment = factor(Treatment, levels = c("Control", "Polymer", "6-PP")),
    Endophyte = factor(Endophyte, levels = c("Nil", "AR37")),
    Water = factor(Water, levels = c("Irrigated", "Drought"))
  )

# Fit the same three-way model used for the omnibus ANOVA at each timepoint.
# Compact-letter displays and pairwise comparisons are then calculated as
# simple effects among Treatment x Endophyte combinations within each Water
# regime. Panel B displays the Drought simple effects only.
make_spad_letters <- function(df_sub) {
  mod <- aov(
    SPAD ~ Treatment * Endophyte * Water,
    data = df_sub
  )

  emm <- emmeans(
    mod,
    ~ Treatment * Endophyte | Water
  )

  emm %>%
    multcomp::cld(Letters = letters, adjust = "tukey") %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(Letters = stringr::str_trim(.group)) %>%
    dplyr::select(Water, Treatment, Endophyte, Letters)
}

make_spad_tukey <- function(df_sub) {
  mod <- aov(
    SPAD ~ Treatment * Endophyte * Water,
    data = df_sub
  )

  emm <- emmeans(
    mod,
    ~ Treatment * Endophyte | Water
  )

  pairs(emm, adjust = "tukey") %>%
    summary(infer = c(TRUE, TRUE)) %>%
    as.data.frame() %>%
    as_tibble()
}

spad_summary <- spad_raw %>%
  group_by(Timepoint, Treatment, Endophyte, Water) %>%
  summarise(
    n = n(),
    Mean = mean(SPAD, na.rm = TRUE),
    SD = sd(SPAD, na.rm = TRUE),
    .groups = "drop"
  )

spad_letters <- spad_raw %>%
  group_by(Timepoint) %>%
  group_modify(~ make_spad_letters(.x)) %>%
  ungroup()

spad_tukey <- spad_raw %>%
  group_by(Timepoint) %>%
  group_modify(~ make_spad_tukey(.x)) %>%
  ungroup()

spad_plot_df <- spad_summary %>%
  left_join(
    spad_letters,
    by = c("Timepoint", "Water", "Treatment", "Endophyte")
  ) %>%
  group_by(Timepoint, Water) %>%
  mutate(label_y = Mean + SD + max(Mean + SD, na.rm = TRUE) * 0.05) %>%
  ungroup()

write_csv(
  spad_plot_df %>% filter(Water == "Drought"),
  file.path(results_dir, "Figure1_PanelB_SPAD_Drought.csv")
)

write_csv(
  spad_letters,
  file.path(results_dir, "Figure1_SPAD_Letters_by_Water.csv")
)

write_csv(
  spad_tukey,
  file.path(results_dir, "Figure1_SPAD_Tukey_by_Water.csv")
)

# Stagger the two endophyte letters vertically to prevent multi-letter CLDs touching.
# This changes label placement only; the original means, SDs and CLDs are retained.
spad_letter_step <- max(spad_plot_df$Mean + spad_plot_df$SD, na.rm = TRUE) * 0.10

panel_b <- spad_plot_df %>%
  filter(Water == "Drought") %>%
  ggplot(aes(x = Treatment, y = Mean, fill = Endophyte)) +
  geom_col(
    position = position_dodge(width = 0.90),
    width = 0.56,
    colour = "black",
    linewidth = 0.4
  ) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.90),
    width = 0.15,
    linewidth = 0.5
  ) +
  geom_text(
    aes(y = label_y + ifelse(Endophyte == "AR37", spad_letter_step, 0),
        label = Letters),
    position = position_dodge(width = 0.90),
    size = figure_mm(FIGURE_TEXT$annotation),
    fontface = "plain",
    vjust = 0,
    show.legend = FALSE
  ) +
  facet_wrap(~ Timepoint, nrow = 1) +
  scale_fill_manual(values = endo_cols) +
  scale_y_continuous(n.breaks = 4, expand = expansion(mult = c(0, 0.18))) +
  labs(
    x = "Seed treatment",
    y = "SPAD",
    fill = "Endophyte"
  ) +
  theme_natcom()

# ─────────────────────────────────────────────────────────────────────────────
# 5. Biomass data
# ─────────────────────────────────────────────────────────────────────────────

parse_biomass_file <- function(path) {
  fn <- basename(path)
  
  tibble(
    Metric = case_when(
      str_detect(fn, regex("^fresh[_ ]?weight", ignore_case = TRUE)) ~ "Fresh weight",
      str_detect(fn, regex("^dry[_ ]?weight", ignore_case = TRUE))   ~ "Dry weight",
      TRUE ~ NA_character_
    ),
    Harvest = case_when(
      str_detect(fn, regex("during[_ ]?drought", ignore_case = TRUE)) ~ "During drought",
      str_detect(fn, regex("post[_ ]?drought", ignore_case = TRUE))   ~ "Post drought",
      TRUE ~ NA_character_
    )
  )
}

read_biomass_file <- function(path) {
  meta <- parse_biomass_file(path)
  
  if (is.na(meta$Metric) || is.na(meta$Harvest)) {
    stop("Could not determine biomass metadata from filename: ", basename(path))
  }
  
  df <- read_excel(path)
  
  names(df) <- names(df) %>%
    tolower() %>%
    gsub("[^a-z0-9]+", "_", .) %>%
    gsub("^_|_$", "", .)
  
  value_col <- case_when(
    "freshweight"  %in% names(df) ~ "freshweight",
    "fresh_weight" %in% names(df) ~ "fresh_weight",
    "dryweight"    %in% names(df) ~ "dryweight",
    "dry_weight"   %in% names(df) ~ "dry_weight",
    TRUE ~ NA_character_
  )
  
  if (is.na(value_col)) {
    stop("Could not find biomass value column in: ", basename(path))
  }
  
  df %>%
    transmute(
      Treatment = str_squish(as.character(treatment)),
      Endophyte = str_squish(as.character(endophyte)),
      Water     = str_squish(as.character(water)),
      Replicate = as.integer(replicate),
      Value     = as.numeric(.data[[value_col]])
    ) %>%
    mutate(
      Metric = meta$Metric,
      Harvest = meta$Harvest
    )
}

biomass_files <- list.files(
  biomass_dir,
  pattern = "\\.xlsx$",
  full.names = TRUE,
  ignore.case = TRUE
)

biomass_files <- biomass_files[
  grepl("fresh|dry|weight", basename(biomass_files), ignore.case = TRUE)
]

if (length(biomass_files) == 0) {
  stop("No biomass .xlsx files found in: ", biomass_dir)
}

biomass_raw <- bind_rows(lapply(biomass_files, read_biomass_file)) %>%
  mutate(
    Treatment = recode(
      Treatment,
      "VOC" = "6-PP",
      "6PP" = "6-PP",
      "6PP plus polymer" = "6-PP",
      "6-PP plus polymer" = "6-PP",
      "6PP+polymer" = "6-PP",
      "6-PP+polymer" = "6-PP"
    ),
    Treatment = factor(Treatment, levels = c("Control", "Polymer", "6-PP")),
    Endophyte = factor(Endophyte, levels = c("Nil", "AR37")),
    Water = factor(Water, levels = c("Irrigated", "Drought")),
    Harvest = factor(Harvest, levels = c("During drought", "Post drought")),
    Metric = factor(Metric, levels = c("Fresh weight", "Dry weight"))
  )

if (any(is.na(biomass_raw$Treatment))) stop("Biomass Treatment contains NA. Check treatment naming.")
if (any(is.na(biomass_raw$Endophyte))) stop("Biomass Endophyte contains NA. Check endophyte naming.")
if (any(is.na(biomass_raw$Water))) stop("Biomass Water contains NA. Check water naming.")
if (any(is.na(biomass_raw$Harvest))) stop("Biomass Harvest contains NA. Check file names.")
if (any(is.na(biomass_raw$Metric))) stop("Biomass Metric contains NA. Check file names.")

write_csv(
  biomass_raw %>% filter(Water == "Drought"),
  file.path(results_dir, "Figure1_Biomass_Drought.csv")
)

make_biomass_letters <- function(df_sub) {
  mod <- aov(Value ~ Treatment, data = df_sub)
  
  emmeans(mod, ~ Treatment) %>%
    multcomp::cld(Letters = letters, adjust = "tukey") %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(Letters = stringr::str_trim(.group)) %>%
    dplyr::select(Treatment, Letters)
}

biomass_letters <- biomass_raw %>%
  filter(Water == "Drought") %>%
  group_by(Metric, Harvest, Endophyte) %>%
  group_modify(~ make_biomass_letters(.x)) %>%
  ungroup()

biomass_label_df <- biomass_raw %>%
  filter(Water == "Drought") %>%
  group_by(Metric, Harvest, Endophyte, Treatment) %>%
  summarise(
    y_pos = max(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(
    biomass_letters,
    by = c("Metric", "Harvest", "Endophyte", "Treatment")
  ) %>%
  group_by(Metric, Harvest, Endophyte) %>%
  mutate(
    y_pos = y_pos + max(y_pos, na.rm = TRUE) * 0.07
  ) %>%
  ungroup()

plot_biomass_panel <- function(metric_name, y_lab) {
  
  df <- biomass_raw %>%
    filter(
      Water == "Drought",
      Metric == metric_name
    )
  
  lab <- biomass_label_df %>%
    filter(Metric == metric_name)
  
  ggplot(
    df,
    aes(
      x = Treatment,
      y = Value,
      fill = Endophyte
    )
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.75),
      width = 0.50,
      outlier.shape = NA,
      colour = "black",
      linewidth = 0.45,
      alpha = 0.85
    ) +
    geom_point(
      colour = "black",
      position = position_jitterdodge(
        jitter.width = 0.06,
        dodge.width = 0.75
      ),
      size = figure_symbol(1.2),
      alpha = 0.60,
      show.legend = FALSE
    ) +
    geom_text(
      data = lab,
      aes(
        x = Treatment,
        y = y_pos,
        label = Letters,
        group = Endophyte
      ),
      position = position_dodge(width = 0.75),
      inherit.aes = FALSE,
      size = figure_mm(FIGURE_TEXT$annotation),
      fontface = "plain",
      vjust = 0,
      show.legend = FALSE
    ) +
    facet_wrap(
      ~ Harvest,
      nrow = 1,
      scales = "free_y"
    ) +
    scale_fill_manual(values = endo_cols) +
    scale_y_continuous(
      n.breaks = 4,
      expand = expansion(mult = c(0, 0.20))
    ) +
    labs(
      x = "Seed treatment",
      y = y_lab,
      fill = "Endophyte"
    ) +
    theme_natcom() +
    theme(
      legend.position = "bottom"
    )
}

panel_c <- plot_biomass_panel("Fresh weight", "Fresh weight (g)")
panel_d <- plot_biomass_panel("Dry weight", "Dry weight (g)")

# ─────────────────────────────────────────────────────────────────────────────
# 6. Assemble Figure 1
# ─────────────────────────────────────────────────────────────────────────────

panel_a <- panel_a +
  ggtitle("A") +
  theme(
    plot.title = element_text(face = "bold", size = figure_pt(FIGURE_TEXT$panel), hjust = 0)
  )

panel_b <- panel_b +
  ggtitle("B") +
  theme(
    plot.title = element_text(face = "bold", size = figure_pt(FIGURE_TEXT$panel), hjust = 0)
  )

panel_c <- panel_c +
  ggtitle("C") +
  theme(
    plot.title = element_text(face = "bold", size = figure_pt(FIGURE_TEXT$panel), hjust = 0)
  )

panel_d <- panel_d +
  ggtitle("D") +
  theme(
    plot.title = element_text(face = "bold", size = figure_pt(FIGURE_TEXT$panel), hjust = 0)
  )

# One shared Endophyte legend, using only the bar layer's fill keys.
legend_plot <- panel_b +
  guides(fill = guide_legend(nrow = 1, title.position = "left")) +
  theme(legend.position = "bottom", legend.box = "horizontal")
common_legend <- figure_legend(legend_plot)

# Restore the original two-by-two panel arrangement.
# Slightly more width for B accommodates its three time-point facets.
top_row <- cowplot::plot_grid(
  panel_a,
  panel_b + theme(legend.position = "none"),
  nrow = 1, rel_widths = c(1.25, 1.75)
)
bottom_row <- cowplot::plot_grid(
  panel_c + theme(legend.position = "none"),
  panel_d + theme(legend.position = "none"),
  nrow = 1, rel_widths = c(1, 1)
)
fig1 <- cowplot::plot_grid(
  top_row, bottom_row, common_legend,
  ncol = 1, rel_heights = c(2.40, 2.30, 0.30)
)

save_manuscript_figure(fig1, manuscript_dir, "Figure1_NatCom")



writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "sessionInfo_Figure1_NatCom.txt")
)

message("Done. Figure 1 saved to: ", manuscript_dir)

