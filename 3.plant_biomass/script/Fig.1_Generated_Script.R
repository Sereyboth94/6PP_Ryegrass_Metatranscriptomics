###############################################################################
# Figure 1 — Nature Communications
# A fungal volatile enhances drought resilience in perennial ryegrass
# through host and microbial modulation
#
# Outputs:
#   Figure1_NatCom.pdf
#   Figure1_NatCom.png
#   Figure1_PanelB_SPAD_Drought.csv
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

project_root <- "C:/Users/soths/OneDrive - Lincoln University/Writing/NatMicro/github"

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

theme_natcom <- function(base_size = 14) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title = element_text(face = "bold", colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.text.x = element_text(angle = 25, hjust = 1),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      strip.text = element_text(face = "bold", colour = "black"),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(colour = "black"),
      panel.grid = element_blank(),
      plot.margin = margin(8, 8, 8, 8)
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
    x = -0.10,
    y = -0.08,
    width = 1.20,
    height = 1.20
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

make_spad_letters <- function(df_sub) {
  mod <- aov(SPAD ~ Treatment * Endophyte, data = df_sub)
  
  emmeans(mod, ~ Treatment * Endophyte) %>%
    multcomp::cld(Letters = letters, adjust = "tukey") %>%
    as.data.frame() %>%
    as_tibble() %>%
    mutate(Letters = stringr::str_trim(.group)) %>%
    dplyr::select(Treatment, Endophyte, Letters)
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
  group_by(Timepoint, Water) %>%
  group_modify(~ make_spad_letters(.x)) %>%
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

panel_b <- spad_plot_df %>%
  filter(Water == "Drought") %>%
  ggplot(aes(x = Treatment, y = Mean, fill = Endophyte)) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.70,
    colour = "black",
    linewidth = 0.4
  ) +
  geom_errorbar(
    aes(ymin = Mean - SD, ymax = Mean + SD),
    position = position_dodge(width = 0.8),
    width = 0.15,
    linewidth = 0.5
  ) +
  geom_text(
    aes(y = label_y, label = Letters),
    position = position_dodge(width = 0.8),
    size = 3.4,
    vjust = 0
  ) +
  facet_wrap(~ Timepoint, nrow = 1) +
  scale_fill_manual(values = endo_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.20))) +
  labs(
    x = "Seed treatment",
    y = "SPAD",
    fill = "Endophyte"
  ) +
  theme_natcom(12)

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
      width = 0.58,
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
      size = 1.5,
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
      size = 3.6,
      vjust = 0
    ) +
    facet_wrap(
      ~ Harvest,
      nrow = 1,
      scales = "free_y"
    ) +
    scale_fill_manual(values = endo_cols) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.20))
    ) +
    labs(
      x = "Seed treatment",
      y = y_lab,
      fill = "Endophyte"
    ) +
    theme_natcom(14) +
    theme(
      legend.position = "right",
      axis.text.x = element_text(
        angle = 25,
        hjust = 1,
        size = 12
      ),
      axis.title = element_text(
        size = 14,
        face = "bold"
      ),
      strip.text = element_text(
        size = 12,
        face = "bold"
      ),
      legend.title = element_text(
        size = 13,
        face = "bold"
      ),
      legend.text = element_text(
        size = 12
      )
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
    plot.title = element_text(face = "bold", size = 18, hjust = 0)
  )

panel_b <- panel_b +
  ggtitle("B") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0)
  )

panel_c <- panel_c +
  ggtitle("C") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0)
  )

panel_d <- panel_d +
  ggtitle("D") +
  theme(
    plot.title = element_text(face = "bold", size = 18, hjust = 0)
  )

fig1 <- (panel_a | panel_b) /
  (panel_c | panel_d) +
  plot_layout(
    widths = c(1.4, 1.6),
    heights = c(1.2, 1)
  )

ggsave(
  filename = file.path(manuscript_dir, "Figure1_NatCom.pdf"),
  plot = fig1,
  width = 16,
  height = 10,
  useDingbats = FALSE,
  limitsize = FALSE
)

ggsave(
  filename = file.path(manuscript_dir, "Figure1_NatCom.png"),
  plot = fig1,
  width = 16,
  height = 10,
  dpi = 600,
  limitsize = FALSE
)

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "sessionInfo_Figure1_NatCom.txt")
)

message("Done. Figure 1 saved to: ", results_dir)

