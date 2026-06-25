###############################################################################
# Table 1 — Germination and endophyte transmission
# A fungal volatile enhances drought resilience in perennial ryegrass
# through host and microbial modulation
#
# Purpose:
#   Analyse seed germination and AR37 endophyte transmission after seed treatment
#   and generate reproducible outputs for manuscript Table 1.
#
# Expected GitHub project structure:
#   2.germinatio_endophyte/
#     data/
#       Combine.xlsx              # germination data; filename can vary
#       Endophyte.xlsx            # endophyte data; filename can vary
#     results/
#   7.manuscript_figure/
#     1.table1/
#
# Outputs:
#   2.germinatio_endophyte/results/
#     Table1_germination_endophyte_results.xlsx
#     Table1_germination_endophyte_manuscript.csv
#     germination_cleaned.csv
#     endophyte_cleaned.csv
#     sessionInfo_Table1_germination_endophyte.txt
#
#   7.manuscript_figure/1.table1/
#     Table1_germination_endophyte_manuscript.xlsx
#     Table1_germination_endophyte_manuscript.csv
###############################################################################

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(tidyr)
  library(emmeans)
  library(multcomp)
  library(writexl)
  library(readr)
})

# ─────────────────────────────────────────────────────────────────────────────
# 1. Paths
# ─────────────────────────────────────────────────────────────────────────────

project_root <- "C:/Users/soths/OneDrive - Lincoln University/Writing/NatMicro/github"

analysis_dir <- file.path(project_root, "2.germinatio_endophyte")
data_dir     <- file.path(analysis_dir, "data")
results_dir  <- file.path(analysis_dir, "results")

table1_dir   <- file.path(project_root, "7.manuscript_figure", "1.table1")

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table1_dir, recursive = TRUE, showWarnings = FALSE)

# Input detection. The script first looks for common filenames, then falls back
# to keyword matching within data_dir. This makes the GitHub workflow portable
# while still failing clearly if files are missing.
find_input_file <- function(data_dir, preferred, keywords, label) {
  preferred_paths <- file.path(data_dir, preferred)
  existing <- preferred_paths[file.exists(preferred_paths)]
  if (length(existing) > 0) return(existing[[1]])
  
  all_xlsx <- list.files(
    data_dir,
    pattern = "\\.xlsx$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  matches <- all_xlsx[stringr::str_detect(tolower(basename(all_xlsx)), keywords)]
  
  if (length(matches) == 1) return(matches[[1]])
  
  if (length(matches) > 1) {
    stop(
      "More than one possible ", label, " file found in data_dir: ",
      paste(basename(matches), collapse = ", "),
      "\nPlease rename the target file or set the file path manually."
    )
  }
  
  stop(
    "No ", label, " .xlsx file found in: ", data_dir,
    "\nExpected one of: ", paste(preferred, collapse = ", "),
    "\nor a filename containing keyword pattern: ", keywords
  )
}

germ_file <- find_input_file(
  data_dir = data_dir,
  preferred = c("Combine.xlsx", "Germination.xlsx", "germination.xlsx"),
  keywords = "germ|combine",
  label = "germination"
)

endo_file <- find_input_file(
  data_dir = data_dir,
  preferred = c("Endophyte.xlsx", "endophyte.xlsx", "Endophyte_transmission.xlsx"),
  keywords = "endo|transmission",
  label = "endophyte"
)

message("Germination input: ", germ_file)
message("Endophyte input:    ", endo_file)

# ─────────────────────────────────────────────────────────────────────────────
# 2. Helper functions
# ─────────────────────────────────────────────────────────────────────────────

clean_names_simple <- function(x) {
  nm <- names(x)
  nm <- tolower(nm)
  nm <- gsub("[^a-z0-9]+", "_", nm)
  nm <- gsub("^_|_$", "", nm)
  names(x) <- nm
  x
}

normalise_seed_treatment <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, stringr::regex("control", ignore_case = TRUE)) ~ "Control",
    stringr::str_detect(x, stringr::regex("polymer|ploymer", ignore_case = TRUE)) ~ "Polymer",
    stringr::str_detect(x, stringr::regex("voc|6[- ]?pp|6pp", ignore_case = TRUE)) ~ "6-PP",
    TRUE ~ NA_character_
  )
}

normalise_seed_lot <- function(x) {
  dplyr::case_when(
    stringr::str_detect(x, stringr::regex("nil", ignore_case = TRUE)) ~ "Nil",
    stringr::str_detect(x, stringr::regex("ar37", ignore_case = TRUE)) ~ "AR37",
    TRUE ~ NA_character_
  )
}

fmt_mean_sd <- function(mean, sd, digits = 0) {
  paste0(round(mean, digits), " ± ", round(sd, digits))
}

extract_anova <- function(model, response_name) {
  as.data.frame(summary(model)[[1]]) |>
    tibble::rownames_to_column("Term") |>
    tibble::as_tibble() |>
    dplyr::mutate(Response = response_name, .before = 1)
}

extract_tukey <- function(model, response_name) {
  # Use emmeans::contrast() rather than emmeans::pairs().
  # pairs() is an S3 method for emmGrid objects, but it is not exported as
  # emmeans::pairs(), which can cause an error in reproducible scripts.
  as.data.frame(
    emmeans::contrast(
      emmeans::emmeans(model, ~ SeedTreatment),
      method = "pairwise",
      adjust = "tukey"
    )
  ) |>
    tibble::as_tibble() |>
    dplyr::mutate(Response = response_name, .before = 1)
}

extract_letters <- function(model) {
  multcomp::cld(
    emmeans::emmeans(model, ~ SeedTreatment),
    adjust = "tukey",
    Letters = letters
  ) |>
    as.data.frame() |>
    tibble::as_tibble() |>
    dplyr::mutate(Letters = stringr::str_trim(.group)) |>
    dplyr::select(SeedTreatment, Letters)
}

# ─────────────────────────────────────────────────────────────────────────────
# 3. Germination test
# ─────────────────────────────────────────────────────────────────────────────

germ_raw <- readxl::read_excel(germ_file) |> clean_names_simple()

required_germ <- c("treatment", "replicate", "germination")
missing_germ <- setdiff(required_germ, names(germ_raw))
if (length(missing_germ) > 0) {
  stop("Missing required germination column(s): ", paste(missing_germ, collapse = ", "))
}

germ <- germ_raw |>
  dplyr::mutate(
    Treatment_raw = as.character(treatment),
    SeedLot = normalise_seed_lot(Treatment_raw),
    SeedTreatment = normalise_seed_treatment(Treatment_raw),
    SeedLot = factor(SeedLot, levels = c("Nil", "AR37")),
    SeedTreatment = factor(SeedTreatment, levels = c("Control", "Polymer", "6-PP")),
    Replicate = as.integer(replicate),
    Germination_percent = as.numeric(germination)
  ) |>
  dplyr::select(
    Treatment_raw,
    SeedLot,
    SeedTreatment,
    Replicate,
    Germination_percent
  )

if (any(is.na(germ$SeedLot))) stop("Germination SeedLot contains NA. Check Treatment naming.")
if (any(is.na(germ$SeedTreatment))) stop("Germination SeedTreatment contains NA. Check Treatment naming.")
if (any(is.na(germ$Germination_percent))) stop("Germination_percent contains NA after conversion.")

germ_summary <- germ |>
  dplyr::group_by(SeedLot, SeedTreatment) |>
  dplyr::summarise(
    n = dplyr::n(),
    Germination_mean = mean(Germination_percent, na.rm = TRUE),
    Germination_sd = sd(Germination_percent, na.rm = TRUE),
    Germination_se = Germination_sd / sqrt(n),
    .groups = "drop"
  )

germ_nil <- dplyr::filter(germ, SeedLot == "Nil")
germ_ar37 <- dplyr::filter(germ, SeedLot == "AR37")

anova_germ_nil <- aov(Germination_percent ~ SeedTreatment, data = germ_nil)
anova_germ_ar37 <- aov(Germination_percent ~ SeedTreatment, data = germ_ar37)

letters_germ_nil <- extract_letters(anova_germ_nil) |> dplyr::mutate(SeedLot = "Nil", .before = 1)
letters_germ_ar37 <- extract_letters(anova_germ_ar37) |> dplyr::mutate(SeedLot = "AR37", .before = 1)
letters_germ <- dplyr::bind_rows(letters_germ_nil, letters_germ_ar37) |>
  dplyr::mutate(SeedLot = factor(SeedLot, levels = c("Nil", "AR37")))

germ_summary <- germ_summary |>
  dplyr::left_join(letters_germ, by = c("SeedLot", "SeedTreatment"))

# ─────────────────────────────────────────────────────────────────────────────
# 4. Endophyte transmission test
# ─────────────────────────────────────────────────────────────────────────────

endo_raw <- readxl::read_excel(endo_file) |> clean_names_simple()

required_endo <- c("treatment", "replicate", "endophyte")
missing_endo <- setdiff(required_endo, names(endo_raw))
if (length(missing_endo) > 0) {
  stop("Missing required endophyte column(s): ", paste(missing_endo, collapse = ", "))
}

endo <- endo_raw |>
  dplyr::mutate(
    Treatment_raw = as.character(treatment),
    SeedTreatment = normalise_seed_treatment(Treatment_raw),
    SeedTreatment = factor(SeedTreatment, levels = c("Control", "Polymer", "6-PP")),
    Replicate = as.integer(replicate),
    Endophyte_positive = as.numeric(endophyte),
    Transmission_percent = Endophyte_positive / 25 * 100,
    No_endophyte_percent = 100 - Transmission_percent
  ) |>
  dplyr::select(Treatment_raw, SeedTreatment, Replicate, Endophyte_positive,
         Transmission_percent, No_endophyte_percent)

if (any(is.na(endo$SeedTreatment))) stop("Endophyte SeedTreatment contains NA. Check Treatment naming.")
if (any(is.na(endo$Endophyte_positive))) stop("Endophyte_positive contains NA after conversion.")

endo_summary <- endo |>
  dplyr::group_by(SeedTreatment) |>
  dplyr::summarise(
    n = dplyr::n(),
    Transmission_mean = mean(Transmission_percent, na.rm = TRUE),
    Transmission_sd = sd(Transmission_percent, na.rm = TRUE),
    Transmission_se = Transmission_sd / sqrt(n),
    No_endophyte_mean = mean(No_endophyte_percent, na.rm = TRUE),
    No_endophyte_sd = sd(No_endophyte_percent, na.rm = TRUE),
    No_endophyte_se = No_endophyte_sd / sqrt(n),
    .groups = "drop"
  )

anova_endo <- aov(Transmission_percent ~ SeedTreatment, data = endo)
letters_endo <- extract_letters(anova_endo)

endo_summary <- endo_summary |>
  dplyr::left_join(letters_endo, by = "SeedTreatment")

# ─────────────────────────────────────────────────────────────────────────────
# 5. Manuscript Table 1
# ─────────────────────────────────────────────────────────────────────────────

# This creates a compact table for the manuscript. Germination is shown for both
# seed lots. Endophyte transmission is applicable only for AR37 seed lots.

germ_for_table <- germ_summary |>
  dplyr::mutate(
    Germination = fmt_mean_sd(Germination_mean, Germination_sd, digits = 0),
    Germination = paste0(Germination, " ", Letters)
  ) |>
  dplyr::select(SeedLot, SeedTreatment, Germination) |>
  tidyr::pivot_wider(
    names_from = SeedLot,
    values_from = Germination,
    names_prefix = "Germination_"
  )

endo_for_table <- endo_summary |>
  dplyr::mutate(
    Endophyte_transmission = fmt_mean_sd(Transmission_mean, Transmission_sd, digits = 0),
    Endophyte_transmission = paste0(Endophyte_transmission, " ", Letters),
    No_endophyte = fmt_mean_sd(No_endophyte_mean, No_endophyte_sd, digits = 0)
  ) |>
  dplyr::select(SeedTreatment, Endophyte_transmission, No_endophyte)

manuscript_table1 <- germ_for_table |>
  dplyr::left_join(endo_for_table, by = "SeedTreatment") |>
  dplyr::arrange(SeedTreatment) |>
  dplyr::transmute(
    `Seed treatment` = as.character(SeedTreatment),
    `Nil germination (%)` = Germination_Nil,
    `AR37 germination (%)` = Germination_AR37,
    `AR37 endophyte transmission (%)` = Endophyte_transmission,
    `AR37 seedlings without detectable endophyte (%)` = No_endophyte
  )

# ─────────────────────────────────────────────────────────────────────────────
# 6. Export outputs
# ─────────────────────────────────────────────────────────────────────────────

readr::write_csv(germ, file.path(results_dir, "germination_cleaned.csv"))
readr::write_csv(endo, file.path(results_dir, "endophyte_cleaned.csv"))
readr::write_csv(manuscript_table1, file.path(results_dir, "Table1_germination_endophyte_manuscript.csv"))
readr::write_csv(manuscript_table1, file.path(table1_dir, "Table1_germination_endophyte_manuscript.csv"))

writexl::write_xlsx(
  list(
    Manuscript_Table1 = manuscript_table1,
    Germination_Cleaned = germ,
    Germination_Summary = germ_summary,
    Germination_ANOVA_Nil = extract_anova(anova_germ_nil, "Germination_Nil"),
    Germination_Tukey_Nil = extract_tukey(anova_germ_nil, "Germination_Nil"),
    Germination_ANOVA_AR37 = extract_anova(anova_germ_ar37, "Germination_AR37"),
    Germination_Tukey_AR37 = extract_tukey(anova_germ_ar37, "Germination_AR37"),
    Endophyte_Cleaned = endo,
    Endophyte_Summary = endo_summary,
    Endophyte_ANOVA = extract_anova(anova_endo, "Endophyte_transmission"),
    Endophyte_Tukey = extract_tukey(anova_endo, "Endophyte_transmission")
  ),
  path = file.path(results_dir, "Table1_germination_endophyte_results.xlsx")
)

writexl::write_xlsx(
  list(Manuscript_Table1 = manuscript_table1),
  path = file.path(table1_dir, "Table1_germination_endophyte_manuscript.xlsx")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(results_dir, "sessionInfo_Table1_germination_endophyte.txt")
)

message("Done. Analysis results saved to: ", results_dir)
message("Done. Manuscript Table 1 saved to: ", table1_dir)
