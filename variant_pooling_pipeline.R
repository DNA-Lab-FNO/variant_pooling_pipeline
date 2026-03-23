#===============================================================================
# variant_pooling_pipeline.R v4.2
# author: Jana Schwarzerova
# description: Variant-Pooling Pipeline for IWBBIO 2026
#===============================================================================

# Clean environment
rm(list = ls())

# ==============================
# Set working directory
# ==============================
setwd("D:/IWBBIO2026/IWBBIO_2026/Pooling")

# ==============================
# Load libraries
# ==============================
library(readxl)
library(dplyr)
library(purrr)
library(stringr)
library(writexl)

# ==============================
# Input / Output paths
# ==============================
input_path <- "Input_Rare_IWBBIO2026/"
output_path <- "output/variant_pooling_output.xlsx"

# ==============================
# 1. List Excel files
# ==============================
files <- list.files(
    path = input_path, 
    pattern = "\\.xlsx$", 
    full.names = TRUE
)

# Remove temp Excel files (~$...)
files <- files[!grepl("^~\\$", basename(files))]

stopifnot(length(files) > 0)

# ==============================
# 2. Load data safely
# ==============================
data_list <- files %>%
    set_names() %>%
    map(
        ~ tryCatch(
            read_excel(.x),
            error = function(e) {
                message(paste("Error reading:", .x))
                return(NULL)
            }
        )
    ) %>%
    compact() %>%
    imap(
        ~ .x %>%
            mutate(SAMPLE_ID = tools::file_path_sans_ext(basename(.y)))
    )

# ==============================
# 3. Merge all data
# ==============================
df <- bind_rows(data_list)

# ==============================
# 4. Data cleaning + normalization
# ==============================
df <- df %>%
    mutate(
        Chr = as.character(Chr),
        Coordinate = as.numeric(Coordinate),
        Reference = toupper(Reference),
        Alternate = toupper(Alternate),
        `Variant Frequency` = as.numeric(
            gsub(",", ".", gsub("%", "", as.character(`Variant Frequency`)))
        ),
        VARIANT_ID = str_c(Chr, Coordinate, Reference, Alternate, sep = "_")
    )

# ==============================
# 5. Filter low-confidence variants
# ==============================
df <- df %>%
    filter(!is.na(`Variant Frequency`) & `Variant Frequency` > 5)

# ==============================
# 6. Count number of samples
# ==============================
n_samples <- df %>%
    distinct(SAMPLE_ID) %>%
    nrow()

# ==============================
# 7. Create variant table
# ==============================
df <- df %>%
    mutate(
        var_occured = as.numeric(str_split(`Variant occurred`, "/", simplify = TRUE)[, 1]),
        total_alleles = as.numeric(str_split(`Variant occurred`, "/", simplify = TRUE)[, 2]),
        allele_fraction = var_occured / total_alleles
    )

variant_table <- df %>%
    group_by(VARIANT_ID) %>%
    summarise(
        Chr = first(Chr),
        Position = first(Coordinate),
        REF = first(Reference),
        ALT = first(Alternate),
        Gene = first(Symbol),
        Impact = case_when(
            any(Impact == "HIGH") ~ "HIGH",
            any(Impact == "MODERATE") ~ "MODERATE",
            any(Impact == "LOW") ~ "LOW",
            TRUE ~ "MODIFIER"
        ),
        Clinical_significance = case_when(
            any(str_detect(`Clinical significance`, "Pathogenic")) ~ "Pathogenic",
            any(str_detect(`Clinical significance`, "Likely pathogenic")) ~ "Likely pathogenic",
            any(str_detect(`Clinical significance`, "Uncertain")) ~ "VUS",
            any(str_detect(`Clinical significance`, "Likely benign")) ~ "Likely benign",
            any(str_detect(`Clinical significance`, "Benign")) ~ "Benign",
            TRUE ~ "Other"
        ),
        occurrence_count = n_distinct(SAMPLE_ID),
        occurrence_freq = occurrence_count / n_samples,
        mean_AF = mean(`Variant Frequency`, na.rm = TRUE),
        max_AF = max(`Variant Frequency`, na.rm = TRUE),
        mean_allele_fraction = mean(allele_fraction, na.rm = TRUE),
        max_allele_fraction = max(allele_fraction, na.rm = TRUE),
        .groups = "drop"
    )

# ==============================
# 8. Create sample-variant table
# ==============================
sample_variant_table <- df %>%
    select(
        SAMPLE_ID,
        VARIANT_ID,
        `Variant Frequency`,
        `Read Depth`,
        Genotype
    )

# ==============================
# 9. Create annotation table
# ==============================
annotation_table <- df %>%
    select(
        Gene = Symbol,
        Chr,
        Genotype,
        Variant_class,
        Coordinate,
        HGVSc,
        HGVSp,
        dbSNP = `VEP dbSNP ID`,
        Consequence,
        ClinVar = `Clinical significance`,
        ReadDepth = `Read Depth`,
        VARIANT_ID
    ) %>%
    distinct()

# ==============================
# 10. Export results
# ==============================
dir.create("output", showWarnings = FALSE)

write_xlsx(
    list(
        Variants = variant_table,
        SampleVariants = sample_variant_table,
        Annotation = annotation_table
    ),
    path = output_path
)

# ==============================
# DONE
# ==============================
cat("Pipeline completed successfully!\n")
cat("Number of samples:", n_samples, "\n")
cat("Number of unique variants:", nrow(variant_table), "\n")