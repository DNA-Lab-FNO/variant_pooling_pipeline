# Variant Pooling Pipeline (R)

Pipeline for loading, cleaning, and merging genomic variant data from Excel files.

**Outputs:**
- **Variants** – summary of unique variants with frequency, impact, and clinical significance  
- **SampleVariants** – variant-level data per sample  
- **Annotation** – detailed variant annotations  

**Usage:**  
1. Place Excel files in `input_path`  
2. Set `output_path` if needed  
3. Run `variant_pooling_pipeline.R` in R  

Generates a single `.xlsx` file with all three tables.

**Author:** Jana Schwarzerova | IWBBIO 2026  
**License:** MIT
