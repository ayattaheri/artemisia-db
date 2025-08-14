# Artemisia Database Development: R Scripts
## Table of Contents
- [Load required packages](#load-required-packages)
- [Load and Filter Metadata by Mapping Rate](#load-and-filter-metadata-by-mapping-rate)
- [Merge Metadata with Alignment Results](#merge-metadata-with-alignment-results)
- [Map Samples to quant.sf Files](#map-samples-to-quantsf-files)
- [Prepare `tx2gene` Mapping from GTF Annotation](#prepare-tx2gene-mapping-from-gtf-annotation)
- [Validation of Quantification File Paths](#validation-of-quantification-file-paths)
- [Import Abundance Estimates Using `tximport`](#import-abundance-estimates-using-tximport)
- [Save Imported Data](#save-imported-data)
- [Transcript-Level Abundance Estimation for DTU Analysis](#transcript-level-abundance-estimation-for-dtu-analysis)
- [Create `SummarizedExperiment` Objects](#create-summarizedexperiment-objects)
- [Transcript-Level `SummarizedExperiment`](#transcript-level-summarizedexperiment)
- [Prepare SingleCellExperiment for Dimensionality Reduction](#prepare-singlecellexperiment-for-dimensionality-reduction)
- [Model Mean-Variance Relationship](#model-mean-variance-relationship)
- [Visualize the Mean-Variance Trend](#visualize-the-mean-variance-trend)
- [Save the plot as a PNG file](#save-the-plot-as-a-png-file)
- [Identify Highly Variable Genes (HVGs)](#identify-highly-variable-genes-hvgs)
- [Perform Principal Component Analysis (PCA)](#perform-principal-component-analysis-pca)
- [Visualize Variance Explained by Principal Components](#visualize-variance-explained-by-principal-components)
- [Save PCA Plots](#save-pca-plots)
- [Save Updated SCE Object](#save-updated-sce-object)
- [Prepare for Non-Linear Dimensionality Reduction](#prepare-for-non-linear-dimensionality-reduction)
- [t-SNE — Perplexity Parameter Tuning](#t-sne-perplexity-parameter-tuning)
- [Visualize All Perplexity Results](#visualize-all-perplexity-results)
- [Save Perplexity Comparison Panel](#save-perplexity-comparison-panel)
- [Select Optimal Perplexity and Finalize t-SNE Plot](#select-optimal-perplexity-and-finalize-t-sne-plot)
- [Save Final t-SNE Plot and Coordinates](#save-final-t-sne-plot-and-coordinates)
- [Update atlas_counts_sce with Optimal t-SNE](#update-atlas_counts_sce-with-optimal-t-sne)
- [# Step 20: UMAP — Nearest Neighbors Parameter Tuning](#step-20-umap-nearest-neighbors-parameter-tuning)
- [Visualize All UMAP Results](#visualize-all-umap-results)
- [Save UMAP Comparison Panel](#save-umap-comparison-panel)
- [Select Optimal `n_neighbors` and Finalize UMAP Plot](#select-optimal-n_neighbors-and-finalize-umap-plot)
- [Save Final UMAP Plot and Coordinates](#save-final-umap-plot-and-coordinates)
- [Update `atlas_counts_sce` with Optimal UMAP](#update-atlas_counts_sce-with-optimal-umap)
- [Export Gene-Level Expression Data to Partitioned Parquet](#export-gene-level-expression-data-to-partitioned-parquet)
    - [Output Schema](#output-schema)
- [Export Gene- and Transcript-Level Expression Data by Plant Part](#export-gene--and-transcript-level-expression-data-by-plant-part)
- [Calculate Gene Expression Specificity Using the Tau (τ) Index](#calculate-gene-expression-specificity-using-the-tau-τ-index)
- [Load Expression Data and Compute Median TPM per Tissue](#load-expression-data-and-compute-median-tpm-per-tissue)
- [Filter for Expressed Genes](#filter-for-expressed-genes)
- [Calculate Tau Index for Expressed Genes](#calculate-tau-index-for-expressed-genes)
- [Classify Genes by Expression Pattern](#classify-genes-by-expression-pattern)
- [Step 30: Identify Tissues of Specific Expression for Tissue-Specific Genes](#step-30-identify-tissues-of-specific-expression-for-tissue-specific-genes)
- [Visualize Gene Classification Distribution](#visualize-gene-classification-distribution)
  - [Bar chart](#bar-chart)
    - [Save Bar Chart](#save-bar-chart)
  - [Pie chart](#pie-chart)
    - [Save Pie Chart](#save-pie-chart)
  - [Donut chart](#donut-chart)
    - [Save Donut Chart](#save-donut-chart)
- [UpSet Plot — Co-Specificity Across Tissues](#upset-plot-co-specificity-across-tissues)
    - [Save UpSet Plot](#save-upset-plot)
- [Heatmap — Expression Profiles of Tissue-Specific Genes](#heatmap-expression-profiles-of-tissue-specific-genes)
    - [Save Heatmap](#save-heatmap)
- [Integrate Multi-Source Functional Annotations](#integrate-multi-source-functional-annotations)
- [Clean and Filter Annotation for the database](#clean-and-filter-annotation-for-the-database)
- [Prepare Functional Term Sets](#prepare-functional-term-sets)
  - [GO Terms (from InterPro)](#go-terms-from-interpro)
  - [InterPro Domains](#interpro-domains)
  - [KEGG Orthology (KO) Terms](#kegg-orthology-ko-terms)
- [Over-Representation Analysis (ORA) for Tissue-Specific Genes](#over-representation-analysis-ora-for-tissue-specific-genes)
- [Visualize Enriched Biological Processes Across Tissues](#visualize-enriched-biological-processes-across-tissues)
  - [Save Functional Heatmap](#save-functional-heatmap)
- [Identify Tissue-Specific Transcription Factors (TFs)](#identify-tissue-specific-transcription-factors-tfs)
  - [Load and Map PlantTFDB Annotations](#load-and-map-planttfdb-annotations)
  - [Count TFs by Family](#count-tfs-by-family)
  - [Identify Tissue-Specific TFs](#identify-tissue-specific-tfs)
  - [Heatmap of Tissue-Specific Transcription Factors](#heatmap-of-tissue-specific-transcription-factors)
  - [Bar Plot of All TF Families (PlantTFDB-Based)](#bar-plot-of-all-tf-families-planttfdb-based)
  - [Identify TFs Using PFAM Domains](#identify-tfs-using-pfam-domains)
  - [Export Interactive TF Family Table](#export-interactive-tf-family-table)
- [Identify Tissue-Specific TFs Using PFAM Domains](#identify-tissue-specific-tfs-using-pfam-domains)
  - [Create Heatmap of PFAM-Based Tissue-Specific TFs](#create-heatmap-of-pfam-based-tissue-specific-tfs)
- [Integrate TF Annotations from PlantTFDB and PFAM](#integrate-tf-annotations-from-planttfdb-and-pfam)
- [Create Final Heatmap — Integrated Tissue-Specific TFs](#create-final-heatmap-integrated-tissue-specific-tfs)
  - [Export Final Tissue-Specific TF Table](#export-final-tissue-specific-tf-table)
- [Identify CRISPR Arrays in the Transcriptome](#identify-crispr-arrays-in-the-transcriptome)
- [Epigenetic Regulators](#epigenetic-regulators)
  - [Identify m⁶A RNA Methylation Regulators](#identify-m⁶a-rna-methylation-regulators)
- [Identify DNA Methylation-Related Genes](#identify-dna-methylation-related-genes)
  - [Identify Histone H3 Modification Genes](#identify-histone-h3-modification-genes)
- [Artemisinin-Related Gene Analysis](#artemisinin-related-gene-analysis)
  - [Identify Artemisinin Biosynthesis-Related Genes](#identify-artemisinin-biosynthesis-related-genes)
  - [Median Expression of Artemisinin-Related Genes Across Tissues](#median-expression-of-artemisinin-related-genes-across-tissues)
  - [Tissue-Specificity of Artemisinin-Related Genes](#tissue-specificity-of-artemisinin-related-genes)
  - [Expression Classification of Artemisinin-Related Genes](#expression-classification-of-artemisinin-related-genes)
- [Final Wrap-Up](#final-wrap-up)
  - [Save Key Objects for Future Use](#save-key-objects-for-future-use)

---
title: "01_creating_the_AEA"
author: "Ayat Taheri"
date: "2024-05-17"
output: 
  html_document:
    toc: true
    toc_depth: 3
    fig_width: 8
    fig_height: 6
    theme: cosmo
editor_options: 
  chunk_output_type: console
---


```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)


# Set seed for reproducibility
set.seed(123)
```

# Load required packages
```{r include=FALSE}
library(here)
library(tidyverse)
library(clusterProfiler)
library(tximport)
library(tximportData)
library(SummarizedExperiment)
library(SingleCellExperiment)
library(scran)
library(scater)
library(arrow)
library(BioNERO)
library(ComplexHeatmap)
library(patchwork)
library(Biostrings)
library(planttfhunter)
library(grid)
library(GOfuncR)
library(RColorBrewer)
library(Rtsne)
library(cluster)
library(factoextra)
library(wordcloud)
library(knitr)
library(kableExtra)
```


# Load and Filter Metadata by Mapping Rate
We begin by loading the full metadata and the Salmon alignment results. Only samples with a mapping rate >50% are retained.

```{r include=FALSE}
# Load full metadata
all_metadata <- read_csv(here("data/metadata/all_metadata.csv")) |>
  select(Run, Project, condition, tissue, developmental_stage, genotype, bio_replicate)

# Load mapping rates from Salmon
mapping_results <- read_csv(here("data/07_salmon/mapping_rates.csv")) |>
  mutate(Run = gsub("_quant$", "", Run, perl = TRUE))  # Remove '_quant' suffix safely

# Filter: Keep only samples with >50% mapping rate
filtered_mapping_results <- mapping_results |>
  filter(mapping_rate > 50)
 
```

# Merge Metadata with Alignment Results
Join the filtered mapping results with the main metadata to create the final sample set.

```{r include=FALSE}
final_metadata <- inner_join(
  all_metadata, 
  select(filtered_mapping_results, Run, mapping_rate), 
  by = "Run"
) |>
  arrange(Project, Run)

# Save cleaned metadata for downstream use
write_csv(final_metadata, here("data/metadata/final_metadata.csv"))
```

# Map Samples to quant.sf Files
Build a data frame linking each sample (by Run) to its quant.sf file path.

```{r include=FALSE}
main_directory <- here("data/07_salmon/quants/")

# Get all project folders (immediate subdirectories)
project_dirs <- list.dirs(main_directory, recursive = FALSE, full.names = TRUE)

sample_info_list <- lapply(project_dirs, function(project_dir) {
  project_name <- basename(project_dir)
  sample_dirs <- list.dirs(project_dir, recursive = FALSE, full.names = TRUE)
  
  # For each sample directory, extract name and path to quant.sf
  lapply(sample_dirs, function(sample_dir) {
    sample_name <- gsub("_quant$", "", basename(sample_dir))
    quant_file <- file.path(sample_dir, "quant.sf")
    
    # Optional: warn if quant.sf doesn't exist
    if (!file.exists(quant_file)) {
      warning("Missing quant.sf for sample: ", sample_name)
    }
    
    data.frame(Sample = sample_name, FilePath = quant_file, stringsAsFactors = FALSE)
  }) |>
    do.call(rbind, args = _)
}) |>
  do.call(rbind, args = _)

# Convert to tibble for consistency
sample_info <- as_tibble(sample_info_list) |>
  filter(Sample %in% final_metadata$Run)  # Only keep filtered samples

# Reorder to match final_metadata
sample_info <- sample_info[match(final_metadata$Run, sample_info$Sample), ]

# Reset row names
rownames(sample_info) <- NULL

# Save mapping
write_csv(sample_info, here("data/07_salmon/sample_to_quant_path.csv"))
```


# Prepare `tx2gene` Mapping from GTF Annotation
To aggregate transcript-level abundance estimates (from Salmon) to the gene level, we need a mapping between transcript IDs and gene IDs — known as `tx2gene`.

We extract this from the GTF file generated during transcriptome assembly.

```{r create-tx2gene}
# Import GTF annotation
gtf_file <- here("data/05_assembly/mikado.loci.gtf")
gtf <- rtracklayer::import(gtf_file)

# Convert to data.frame for easier manipulation
gtf_df <- as.data.frame(gtf)

# Inspect feature types in the GTF
cat("Feature types in GTF:\n")
print(
  gtf_df %>%
    group_by(type) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(desc(count))
)

# Extract transcript-to-gene mapping for coding and non-coding RNAs
tx2gene <- gtf_df %>%
  filter(type %in% c("mRNA", "ncRNA")) %>%
  select(transcript_id, gene_id) %>%
  distinct() %>%
  rename(TXNAME = transcript_id, GENEID = gene_id)

# Validate required columns
if (!all(c("TXNAME", "GENEID") %in% colnames(tx2gene))) {
  stop("tx2gene must contain columns 'TXNAME' and 'GENEID'.")
}

# Save for reuse
write_tsv(tx2gene, here("data/05_assembly/tx2gene.txt"))
```

# Validation of Quantification File Paths
```{r}
# Load sample-to-file mapping (created in previous step)
sample_info <- read_csv(here("data/07_salmon/sample_to_quant_path.csv"))
files <- sample_info$FilePath

# Check existence
missing_files <- files[!file.exists(files)]
if (length(missing_files) > 0) {
  stop("Missing quant.sf files:\n", paste(missing_files, collapse = "\n"))
}

cat("All", length(files), "quant.sf files found.\n")
```

# Import Abundance Estimates Using `tximport`
We use the bias-corrected counts without an offset, specifically setting countsFromAbundance = "lengthScaledTPM".

Why `lengthScaledTPM`?
  - Converts transcript-level TPMs into gene-level counts.
  - Accounts for:
      - Library size (via scaling)
      - Transcript length differences
  - Produces scaled counts suitable for tools like DESeq2 and edgeR.
  - Avoids need for an offset matrix in downstream models.

```{r include=FALSE}
# Run tximport
txi <- tximport(
  files = files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "lengthScaledTPM"
)

# Add sample names (from Run IDs)
colnames(txi$abundance) <- sample_info$Sample
colnames(txi$counts)     <- sample_info$Sample
colnames(txi$length)     <- sample_info$Sample

# Confirm import
cat("tximport completed.\n")
cat("Number of genes detected:", nrow(txi$counts), "\n")
cat("Number of samples:", ncol(txi$counts), "\n")

# Show structure
names(txi)
head(txi$counts) %>% kable(caption = "Head of gene-level counts (lengthScaledTPM)")
```



# Save Imported Data

```{r}
output_path <- here("data/07_salmon/gene_level_abundances.rds")
saveRDS(txi, file = output_path)
```

# Transcript-Level Abundance Estimation for DTU Analysis
For **Differential Transcript Usage (DTU)** analysis, we need transcript-level abundance estimates that account for gene-wise structure. We use `tximport` with `countsFromAbundance = "dtuScaledTPM"`:

- Scales by **median transcript length per gene**
- Then adjusts for **library size**
- Preserves transcript-level structure (`txOut = TRUE`)
- Produces estimates suitable for tools like `DRIMSeq`, `IsoformSwitchAnalyzeR`, or `SUPPA2`

```{r import-transcript-level}
# Run tximport for transcript-level (DTU-ready) data
txi.tx <- tximport(
  files = files,
  type = "salmon",
  tx2gene = tx2gene,
  countsFromAbundance = "dtuScaledTPM",
  txOut = TRUE  # Keep transcript-level output
)

# Assign sample names
colnames(txi.tx$abundance) <- sample_info$Sample
colnames(txi.tx$counts)     <- sample_info$Sample
colnames(txi.tx$length)     <- sample_info$Sample


# Save for downstream use
saveRDS(txi.tx, here("data/07_salmon/transcript_level_abundances.rds"))
```

# Create `SummarizedExperiment` Objects
We now wrap both gene- and transcript-level data into `SummarizedExperiment` containers for standardized access and interoperability with Bioconductor tools.
```{r}
# Ensure colData is properly aligned
coldata <- final_metadata %>%
  dplyr::select(Run, Project, condition, tissue, developmental_stage, genotype, bio_replicate, mapping_rate) %>%
  dplyr::rename(sample = Run) %>%
  column_to_rownames("sample")

# Create gene-level SE
se_atlas_gene <- SummarizedExperiment(
  assays = list(
    gene_TPM = txi$abundance,
    gene_counts = txi$counts
  ),
  colData = coldata
)

# Round counts (they should already be numeric, but rounding ensures no floating point)
assay(se_atlas_gene, "gene_counts") <- round(assay(se_atlas_gene, "gene_counts"))

# Add row names (gene IDs) from txi
rownames(se_atlas_gene) <- rownames(txi$counts)

# Save
saveRDS(se_atlas_gene, here("data/08_dge/se_atlas_gene.rds"))
```

# Transcript-Level `SummarizedExperiment`
```{r}
se_atlas_transcript <- SummarizedExperiment(
  assays = list(
    tx_TPM = txi.tx$abundance,
    tx_counts = txi.tx$counts
  ),
  colData = coldata
)

# Add row names (transcript IDs)
rownames(se_atlas_transcript) <- rownames(txi.tx$counts)

# Save
saveRDS(se_atlas_transcript, here("data/08_dge/se_atlas_transcript.rds"))
```

# Prepare SingleCellExperiment for Dimensionality Reduction
Even for bulk RNA-seq, `SingleCellExperiment` (SCE) is useful for normalization, variance modeling, and visualization.

We use gene counts and log-normalized values.

```{r}
atlas_counts_sce <- SingleCellExperiment(
  assays = list(
    counts = assay(se_atlas_gene, "gene_counts"),
    logcounts = log2(assay(se_atlas_gene, "gene_counts") + 1)
  ),
  colData = colData(se_atlas_gene)
)

# Ensure row and column names are set
rownames(atlas_counts_sce) <- rownames(se_atlas_gene)
colnames(atlas_counts_sce) <- colnames(se_atlas_gene)

# Save SCE object
saveRDS(atlas_counts_sce, here("data/08_dge/atlas_counts_sce.rds"))

```

# Model Mean-Variance Relationship
We model the mean-variance trend to identify highly variable features (genes), even in bulk data.
```{r include=FALSE}
# Estimate biological and technical variance
mean_var_model <- modelGeneVar(atlas_counts_sce)

# Extract fit data
fit_data <- data.frame(
  mean = mean_var_model$mean,
  total = mean_var_model$total,
  bio = mean_var_model$bio,
  trend = mean_var_model$trend(mean_var_model$mean)
)

```

# Visualize the Mean-Variance Trend
```{r}
p_fit_mean_var <- ggplot(fit_data, aes(x = mean, y = total)) +
  geom_point(alpha = 0.07, color = "gray30") +
  geom_line(aes(y = trend), color = "steelblue3", linewidth = 1) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    title = "Per-Gene Mean-Variance Relationship",
    subtitle = "Log-normalized counts; trend line shows expected technical noise",
    x = "Mean of Log-Expression",
    y = "Variance of Log-Expression"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

# Display plot
print(p_fit_mean_var)

```
# Save the plot as a PNG file
```{r include=FALSE}
# Ensure directory exists
dir.create(here("data/figures"), showWarnings = FALSE, recursive = TRUE)

ggsave(
  here("data/figures/mean_var_scatterplot.png"),
  plot = p_fit_mean_var,
  width = 8, height = 6, dpi = 300
)
```


# Identify Highly Variable Genes (HVGs)
We extract the top 3000 genes with the highest **biological component of variation** using the mean-variance model. These genes capture biological heterogeneity (e.g., due to condition, tissue, or development) and are ideal for downstream dimensionality reduction.

```{r identify-hvgs}
# Extract top 3000 highly variable genes
hvg_genes <- getTopHVGs(mean_var_model, n = 3000)

```

#The object `hvg` is a character vector containing the IDs of the 
#top 3000 genes with the highest biological components.


# Perform Principal Component Analysis (PCA)
We compute PCA using only the top HVGs to focus on biologically meaningful variation.

```{r}
# Perform PCA on logcounts, using only HVGs
atlas_counts_sce <- fixedPCA(
  atlas_counts_sce,
  subset.row = hvg_genes,
  exprs_values = "logcounts",
  rank = 20  # Compute first 20 PCs
)
```

# Visualize Variance Explained by Principal Components

```{r}
# Extract percent variance explained
percent_var <- attr(reducedDim(atlas_counts_sce, "PCA"), "percentVar")

# Create data frame
pca_data <- data.frame(
  PC = factor(1:length(percent_var), levels = 1:length(percent_var)),
  Individual = round(percent_var, 2),
  Cumulative = round(cumsum(percent_var), 2)
)

# Plot 1: Individual variance per PC
p_pca_individual <- ggplot(pca_data, aes(x = PC, y = Individual)) +
  geom_col(fill = "grey70", width = 0.6) +
  geom_text(aes(label = Individual), vjust = -0.5, size = 3, color = "black") +
  labs(
    title = "Variance Explained by Each Principal Component",
    x = "Principal Component (PC)",
    y = "Variance Explained (%)"
  ) +
  coord_flip() +
  theme_minimal() +
  scale_y_continuous(limits = c(0, max(pca_data$Individual) * 1.1))

# Plot 2: Cumulative variance
p_pca_cumulative <- ggplot(pca_data, aes(x = as.numeric(PC), y = Cumulative)) +
  geom_line(color = "grey50", linewidth = 1) +
  geom_point(size = 2.5, color = "#69b3a2") +
  geom_hline(yintercept = 90, linetype = "dashed", color = "red", alpha = 0.7) +
  annotate("text", x = 5, y = 93, label = "90% threshold", color = "red", size = 3.5) +
  labs(
    title = "Cumulative Variance Explained by Principal Components",
    x = "Number of Principal Components",
    y = "Cumulative Variance Explained (%)"
  ) +
  theme_minimal()

# Display both plots side by side
gridExtra::grid.arrange(p_pca_individual, p_pca_cumulative, nrow = 1)

```

# Save PCA Plots
```{r}
# Ensure figures directory exists
dir.create(here("data/figures"), showWarnings = FALSE, recursive = TRUE)

# Save individual variance plot

ggsave(
  here("data/figures/p_pca_percent_var.png"),
  plot = p_pca_individual,
  width = 8, height = 5, dpi = 300
)

# Save cumulative variance plot
ggsave(
  here("data/figures/p_cumulative_var.png"),
  plot = p_pca_cumulative,
  width = 8, height = 5, dpi = 300
)

```

# Save Updated SCE Object
```{r}
saveRDS(atlas_counts_sce, here("data/08_dge/atlas_counts_sce.rds"))

```


# Prepare for Non-Linear Dimensionality Reduction

Based on the PCA results, the **top 12 principal components explain over 70% of the variance**, capturing major sources of biological variation.

We will use these 12 PCs as input for **t-SNE** (and later UMAP), which improves stability and reduces noise.

# t-SNE — Perplexity Parameter Tuning
t-SNE performance depends heavily on the perplexity parameter (roughly a smooth measure of effective number of neighbors). We test several values to find the best balance between:

Local structure preservation (low perplexity)
Global structure (high perplexity)
We evaluate perplexities: `5, 10, 15, 20, 25, 30`.


```{r}
perplexities <- c(5, 10, 15, 20, 25, 30)

# Generate t-SNE plots for each perplexity
p_tsne_list <- lapply(perplexities, function(perm) {
  # Run t-SNE using top 12 PCs
  tsne_result <- runTSNE(
    atlas_counts_sce,
    dimred = "PCA",
    n_dimred = 12,
    perplexity = perm,
    max_iter = 1000,
    theta = 0.5
  )
  
  # Plot, colored by 'Part' (or any batch-like factor)
  p <- plotReducedDim(tsne_result, dimred = "TSNE", colour_by = "Part") +
    labs(
      title = paste0("t-SNE | Perplexity = ", perm),
      x = "t-SNE 1", y = "t-SNE 2"
    ) +
    theme_minimal()
  
  return(p)
})

# Name plots for easy reference
names(p_tsne_list) <- paste0("perplexity_", perplexities)
```
# Visualize All Perplexity Results

```{r include=FALSE}
# Combine all plots
p_tsne_panel <- wrap_plots(p_tsne_list, ncol = 3) + 
  plot_layout(guides = "collect") & 
  scale_color_d3(palette = "category20") & 
  theme(
    legend.position = "right",
    plot.title = element_text(size = 10)
  ) +
  labs(color = "Part")

# Display
print(p_tsne_panel)

```

# Save Perplexity Comparison Panel
```{r include=FALSE}
# Ensure directory exists
dir.create(here("data/figures"), showWarnings = FALSE, recursive = TRUE)

ggsave(
  here("data/figures/p_tsne_all_perplexities.png"),
  plot = p_tsne_panel,
  width = 12, height = 8, dpi = 300
)


```

# Select Optimal Perplexity and Finalize t-SNE Plot
After visual inspection, `perplexity = 25` provides the best balance:

- Clear separation of sample groups
- No excessive fragmentation or over-clumping
```{r}
# Extract optimal plot (index 5 = perplexity 25)
p_tsne_optimal <- p_tsne_list[[5]]

# Improve title and styling
p_tsne_optimal <- p_tsne_optimal +
  labs(
    title = "t-SNE: Transcriptome-Wide Sample Relationships",
    subtitle = "Computed from top 12 PCs; Perplexity = 25",
    caption = "Colored by sequencing 'Part' (technical batch)"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# Display
print(p_tsne_optimal)

```

# Save Final t-SNE Plot and Coordinates
```{r include=FALSE}
# Save plot
ggsave(
  here("data/figures/p_tsne_optimal_perplexity.png"),
  plot = p_tsne_optimal,
  width = 8, height = 6, dpi = 300
)

# Extract coordinates for external use
tsne_coords <- reducedDim(atlas_counts_sce, "TSNE") %>%
  as.data.frame() %>%
  rownames_to_column("Run")

write_csv(tsne_coords, here("data/dimred/tsne_coordinates.csv"))
```
# Update atlas_counts_sce with Optimal t-SNE
```{r}
# Re-run t-SNE with optimal perplexity and attach to SCE
atlas_counts_sce <- runTSNE(
  atlas_counts_sce,
  dimred = "PCA",
  n_dimred = 12,
  perplexity = 25,
  name = "TSNE"
)

# Save updated SCE
saveRDS(atlas_counts_sce, here("data/08_dge/atlas_counts_sce.rds"))
```



# # Step 20: UMAP — Nearest Neighbors Parameter Tuning
We now apply **Uniform Manifold Approximation and Projection (UMAP)** using the top 12 principal components (PCs), as before.

UMAP’s `n_neighbors` parameter controls:
- **Local vs. global structure**: Low values emphasize local neighborhoods; high values preserve global topology.
- We test: `10, 20, 30, 40, 50, 60`


```{r}
# Define parameter values
n_neighbors <- c(10, 20, 30, 40, 50, 60)

# Generate UMAP plots for each n_neighbors
p_umap_list <- lapply(n_neighbors, function(nn) {
  # Run UMAP using top 12 PCs (corrected from 11 to 12 for consistency)
  umap_result <- runUMAP(
    atlas_counts_sce,
    dimred = "PCA",
    n_dimred = 12,         # Must match t-SNE input
    n_neighbors = nn,
    min_dist = 0.5,
    metric = "euclidean",
    max_iter = 500
  )
  
  # Plot, colored by 'Part' (technical batch)
  p <- plotReducedDim(umap_result, dimred = "UMAP", colour_by = "Part") +
    labs(
      title = paste0("UMAP | Nearest Neighbors = ", nn),
      x = "UMAP 1", y = "UMAP 2"
    ) +
    theme_minimal()
  
  return(p)
})

# Name plots for clarity
names(p_umap_list) <- paste0("nn_", n_neighbors)
```

# Visualize All UMAP Results
```{r include=FALSE}
# Combine all plots
p_umap_panel <- wrap_plots(p_umap_list, ncol = 3) +
  plot_layout(guides = "collect") &
  scale_color_d3(palette = "category20") &
  theme(
    legend.position = "right",
    plot.title = element_text(size = 10)
  ) +
  labs(color = "Part")

# Display
print(p_umap_panel)
```

# Save UMAP Comparison Panel
```{r include=FALSE}
# Ensure directory exists
dir.create(here("data/figures"), showWarnings = FALSE, recursive = TRUE)

ggsave(
  here("data/figures/p_umap_all_nneighbors_panel.png"),
  plot = p_umap_panel,
  width = 12, height = 8, dpi = 300
)
```

# Select Optimal `n_neighbors` and Finalize UMAP Plot

After visual inspection, `n_neighbors = 60` was selected because it:

- Preserves global structure of sample groups
- Avoids fragmentation
- Shows clean separation between known batches or conditions

```{r}
# Extract optimal plot (index 6 = n_neighbors = 60)
p_umap_optimal <- p_umap_list[[6]]

# Improve title and styling
p_umap_optimal <- p_umap_optimal +
  labs(
    title = "UMAP: Transcriptome-Wide Sample Relationships",
    subtitle = "Computed from top 12 PCs; n_neighbors = 60",
    caption = "Colored by sequencing 'Part' (technical batch)"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# Display
print(p_umap_optimal)
```

# Save Final UMAP Plot and Coordinates
```{r include=FALSE}
# Save plot
ggsave(
  here("data/figures/p_umap_optimal_nneighbors.png"),
  plot = p_umap_optimal,
  width = 8, height = 6, dpi = 300
)

# Extract coordinates
umap_coords <- reducedDim(atlas_counts_sce, "UMAP") %>%
  as.data.frame() %>%
  rownames_to_column("Run")

write_csv(umap_coords, here("data/dimred/umap_coordinates.csv"))
```

# Update `atlas_counts_sce` with Optimal UMAP
```{r}
# Re-run UMAP with optimal parameter and attach to SCE
atlas_counts_sce <- runUMAP(
  atlas_counts_sce,
  dimred = "PCA",
  n_dimred = 12,
  n_neighbors = 60,
  min_dist = 0.5,
  name = "UMAP"
)

# Save updated SCE
saveRDS(atlas_counts_sce, here("data/08_dge/atlas_counts_sce.rds"))
```
# Export Gene-Level Expression Data to Partitioned Parquet

To enable **efficient, scalable access** to the transcriptome atlas (e.g., for web apps, downstream analysis, or public archives), we export the gene-level expression matrix in **long format** as a **partitioned Apache Parquet dataset**.

### Output Schema

Each file contains a long-format data frame with:
- `Gene`: gene ID (character)
- `Sample`: sample name (factor)
- `TPM`: gene-level TPM (numeric)
- `Count`: bias-corrected counts (`lengthScaledTPM`) (numeric)
- `BioProject`: study/project ID (factor)
- `Part`: plant tissue/part (factor)

```{r include=FALSE}
# --- 1. Extract TPM and Count matrices in long format ---

# TPM
exp_tpm <- assay(se_atlas_gene, "gene_TPM") %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "TPM"
  ) %>%
  select(Gene, Sample, TPM)

# Counts
exp_counts <- assay(se_atlas_gene, "gene_counts") %>%
  as.data.frame() %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Count"
  ) %>%
  select(Gene, Sample, Count)

# --- 2. Validate alignment of gene-sample order ---
if (!identical(exp_tpm[, c("Gene", "Sample")], exp_counts[, c("Gene", "Sample")])) {
  stop("Gene and Sample order differs between TPM and Count matrices. Reorder before merging.")
}

# Combine
exp_final <- exp_tpm %>%
  mutate(Count = exp_counts$Count)

# --- 3. Add sample metadata (BioProject, Part) ---
sample_metadata <- colData(se_atlas_gene) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample")

sample_info <- sample_metadata %>%
  select(Sample, Project, Part) %>%
  rename(BioProject = Project)

# Join metadata
exp_final_with_meta <- exp_final %>%
  left_join(sample_info, by = "Sample") %>%
  mutate(
    Gene = as.character(Gene),
    Sample = as.factor(Sample),
    BioProject = as.factor(BioProject),
    Part = as.factor(Part),
    TPM = as.numeric(TPM),
    Count = as.numeric(Count)
  ) %>%
  # Reorder columns
  select(Gene, Sample, TPM, Count, BioProject, Part)

# --- 4. Write to partitioned Parquet dataset ---
parquet_output_path <- here("data/parquet_dir")

# Ensure directory is clean or removed first (optional)
if (dir_exists(parquet_output_path)) {
  dir_delete(parquet_output_path)
}
dir_create(parquet_output_path)

# Write dataset, partitioned by BioProject and Part
write_dataset(
  exp_final_with_meta,
  path = parquet_output_path,
  format = "parquet",
  partitioning = c("BioProject", "Part"),
  compression = "ZSTD",  # High compression, good speed
  overwrite = TRUE
)
```

# Export Gene- and Transcript-Level Expression Data by Plant Part
To support **interactive exploration**, **web applications**, or **focused tissue-specific analyses**, we export expression matrices **separately for each plant part (tissue)**.

Each tissue gets two files:
- `{Part}_TPM.tsv` — gene/transcript-level TPMs
- `{Part}_count.tsv` — bias-corrected counts

For **transcript-level** data:
- `{Part}_TPM_tx.tsv`
- `{Part}_count_tx.tsv`

All files are saved in **TSV format** for easy reading in R, Python, or Excel.
```{r}
# --- 1. Get Sample List by Plant Part (Tissue) ---
samples_by_part <- colData(se_atlas_gene) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample") %>%
  select(Sample, Part) %>%
  split(x = .$Sample, f = .$Part)

cat("Exporting expression data for", length(samples_by_part), "plant parts.\n")

# --- 2. Define Output Directory ---
outdir <- here("data", "app_data", "expression_by_body_part")
if (!dir.exists(outdir)) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
}

# --- 3. Helper Function: Export Matrix List ---
export_matrices <- function(
  assay_data,
  sample_groups,
  suffix = "",
  value_type = "TPM"
) {
  lapply(seq_along(sample_groups), function(i) {
    part <- names(sample_groups)[i]
    samples <- sample_groups[[i]]
    
    # Subset matrix and convert to data frame
    mat <- assay_data[, samples, drop = FALSE] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("Gene" %>% paste("ID", suffix, sep = "_"))
    
    # Build filename
    file <- file.path(outdir, paste0(part, "_", tolower(value_type), suffix, ".tsv"))
    
    # Write
    write_tsv(mat, file)
    cat("Saved:", file, "\n")
  })
}

# --- 4. Export Gene-Level Expression ---
tpm_gene <- assay(se_atlas_gene, "gene_TPM")
counts_gene <- assay(se_atlas_gene, "gene_counts")

export_matrices(tpm_gene, samples_by_part, suffix = "", value_type = "TPM")
export_matrices(counts_gene, samples_by_part, suffix = "", value_type = "count")

# --- 5. Export Transcript-Level Expression ---
tpm_tx <- assay(se_atlas_transcript, "tx_TPM")
counts_tx <- assay(se_atlas_transcript, "tx_counts")

# Update gene ID column name
export_matrices(tpm_tx, samples_by_part, suffix = "_tx", value_type = "TPM")
export_matrices(counts_tx, samples_by_part, suffix = "_tx", value_type = "count")


```

# Calculate Gene Expression Specificity Using the Tau (τ) Index
The **tau (τ) index** measures **tissue specificity** of gene expression:

- τ = 0: uniform expression across tissues (housekeeping-like)
- τ ≈ 1: highly tissue-specific expression

It is defined as:

$$
\tau = \frac{1}{n-1} \sum_{i=1}^{n} \left(1 - \frac{e_i}{\max(e)}\right)
$$

where $e_i$ is the expression in tissue $i$.

We compute τ using **median TPM per tissue**, focusing on **protein-coding tissues**.
```{r}
#' Calculate Tau (τ) Index for Tissue Specificity
#'
#' @param x A numeric vector of gene expression values (e.g., median TPM per tissue)
#' @return Tau index (0 = uniform, 1 = specific)
#'
#' @details Genes with all zero/low values return NA. Non-negative values required.
calculate_tau <- function(x) {
  if (all(is.na(x)) || all(x < 1, na.rm = TRUE)) {
    return(NA_real_)
  }
  if (min(x, na.rm = TRUE) < 0) {
    warning("Expression values < 0 detected. Skipping gene.")
    return(NA_real_)
  }
  max_val <- max(x, na.rm = TRUE)
  if (max_val == 0) return(0)
  
  scaled <- 1 - (x / max_val)
  sum(scaled, na.rm = TRUE) / (length(x) - 1)
}
```

# Load Expression Data and Compute Median TPM per Tissue

```{r}
# Define tissues of interest
target_parts <- c("Flower", "Leaf", "Petiole", "Root", "Stem")

# Open dataset
db <- open_dataset(here("data/parquet_dir"))

# Fetch all unique genes
all_genes <- db %>%
  select(Gene) %>%
  distinct() %>%
  collect() %>%
  pull(Gene)

# Split genes into chunks for memory-efficient processing
chunk <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))
gene_chunks <- chunk(all_genes, n = 100)

# Process each chunk: compute median TPM per gene per tissue
cat("Computing median TPM per tissue...\n")

median_per_part_list <- map_dfr(gene_chunks, ~ {
  message("Processing gene chunk ", which(gene_chunks == .x), "/", length(gene_chunks))
  
  db %>%
    select(Gene, Part, TPM) %>%
    filter(Part %in% target_parts, Gene %in% .x) %>%
    group_by(Gene, Part) %>%
    summarise(Median_TPM = median(TPM, na.rm = TRUE), .groups = "drop") %>%
    collect()
})

# Pivot to wide format: genes × tissues
median_per_part_wide <- median_per_part_list %>%
  pivot_wider(names_from = Part, values_from = Median_TPM, values_fill = 0) %>%
  na.omit() %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# Convert to data frame for export
median_per_part_df <- as.data.frame(median_per_part_wide)

# Save raw and processed data
write_csv(median_per_part_list, here("data/metadata/median_per_part_long.csv"))
write_csv(median_per_part_df, here("data/metadata/median_per_part.csv"))
```

# Filter for Expressed Genes
We define "expressed" genes as those with median TPM ≥ 1 in at least one tissue.

```{r}
# Identify genes with median TPM < 1 in *all* tissues
low_expression <- apply(median_per_part_wide, 1, function(x) all(x < 1))

# Keep only expressed genes
final_median_per_part <- median_per_part_wide[!low_expression, ]

# Save filtered matrix
final_median_per_part_df <- as.data.frame(final_median_per_part)
write_csv(final_median_per_part_df, here("data/metadata/final_median_per_part.csv"))

```
# Calculate Tau Index for Expressed Genes

```{r}
# Apply tau to log2-transformed values (optional, improves dynamic range)
tau_values <- apply(final_median_per_part, 1, function(x) calculate_tau(log2(x + 1)))

# Remove genes with NA tau
tau_values <- tau_values[!is.na(tau_values)]

# Save
tau_df <- data.frame(
  Gene = names(tau_values),
  Tau = as.numeric(tau_values)
) %>%
  arrange(desc(Tau))

write_csv(tau_df, here("data/metadata/genes_tau_index.csv"))
```
# Classify Genes by Expression Pattern

We classify genes into four categories:

1. Null expression: Not expressed in any tissue (median TPM < 1 everywhere)
2. Weak expression: Expressed, but median TPM < 5 in all tissues
3. Broadly expressed: τ < 0.85
4. Tissue-specific: τ ≥ 0.85

```{r include=FALSE}
# Classify expressed genes
gene_classes <- final_median_per_part_df %>%
  rownames_to_column("Gene") %>%
  as_tibble() %>%
  mutate(
    Max_Median_TPM = pmax(Flower, Leaf, Petiole, Root, Stem, na.rm = TRUE),
    Weak = Max_Median_TPM < 5,
    Tau = tau_df$Tau[match(Gene, tau_df$Gene)],
    Class = case_when(
      Tau < 0.85 ~ "Broadly expressed",
      Tau >= 0.85 ~ "Tissue-specific"
    )
  ) %>%
  mutate(
    Class = ifelse(Weak, "Weakly expressed", Class)
  ) %>%
  select(Gene, Class, Tau, Max_Median_TPM, everything()) %>%
  # Add null expression genes
  bind_rows(
    data.frame(Gene = rownames(median_per_part_wide)[low_expression]) %>%
      mutate(Class = "Null expression", Tau = NA_real_, Max_Median_TPM = NA_real_)
  )

# Summary
class_summary <- gene_classes %>%
  count(Class) %>%
  mutate(Percent = round(n / sum(n) * 100, 1))

kable(class_summary, caption = "Gene Expression Classification")

# Save classification
write_csv(gene_classes, here("data/metadata/gene_expression_classes.csv"))

```


# Step 30: Identify Tissues of Specific Expression for Tissue-Specific Genes
We now determine **in which plant parts** each **tissue-specific gene** is most highly expressed.

A gene is considered "specifically expressed" in a tissue if:
- It is classified as **tissue-specific (τ ≥ 0.85)**
- And has **median TPM > 5** in that tissue (biologically meaningful expression)

```{r}
# Ensure required object exists
if (!exists("gene_classes")) {
  gene_classes <- read_csv(here("data/metadata/gene_expression_classes.csv"))
}

# Use the long-format median expression data
median_per_part_long <- read_csv(here("data/metadata/median_per_part_long.csv"))

# Rename for clarity
classified_genes <- gene_classes %>%
  select(Gene, Class, Tau) %>%
  rename(Classification = Class)

# Add tissue-level median expression
classified_genes_long <- median_per_part_long %>%
  left_join(classified_genes, by = "Gene") %>%
  filter(Classification == "Tissue-specific") %>%
  mutate(Median = Median_TPM)  # For backward compatibility

# Identify the tissue(s) where each specific gene is expressed (Median > 5)
specific_genes_and_parts <- classified_genes_long %>%
  filter(Median > 5) %>%
  group_by(Gene) %>%
  summarise(
    Specific_parts = str_c(Part, collapse = ", "),
    .groups = "drop"
  )

# Handle genes with no tissue >5 TPM
all_specific_genes <- unique(classified_genes_long$Gene)
missing_parts <- setdiff(all_specific_genes, specific_genes_and_parts$Gene)
if (length(missing_parts) > 0) {
  missing_df <- data.frame(Gene = missing_parts, Specific_parts = "Low in all (≤5 TPM)")
  specific_genes_and_parts <- bind_rows(specific_genes_and_parts, missing_df)
}

# Final classification table
final_classified_genes <- classified_genes %>%
  mutate(Classification = fct_relevel(
    Classification,
    "Tissue-specific", "Broadly expressed", "Weakly expressed", "Null expression"
  )) %>%
  left_join(specific_genes_and_parts, by = "Gene") %>%
  arrange(Classification, Gene)

# Save
write_csv(final_classified_genes, here("data/metadata/final_classified_genes.csv"))
saveRDS(final_classified_genes, here("data/metadata/final_classified_genes.rds"))
```

# Visualize Gene Classification Distribution
## Bar chart
```{r}
# Summarize
class_summary <- final_classified_genes %>%
  count(Classification) %>%
  mutate(Classification = fct_relevel(
    Classification,
    "Tissue-specific", "Broadly expressed", "Weakly expressed", "Null expression"
  ))

# Define consistent color palette
class_colors <- c(
  "Tissue-specific"    = "#c60039",
  "Broadly expressed"  = "#ff8d1a",
  "Weakly expressed"   = "#57c785",
  "Null expression"    = "#2a7b9b"
)

p_bar <- class_summary %>%
  ggplot(aes(x = reorder(Classification, n), y = n, fill = Classification)) +
  geom_col() +
  geom_text(
    aes(label = scales::comma(n)),
    hjust = -0.1,
    size = 3.5,
    color = "black"
  ) +
  coord_flip() +
  labs(
    title = "Distribution of Genes Across Expression Categories",
    subtitle = "Genes classified by expression breadth (tau index) and level",
    x = "Expression Category",
    y = "Number of Genes"
  ) +
  scale_y_continuous(
    limits = c(0, max(class_summary$n) * 1.1),
    expand = expansion(mult = c(0, 0.1))
  ) +
  scale_fill_manual(values = class_colors) +
  theme_classic() +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold"),
    axis.title = element_text(size = 11)
  )

print(p_bar)
```
### Save Bar Chart
```{r}
dir.create(here("data/figures"), showWarnings = FALSE, recursive = TRUE)
ggsave(
  here("data/figures/p_genes_per_group_bar.png"),
  plot = p_bar,
  width = 8, height = 6, dpi = 300, bg = "white"
)

```



## Pie chart
```{r include=FALSE}
p_pie <- class_summary %>%
  ggplot(aes(x = "", y = n, fill = Classification)) +
  geom_bar(stat = "identity", width = 1, color = "white", linewidth = 0.3) +
  coord_polar("y", start = 0) +
  geom_text(
    aes(label = paste0(round(n / sum(n) * 100, 1), "%")),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "white",
    fontface = "bold"
  ) +
  labs(
    title = "Proportion of Genes by Expression Category",
    fill = "Category"
  ) +
  scale_fill_manual(values = class_colors) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold")
  )

print(p_pie)
```

### Save Pie Chart
```{r}
ggsave(
  here("data/figures/p_genes_per_group_pie.png"),
  plot = p_pie,
  width = 7, height = 7, dpi = 300, bg = "white"
)

```


## Donut chart
```{r include=FALSE}
# Summarize classification
class_summary <- final_classified_genes %>%
  mutate(Classification = fct_relevel(
    Classification,
    "Tissue-specific", "Broadly expressed", "Weakly expressed", "Null expression"
  )) %>%
  count(Classification) %>%
  mutate(
    percentage = n / sum(n) * 100,
    # Fix rounding error in last category
    rounded_percentage = ifelse(
      row_number() == n(),
      100 - sum(round(percentage[-n()], 1)),
      round(percentage, 1)
    )
  )

# Define hole size
hole_size <- 1.8

# Create donut chart
p_donut <- ggplot(class_summary, aes(x = hole_size, y = n, fill = Classification)) +
  geom_col(width = 1, color = "white", linewidth = 0.3) +
  coord_polar(theta = "y") +
  xlim(c(0.2, hole_size + 0.5)) +
  geom_text(
    aes(label = paste0(rounded_percentage, "%")),
    position = position_stack(vjust = 0.5),
    size = 3.5,
    color = "white",
    fontface = "bold"
  ) +
  theme_void() +
  labs(
    title = "Proportion of Genes by Expression Category",
    fill = "Category"
  ) +
  scale_fill_manual(values = c(
    "Tissue-specific"    = "#c60039",
    "Broadly expressed"  = "#ff8d1a",
    "Weakly expressed"   = "#57c785",
    "Null expression"    = "#2a7b9b"
  )) +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 10)
  )

print(p_donut)
```

### Save Donut Chart
```{r}
dir.create(here("data/figures"), showWarnings = FALSE, recursive = TRUE)

ggsave(
  here("data/figures/p_genes_per_group_donut.png"),
  plot = p_donut,
  width = 7, height = 7, dpi = 300, bg = "white"
)
```
# UpSet Plot — Co-Specificity Across Tissues
The UpSet plot shows how many genes are specific to one tissue vs. shared across multiple tissues.

This reveals:

- How many genes are truly single-tissue-specific
- How many are multi-tissue-specific
```{r}
# Extract tissue-specific genes with defined specific parts
specific_genes_long <- final_classified_genes %>%
  filter(Classification == "Tissue-specific", !is.na(Specific_parts)) %>%
  select(Gene, Specific_parts) %>%
  # Fix tissue capitalization
  mutate(Specific_parts = str_to_title(Specific_parts)) %>%
  # Split comma-separated tissues
  separate_longer_delim(Specific_parts, delim = ",")

# Trim whitespace
specific_genes_long$Specific_parts <- trimws(specific_genes_long$Specific_parts)

# Create list: tissue → genes
specific_genes_list <- split(
  specific_genes_long$Gene,
  specific_genes_long$Specific_parts
)

# Create combination matrix
comb_matrix <- make_comb_mat(specific_genes_list)
set_sizes <- set_size(comb_matrix)
comb_sizes <- comb_size(comb_matrix)
degrees <- comb_degree(comb_matrix)

# Custom color palette by degree (number of tissues)
custom_palette <- c("#d51438", "#05baab", "#f85433", "#345b81", "#57b87c")[degrees]

# Top annotation: Intersection sizes
top_annotation <- HeatmapAnnotation(
  "Intersection Size" = anno_barplot(
    comb_sizes,
    ylim = c(0, max(comb_sizes) * 1.1),
    border = FALSE,
    gp = gpar(fill = "#3b3b3b"),
    height = unit(5, "cm"),
    width = unit(2.1, "cm"),
    add_numbers = TRUE,
    numbers_rot = 0
  ),
  annotation_name_side = "left",
  annotation_name_rot = 90,
  annotation_name_gp = gpar(fontsize = 10)
)

# Left annotation: Genes per tissue
left_annotation <- rowAnnotation(
  "Genes per Tissue" = anno_barplot(
    -set_sizes,
    baseline = 0,
    axis_param = list(
      at = seq(0, -max(set_sizes), length.out = 6)[-1],
      labels = scales::comma(seq(max(set_sizes), 0, length.out = 6)[-6]),
      labels_rot = 0
    ),
    border = FALSE,
    gp = gpar(fill = "#ffa501"),
    width = unit(3.5, "cm")
  ),
  set_name = anno_text(
    set_name(comb_matrix),
    location = 0.2,
    just = "left",
    width = max_text_width(set_name(comb_matrix)) + unit(1.5, "mm"),
    gp = gpar(fontsize = 10)
  )
)

# Create UpSet plot
p_upset <- UpSet(
  comb_matrix,
  set_order = order(set_sizes),
  comb_order = order(degrees, -comb_sizes),
  comb_col = custom_palette,
  top_annotation = top_annotation,
  left_annotation = left_annotation,
  show_row_names = FALSE
)

# Display
draw(p_upset)

```

### Save UpSet Plot
```{r}
ggsave(
  here("data/figures/p_upset_1.png"),
  as.ggplot(p_upset),
  width = 8, height = 6, dpi = 300, bg = "white"
)
```

# Heatmap — Expression Profiles of Tissue-Specific Genes
We visualize the median expression (log2-TPM) of all tissue-specific genes across plant parts.

This reveals expression clusters and confirms specificity.

```{r include=FALSE}
# Extract tissue-specific genes
ts_genes <- final_classified_genes %>%
  filter(Classification == "Tissue-specific") %>%
  pull(Gene)

# Load median expression matrix
median_per_part_df <- read_csv(here("data/metadata/median_per_part.csv"), show_col_types = FALSE)

# Fix column names
colnames(median_per_part_df) <- str_to_title(colnames(median_per_part_df))

# Subset and log-transform
exp_matrix <- median_per_part_df %>%
  filter(Gene %in% ts_genes) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# Log2-transform: log2(TPM + 1)
exp_matrix_log <- log2(exp_matrix + 1)

# Define color scale
breaks <- c(
  seq(min(exp_matrix_log), 2, length.out = 50),
  seq(2, max(exp_matrix_log), length.out = 51)
)
color_palette <- rev(colorRampPalette(brewer.pal(9, "RdBu"))(100))

# Create heatmap
p_heatmap <- pheatmap(
  t(exp_matrix_log),
  main = "Median Expression of Tissue-Specific Genes",
  fontsize_row = 8,
  fontsize_col = 10,
  show_rownames = TRUE,
  show_colnames = FALSE,
  color = color_palette,
  breaks = breaks,
  legend = TRUE,
  border_color = NA
)
```
### Save Heatmap
```{r}
ggsave(
  here("data/figures/p_heatmap_median.png"),
  as.ggplot(p_heatmap),
  width = 8, height = 6, dpi = 300, bg = "white"
)
```


# Integrate Multi-Source Functional Annotations
We combine functional annotations from:
- **eggNOG-mapper**: Orthology, COG, KEGG KO
- **InterProScan**: Protein domains (InterPro), GO terms
- **BLASTx vs. Arabidopsis & UniProt**: Homology-based IDs
- **tx2gene mapping**: Link transcripts to genes

This creates a unified annotation table for downstream analysis.

```{r include=FALSE}
# --- 1. Load eggNOG-mapper annotations ---
eggnog <- read_tsv(
  here("data/08_functional_annotation/eggnog-mapper/hmmer.emapper.annotations"),
  skip = 4
) %>%
  rename(TXNAME = `#query`) %>%
  filter(!str_starts(TXNAME, "#"), row_number() <= n() - 3)

# --- 2. Load BLASTx vs. Arabidopsis ---
arabidopsis_blast <- read_tsv(
  here("data/08_functional_annotation/arabidopsis_blastx/transcriptome_vs_arabidopsis.outfmt6"),
  col_names = FALSE
) %>%
  rename(TXNAME = X1, Arabidopsis_ID = X2) %>%
  select(TXNAME, Arabidopsis_ID)

# --- 3. Load InterProScan ---
interpro <- read_tsv(
  here("data/08_functional_annotation/interproscan/mikado_interproscan.tsv"),
  col_names = FALSE
) %>%
  mutate(TXNAME = str_remove(X1, "_\\d+$")) %>%
  rename(
    Interpro_query = X1,
    InterPro_ID = X12,
    InterPro_Description = X13,
    GO_id = X14
  ) %>%
  select(TXNAME, InterPro_ID, InterPro_Description, GO_id)

# --- 4. Load BLASTx vs. UniProt ---
uniprot_blast <- read_tsv(
  here("data/08_functional_annotation/uniprot_blastx/transcriptome_vs_uniprot.outfmt6"),
  col_names = FALSE
) %>%
  rename(TXNAME = X1, Uniprot_full = X2) %>%
  mutate(Uniprot_ID = str_extract(Uniprot_full, "(?<=\\|)[^|]+(?=\\|)")) %>%
  select(TXNAME, Uniprot_ID)

# --- 5. Merge all sources ---
unified <- eggnog %>%
  select(TXNAME, eggNOG_OGs, COG_category, Description, KEGG_ko, Preferred_name) %>%
  full_join(interpro, by = "TXNAME") %>%
  full_join(arabidopsis_blast, by = "TXNAME") %>%
  full_join(uniprot_blast, by = "TXNAME") %>%
  full_join(tx2gene %>% rename(TXNAME = TXNAME), by = "TXNAME") %>%
  relocate(GENEID, .before = TXNAME) %>%
  arrange(GENEID)

# Save full annotation
write_csv(unified, here("data/08_functional_annotation/unified_annotation.csv"))

```
# Clean and Filter Annotation for the database

```{r}
# Extract relevant columns
unified_filtered <- unified %>%
  select(
    Gene = GENEID,
    InterPro_ID,
    InterPro_Description,
    eggNOG_OGs,
    COG_category,
    COG_Description = Description,
    Preferred_name,
    KEGG_ko,
    Arabidopsis_ID,
    Uniprot_ID
  )

# Save filtered version
write_csv(
  unified_filtered,
  here("data/08_functional_annotation/unified_annotation_filtered.csv"),
  na = ""
)

```

# Prepare Functional Term Sets
We prepare three types of functional terms for enrichment analysis.

## GO Terms (from InterPro)
```{r}
# Extract GO terms
go_terms <- unified %>%
  filter(!is.na(GO_id), GO_id != "-") %>%
  select(Gene = GENEID, GO_id) %>%
  mutate(GO_id = str_split(GO_id, "\\|")) %>%
  unnest(GO_id) %>%
  mutate(GO_id = str_remove(GO_id, "\\(InterPro\\)")) %>%
  filter(GO_id != "") %>%
  distinct()

# Get GO term descriptions
go_names <- get_names(go_terms$GO_id) %>%
  drop_na() %>%
  rename(GO_ID = go_id, GO_Description = go_name)

# Combine
unified_go <- left_join(
  go_terms %>% rename(GO_ID = GO_id),
  go_names,
  by = "GO_ID"
) %>%
  distinct()

# Save
write_csv(unified_go, here("data/08_functional_annotation/unified_go.csv"))

# Format for `clusterProfiler`
go_list <- list(
  term2gene = unified_go %>% select(GO = GO_ID, Gene = Gene),
  term2name = unified_go %>% select(GO = GO_ID, Description = GO_Description)
)

```
## InterPro Domains

```{r}
interpro_terms <- unified %>%
  filter(
    !is.na(InterPro_ID), InterPro_ID != "-",
    !is.na(InterPro_Description), InterPro_Description != "-"
  ) %>%
  select(Gene = GENEID, InterPro_ID, InterPro_Description)

interpro_list <- list(
  term2gene = interpro_terms %>% select(InterPro = InterPro_ID, Gene = Gene),
  term2name = interpro_terms %>% select(InterPro = InterPro_ID, Description = InterPro_Description)
)

```
## KEGG Orthology (KO) Terms
We extract KO IDs from eggNOG and map them to descriptions.

```{r}
# Extract KO IDs
kegg_ko <- unified %>%
  filter(!is.na(KEGG_ko), KEGG_ko != "-") %>%
  select(Gene = GENEID, KEGG_ko) %>%
  separate_longer_delim(KEGG_ko, delim = ",") %>%
  distinct()

# Save for external lookup (e.g., KEGG website)
write_tsv(
  kegg_ko$KEGG_ko %>% unique(),
  here("data/08_functional_annotation/kegg_ko_ids.txt"),
  col_names = FALSE
)

# Load manually retrieved KO descriptions
# (Obtained via KEGG website: https://www.genome.jp/kegg/ko.html)
ko_names <- read_tsv(
  here("data/08_functional_annotation/ko_names.txt"),
  col_names = c("KEGG_ko", "ko_Description")
) %>%
  distinct()

# Join
kegg_ko <- kegg_ko %>%
  left_join(ko_names, by = "KEGG_ko") %>%
  filter(!is.na(ko_Description))

# Save
write_csv(kegg_ko, here("data/08_functional_annotation/kegg_ko.csv"))

# Format for `clusterProfiler`
kegg_list <- list(
  term2gene = kegg_ko %>% select(KEGG = KEGG_ko, Gene = Gene),
  term2name = kegg_ko %>% select(KEGG = KEGG_ko, Description = ko_Description)
)
  
```
# Over-Representation Analysis (ORA) for Tissue-Specific Genes
We test for enriched functions in each tissue-specific gene set, using expressed genes as background.

```{r}
# Background: all expressed genes
background_genes <- final_classified_genes %>%
  filter(Classification != "Null expression") %>%
  pull(Gene) %>%
  unique()

# Perform ORA for each tissue
enrichment_results <- lapply(specific_genes_list, function(gene_set) {
  # GO enrichment
  go_enrich <- enricher(
    gene = gene_set,
    universe = background_genes,
    TERM2GENE = go_list$term2gene,
    TERM2NAME = go_list$term2name,
    pvalueCutoff = 1,
    qvalueCutoff = 0.05,
    minGSSize = 5
  )

  # InterPro enrichment
  interpro_enrich <- enricher(
    gene = gene_set,
    universe = background_genes,
    TERM2GENE = interpro_list$term2gene,
    TERM2NAME = interpro_list$term2name,
    pvalueCutoff = 1,
    qvalueCutoff = 0.05,
    minGSSize = 5
  )

  # KEGG enrichment
  kegg_enrich <- enricher(
    gene = gene_set,
    universe = background_genes,
    TERM2GENE = kegg_list$term2gene,
    TERM2NAME = kegg_list$term2name,
    pvalueCutoff = 1,
    qvalueCutoff = 0.05,
    minGSSize = 5
  )

  # Combine
  bind_rows(
    as.data.frame(go_enrich) %>% mutate(Category = "GO"),
    as.data.frame(interpro_enrich) %>% mutate(Category = "InterPro"),
    as.data.frame(kegg_enrich) %>% mutate(Category = "KEGG")
  )
})

# Combine all results
enrichment_df <- bind_rows(enrichment_results, .id = "Part") %>%
  filter(qvalue <= 0.05) %>%
  arrange(Part, qvalue)

# Save
write_tsv(
  enrichment_df,
  here("data/08_functional_annotation/enrichment_df.tsv")
)

```

# Visualize Enriched Biological Processes Across Tissues
We create a **presence/absence heatmap** of the **top 15 enriched terms** (GO, InterPro, KEGG) for each **tissue-specific gene set**.

This reveals:
- Which biological processes are uniquely or commonly enriched
- Functional specialization of each plant part
```{r}
# Extract top 15 most significant terms per tissue
top_terms <- enrichment_df %>%
  group_by(Part) %>%
  slice_min(qvalue, n = 15) %>%  # Use qvalue (FDR) instead of p.adjust
  ungroup() %>%
  mutate(Description = str_to_title(Description))  # Capitalize first letter of each word

# Create list: tissue → enriched terms
terms_list <- split(top_terms$Description, top_terms$Part)

# Convert to binary presence/absence matrix
pam <- list_to_matrix(terms_list)

# Define colors: absent = light blue, present = dark blue
custom_colors <- c("0" = "#ACC8E5", "1" = "#112A46")

# Create heatmap
p_heatmap_terms <- pheatmap(
  pam,
  color = custom_colors,
  breaks = c(0, 0.5, 1),
  name = "Term Enrichment",
  main = "Top Enriched Functions in Tissue-Specific Genes",
  fontsize_row = 6,
  fontsize_col = 8,
  cellwidth = unit(1.1, "cm"),
  border_color = "white",
  show_rownames = TRUE,
  show_colnames = TRUE,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)

# Adjust dendrogram size
p_heatmap_terms@row_dend_param$width <- unit(5, "mm")
p_heatmap_terms@column_dend_param$height <- unit(5, "mm")

```

## Save Functional Heatmap
```{r}
dir.create(here("data/figures"), showWarnings = FALSE, recursive = TRUE)

ggsave(
  here("data/figures/p_heatmap_terms_pav.png"),
  as.ggplot(p_heatmap_terms),
  width = 14, height = 10, dpi = 300, bg = "white"
)


```

# Identify Tissue-Specific Transcription Factors (TFs)

We identify transcription factors (TFs) within tissue-specific genes using annotations from PlantTFDB.

## Load and Map PlantTFDB Annotations
```{r}
# Load TF list from PlantTFDB
plant_tfs <- read_tsv(
  here("data/08_functional_annotation/transcription_factor/all_tf_list.txt"),
  col_names = TRUE
)

# Load BLAST results (transcriptome vs. PlantTFDB)
blast_results <- read_tsv(
  here("data/08_functional_annotation/transcription_factor/transcriptome_vs_all_tf.outfmt6"),
  col_names = FALSE
) %>%
  rename(TF_ID = X2, TXNAME = X1)

# Join to get TF annotations for transcripts
tf_annotations <- plant_tfs %>%
  select(TF_ID, Family) %>%
  inner_join(blast_results, by = "TF_ID") %>%
  select(TXNAME, TF_ID, Family)

# Map to genes using tx2gene
tf_gene_mapping <- tf_annotations %>%
  inner_join(tx2gene, by = "TXNAME") %>%
  select(Gene = GENEID, TF_ID, Family) %>%
  distinct()

# Save full TF list
write_csv(
  tf_gene_mapping,
  here("data/08_functional_annotation/transcription_factor/artemisia_all_tfdb.csv")
)
```
##  Count TFs by Family

```{r}
# Count TFs per family
tf_family_counts <- tf_gene_mapping %>%
  group_by(Family) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(desc(Count))

# Save
write_csv(
  tf_family_counts,
  here("data/08_functional_annotation/transcription_factor/artemisia_all_tfdb_count.csv")
)

# Optional: Bar plot of TF families
# ggplot(tf_family_counts, aes(x = reorder(Family, Count), y = Count)) +
#   geom_col(fill = "#57c785") + coord_flip() + theme_minimal()
```

## Identify Tissue-Specific TFs

```{r}
# Get tissue-specific genes
tissue_specific_genes <- final_classified_genes %>%
  filter(Classification == "Tissue-specific") %>%
  pull(Gene)

# Find TFs among them
tissue_specific_tfs <- tf_gene_mapping %>%
  filter(Gene %in% tissue_specific_genes)

# Add tissue context
tissue_specific_tfs_with_part <- tissue_specific_tfs %>%
  inner_join(
    final_classified_genes %>%
      filter(Gene %in% tissue_specific_tfs$Gene) %>%
      select(Gene, Specific_parts),
    by = "Gene"
  )

# Save
write_csv(
  tissue_specific_tfs_with_part,
  here("data/08_functional_annotation/transcription_factor/tissue_specific_tfs.csv")
)
```


## Heatmap of Tissue-Specific Transcription Factors
We visualize the **abundance of transcription factors (TFs)** among **tissue-specific genes**, showing which TF families dominate in each plant part.

```{r include=FALSE}
# Create TF count matrix: Family × Tissue
tf_counts_matrix <- tissue_specific_tfs_with_part %>%
  separate_rows(Specific_parts, sep = ", ") %>%
  mutate(Specific_parts = str_trim(Specific_parts)) %>%
  group_by(Specific_parts, Family) %>%
  summarise(count = n(), .groups = "drop") %>%
  pivot_wider(names_from = Family, values_from = count, fill = 0)

# Log2-transform for visualization
tf_log_matrix <- log2(as.matrix(tf_counts_matrix[, -1]) + 1)
rownames(tf_log_matrix) <- tf_counts_matrix$Specific_parts

# Define color palette
pal <- colorRampPalette(brewer.pal(9, "BuGn"))(100)[1:70]

# Create heatmap
p_heatmap_tfs <- pheatmap(
  tf_log_matrix,
  color = pal,
  display_numbers = as.matrix(tf_counts_matrix[, -1]),
  border_color = "gray90",
  name = "Log2 Count",
  main = "Transcription Factors in Tissue-Specific Genes",
  cellwidth = unit(0.6, "cm"),
  cellheight = unit(0.6, "cm"),
  fontsize_col = 7,
  fontsize_row = 8,
  angle_col = 45,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)

# Save the heatmap
dir.create(here("data/figures"), showWarnings = FALSE, recursive = TRUE)

ggsave(
  here("data/figures/p_heatmap_specific_tfs.png"),
  as.ggplot(p_heatmap_tfs),
  width = 10, height = 6, dpi = 300, bg = "white"
)

```

## Bar Plot of All TF Families (PlantTFDB-Based)

```{r}
# Load full TF count data
artemisia_all_tfdb_count <- read_csv(
  here("data/08_functional_annotation/transcription_factor/artemisia_all_tfdb_count.csv")
)

p_tf_bar <- ggplot(artemisia_all_tfdb_count, aes(x = reorder(Family, Count), y = Count)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Transcription Factor Family Distribution",
    subtitle = "Based on PlantTFDB annotations",
    x = "TF Family",
    y = "Number of TFs"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8))

print(p_tf_bar)

# Save the plot
ggsave(
  here("data/figures/planttfdb_all_bar.png"),
  p_tf_bar,
  width = 9, height = 7, dpi = 300, bg = "white"
)
```


## Identify TFs Using PFAM Domains
We cross-validate TF annotations using PFAM domain signatures from InterProScan.

```{r}
# Load PFAM-to-TF family mapping
pfam_map <- read_tsv(
  here("data/08_functional_annotation/transcription_factor/tf_family_pfam.ids"),
  col_names = c("Pfam_ID", "Family")
)

# Load InterProScan results
interpro_all <- read_tsv(
  here("data/08_functional_annotation/interproscan/mikado_interproscan.tsv"),
  col_names = FALSE
) %>%
  rename(
    Query = X1, Tool = X4, Pfam_ID = X5, Description = X6,
    Evalue = X9
  )

# Filter for high-confidence PFAM hits
pfam_hits <- interpro_all %>%
  filter(Tool == "Pfam", Evalue < 1e-5) %>%
  select(Query, Pfam_ID, Description)

# Extract TF-related hits
pfam_tfs <- pfam_hits %>%
  semi_join(pfam_map, by = "Pfam_ID") %>%
  mutate(TXNAME = str_remove(Query, "_\\d+$"))

# Map to genes
pfam_tfs_gene <- pfam_tfs %>%
  left_join(tx2gene %>% rename(TXNAME = TXNAME), by = "TXNAME") %>%
  select(Gene = GENEID, Pfam_ID, Family = Description) %>%
  left_join(pfam_map, by = "Pfam_ID") %>%
  select(Gene, Pfam_ID, Family = Family.y) %>%
  distinct()

# Save
write_csv(
  pfam_tfs_gene,
  here("data/08_functional_annotation/transcription_factor/pfam_tf_family.csv")
)

# Count by family
pfam_tf_counts <- pfam_tfs_gene %>%
  group_by(Family) %>%
  summarise(Count = n(), .groups = "drop") %>%
  arrange(desc(Count))

write_csv(
  pfam_tf_counts,
  here("data/08_functional_annotation/transcription_factor/pfam_tf_family_count.csv")
)

```

## Export Interactive TF Family Table

```{r}
# Use PlantTFDB-based counts
tf_family_table <- artemisia_all_tfdb_count %>%
  arrange(desc(Count)) %>%
  rename(`TF Family` = Family, Count = Count)

# Create styled HTML table
table_html <- kable(
  tf_family_table,
  format = "html",
  col.names = c("TF Family", "Count"),
  caption = "Transcription Factor Family Counts in the Transcriptome"
) %>%
  kable_styling(
    bootstrap_options = c("striped", "hover", "condensed", "responsive"),
    full_width = FALSE
  )

# Save
save_kable(table_html, file = here("docs/tf_family_count_table.html"))

```

# Identify Tissue-Specific TFs Using PFAM Domains
We identify **transcription factors (TFs)** based on **diagnostic PFAM domains**, independent of PlantTFDB, for cross-validation.
```{r}

# Load PFAM-to-TF mapping
pfam_map <- read_tsv(
  here("data/08_functional_annotation/transcription_factor/tf_family_pfam.ids"),
  col_names = c("Pfam_ID", "Family")
)

# Load filtered PFAM hits (from InterProScan, E-value < 1e-5)
filtered_pfam <- read_tsv(
  here("data/08_functional_annotation/interproscan/mikado_interproscan.tsv"),
  col_names = FALSE
) %>%
  rename(Query = X1, Tool = X4, Pfam_ID = X5, Evalue = X9) %>%
  filter(Tool == "Pfam", Evalue < 1e-5) %>%
  select(Query, Pfam_ID) %>%
  # Remove isoform suffix
  mutate(TXNAME = str_remove(Query, "_\\d+$"))

# Join with TF family map
pfam_tfs <- filtered_pfam %>%
  inner_join(pfam_map, by = "Pfam_ID") %>%
  inner_join(tx2gene %>% rename(TXNAME = TXNAME), by = "TXNAME") %>%
  select(Gene = GENEID, Family) %>%
  distinct()

# Save
write_csv(
  pfam_tfs,
  here("data/08_functional_annotation/transcription_factor/pfam_tf_family.csv")
)

```

## Create Heatmap of PFAM-Based Tissue-Specific TFs
Visualize TF family counts across tissues using PFAM-only annotations.

```{r include=FALSE}
# Get tissue-specific genes
tissue_specific_genes <- final_classified_genes %>%
  filter(Classification == "Tissue-specific") %>%
  select(Gene, Specific_parts)

# Combine with PFAM TFs
pfam_tf_enriched <- inner_join(tissue_specific_genes, pfam_tfs, by = "Gene") %>%
  mutate(Specific_parts = str_to_title(Specific_parts)) %>%
  separate_longer_delim(Specific_parts, delim = ",") %>%
  mutate(Specific_parts = str_trim(Specific_parts)) %>%
  group_by(Specific_parts) %>%
  count(Family, .drop = FALSE) %>%
  pivot_wider(names_from = Family, values_from = n, fill = 0) %>%
  column_to_rownames("Specific_parts") %>%
  as.matrix()

# Replace NA with 0
pfam_tf_enriched[is.na(pfam_tf_enriched)] <- 0

# Define color palette
pal <- colorRampPalette(brewer.pal(9, "BuGn"))(100)[1:70]

# Create heatmap
p_heatmap_pfam_tfs <- pheatmap(
  log2(pfam_tf_enriched + 1),
  color = pal,
  display_numbers = pfam_tf_enriched,
  border_color = "gray90",
  name = "Log2 Count",
  main = "PFAM-Based TFs in Tissue-Specific Genes",
  cellwidth = unit(0.6, "cm"),
  cellheight = unit(0.6, "cm"),
  fontsize_col = 7,
  fontsize_row = 8,
  angle_col = 45,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)

# Save PFAM TF Heatmap
ggsave(
  here("data/figures/p_heatmap_pfam_specific_tfs.png"),
  as.ggplot(p_heatmap_pfam_tfs),
  width = 12, height = 6, dpi = 300, bg = "white"
)

```

# Integrate TF Annotations from PlantTFDB and PFAM
We create a consensus TF annotation:

- Priority: PlantTFDB (curated database)
- Fallback: PFAM domain (functional signature)

This ensures maximum coverage with high confidence.

```{r}
# Load PlantTFDB-based TFs
plant_tfs <- read_csv(
  here("data/08_functional_annotation/transcription_factor/artemisia_all_tfdb.csv")
) %>%
  select(Gene, Family)

# Combine: use PlantTFDB first, then PFAM
combined_tfs <- full_join(plant_tfs, pfam_tfs, by = "Gene", suffix = c(".plant", ".pfam")) %>%
  mutate(Family = coalesce(Family.plant, Family.pfam)) %>%
  select(Gene, Family) %>%
  distinct()

# Save combined TF list
write_csv(
  combined_tfs,
  here("data/08_functional_annotation/transcription_factor/combined_tfs.csv")
)

```

# Create Final Heatmap — Integrated Tissue-Specific TFs
Generate a high-confidence heatmap using integrated TF annotations.

```{r include=FALSE}
# Get tissue-specific TFs
combined_tf_data <- inner_join(
  final_classified_genes %>% filter(Classification == "Tissue-specific"),
  combined_tfs,
  by = "Gene"
) %>%
  select(Gene, Specific_parts, Family) %>%
  mutate(Specific_parts = str_to_title(Specific_parts)) %>%
  separate_longer_delim(Specific_parts, delim = ",") %>%
  mutate(Specific_parts = str_trim(Specific_parts))

# Create count matrix
combined_tf_matrix <- combined_tf_data %>%
  group_by(Specific_parts, Family) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = Family, values_from = n, fill = 0) %>%
  column_to_rownames("Specific_parts") %>%
  as.matrix()

# Replace NA
combined_tf_matrix[is.na(combined_tf_matrix)] <- 0

# Save for web (Plotly)
saveRDS(
  combined_tf_matrix,
  here("data/08_functional_annotation/transcription_factor/combined_tf_counts_matrix.rds")
)

# Create heatmap
pal <- colorRampPalette(brewer.pal(9, "BuGn"))(100)[1:70]

p_heatmap_combined <- pheatmap(
  log2(combined_tf_matrix + 1),
  color = pal,
  display_numbers = combined_tf_matrix,
  border_color = "gray90",
  name = "Log2 Count",
  main = "Integrated TFs in Tissue-Specific Genes",
  cellwidth = unit(0.6, "cm"),
  cellheight = unit(0.6, "cm"),
  fontsize_col = 7,
  fontsize_row = 8,
  angle_col = 45,
  cluster_rows = TRUE,
  cluster_cols = TRUE
)

# Save Combined TF Heatmap

ggsave(
  here("data/figures/p_heatmap_combined_tfs.png"),
  as.ggplot(p_heatmap_combined),
  width = 12, height = 6, dpi = 300, bg = "white"
)

```

## Export Final Tissue-Specific TF Table

```{r}
# Add tissue context
combined_tf_counts <- combined_tf_data %>%
  arrange(Specific_parts, Family, Gene)

# Save
write_csv(
  combined_tf_counts,
  here("data/08_functional_annotation/transcription_factor/combined_tissue-specific-tfs.csv")
)
```


# Identify CRISPR Arrays in the Transcriptome
We analyze CRISPR loci identified by **CasFinder** to explore genome defense systems.

```{r}
# Load CRISPR predictions from CasFinder
raw_crispr <- read_tsv(
  here("data/08_functional_annotation/crispr/Crisprs_REPORT.tsv")
)

# Select and clean relevant columns
crispr_clean <- raw_crispr %>%
  select(
    Sequence = Sequence,
    CRISPR_Id,
    CRISPR_Start,
    CRISPR_End,
    CRISPR_Length,
    `Potential_Orientation (AT%)`,
    Consensus_Repeat
  ) %>%
  # Fix sequence naming (replace _ with . for compatibility)
  mutate(
    Sequence = sub("_(\\d+)$", ".\\1", Sequence),
    Sequence = gsub("mikado_Super", "mikado.Super", Sequence)
  )

# Map to genes using tx2gene
crispr_genes <- crispr_clean %>%
  left_join(tx2gene %>% rename(Sequence = TXNAME), by = "Sequence") %>%
  select(
    gene_id = GENEID,
    CRISPR_Start,
    CRISPR_End,
    CRISPR_Length,
    `Potential_Orientation (AT%)`,
    Consensus_Repeat
  ) %>%
  distinct()

# Save
write_csv(
  crispr_genes,
  here("data/08_functional_annotation/crispr/crispr_final.csv")
)
```

# Epigenetic Regulators
## Identify m⁶A RNA Methylation Regulators
We identify writers, readers, and erasers of N⁶-methyladenosine (m⁶A) — a key post-transcriptional regulator.

```{r}
# Define Pfam signatures for m⁶A regulators
m6a_families <- tibble(
  Type = c("Writer", "Writer", "Writer", "Writer", "Reader", "Eraser"),
  Gene_Family = c("MTA70", "WTAP", "HAKAI", "VIRILIZER", "YTH", "AlkB"),
  pfam_id = c("PF05063", "PF17098", "PF18408", "PF15912", "PF04146", "PF13532")
)

# Extract from InterProScan results
m6a_hits <- pfam_filtered_ids %>%
  rename(pfam_id = X5) %>%
  inner_join(m6a_families, by = "pfam_id") %>%
  select(TXNAME, Type, Gene_Family, pfam_id)

# Map to genes
m6a_genes <- m6a_hits %>%
  left_join(tx2gene %>% rename(TXNAME = TXNAME), by = "TXNAME") %>%
  select(gene_id = GENEID, Type, Gene_Family, pfam_id) %>%
  distinct()

# Save
write_csv(
  m6a_genes,
  here("data/08_functional_annotation/methylation/m6a/m6a_genes.csv")
)

# Summary
m6a_summary <- m6a_genes %>%
  group_by(Type) %>%
  summarise(Count = n(), .groups = "drop")

write_tsv(
  m6a_summary,
  here("data/08_functional_annotation/methylation/m6a/m6a_statistics.txt")
)  
```

# Identify DNA Methylation-Related Genes
We identify 5-methylcytosine (5mC) regulators, including DNA methyltransferases and Tet-like enzymes.

```{r}
# Define Pfam IDs for DNA methylation machinery
dna_methylation_families <- tibble(
  Short_name = c("DNA_methylase", "DNMT1-RFD", "TP_methylase"),
  Name = c(
    "C-5 cytosine-specific DNA methylase",
    "Cytosine specific DNA methyltransferase replication foci domain",
    "Tetrapyrrole (Corrin/Porphyrin) Methylases"
  ),
  Pfam_id = c("PF00145", "PF12047", "PF00590")
)

# Extract and map
five_methyl_genes <- pfam_filtered_ids %>%
  rename(Pfam_id = X5) %>%
  inner_join(dna_methylation_families, by = "Pfam_id") %>%
  select(TXNAME, Pfam_id, Short_name, Name) %>%
  left_join(tx2gene %>% rename(TXNAME = TXNAME), by = "TXNAME") %>%
  select(gene_id = GENEID, Pfam_id, Short_name, Name) %>%
  distinct()

# Save
write_csv(
  five_methyl_genes,
  here("data/08_functional_annotation/methylation/five_methyladenosine/five_methyl_genes.csv")
)
```

## Identify Histone H3 Modification Genes
We identify PHD and SET domain-containing proteins involved in histone H3 methylation.

```{r}
# Define Pfam IDs for histone modifiers
histone_families <- tibble(
  Short_name = c("PHD", "SET"),
  Name = c("PHD-finger", "SET domain"),
  Pfam_id = c("PF00628", "PF00856")
)

# Extract and map
histone_h3_genes <- pfam_filtered_ids %>%
  rename(Pfam_id = X5) %>%
  inner_join(histone_families, by = "Pfam_id") %>%
  select(TXNAME, Pfam_id, Short_name, Name) %>%
  left_join(tx2gene %>% rename(TXNAME = TXNAME), by = "TXNAME") %>%
  select(gene_id = GENEID, Pfam_id, Short_name, Name) %>%
  distinct()

# Save
write_csv(
  histone_h3_genes,
  here("data/08_functional_annotation/methylation/histone_h3/histone_h3_genes.csv")
)
```

# Artemisinin-Related Gene Analysis
## Identify Artemisinin Biosynthesis-Related Genes
We identify genes involved in **artemisinin biosynthesis** by integrating:
- Known biosynthetic genes from literature (e.g., *ADS*, *CYP71AV1*, *DBR2*, *ALDH1*)
- Homology-based annotation via BLAST
- UniProt and gene name mapping


```{r}
# Load reference list of artemisinin-related genes
art_ref <- read_csv(
  here("data/08_functional_annotation/artemisinin_related_genes/art_related_genes_hort.csv"),
  col_names = FALSE
) %>%
  rename(id = X2)

# Load BLAST results (transcriptome vs. reference sequences)
art_blast <- read_tsv(
  here("data/08_functional_annotation/artemisinin_related_genes/mikado_vs_artemisinin_hort.outfmt"),
  col_names = FALSE
) %>%
  rename(id = X1, TXNAME = X2)

# Join to find matches
art_matches <- inner_join(art_ref, art_blast, by = "id") %>%
  select(id, TXNAME)

# Clean UniProt IDs (extract primary ID between | |)
uniprot_clean <- uniprot %>%
  mutate(uniprot_id = str_extract(Uniprot_ID, "(?<=\\|)[^|]+(?=\\|)")) %>%
  rename(id = uniprot_id)

# Add UniProt annotation
art_with_annotation <- left_join(art_matches, uniprot_clean, by = "id") %>%
  mutate(TXNAME_final = coalesce(TXNAME.x, TXNAME.y)) %>%
  select(id, TXNAME = TXNAME_final) %>%
  na.omit()

# Map to genes
art_genes <- art_with_annotation %>%
  left_join(tx2gene %>% rename(TXNAME = TXNAME), by = "TXNAME") %>%
  select(GENEID, TXNAME, id) %>%
  rename(annotation = id) %>%
  distinct()

# Save
write_csv(
  art_genes,
  here("data/08_functional_annotation/artemisinin_related_genes/artemisinin_genes.csv")
)

```

## Median Expression of Artemisinin-Related Genes Across Tissues
We examine median TPM of artemisinin genes in each plant part.

```{r}

# Load median expression per gene per tissue
genes_median <- read_csv(
  here("data/metadata/median_per_part_long.csv"),
  col_types = cols(Median = col_double())
)

# Join with artemisinin genes
art_expression <- art_genes %>%
  rename(Gene = GENEID) %>%
  inner_join(genes_median, by = "Gene") %>%
  select(Gene, TXNAME, annotation, Part, Median_TPM = Median)

# Save
write_csv(
  art_expression,
  here("data/08_functional_annotation/artemisinin_related_genes/artemisinin_gene_expression.csv")
)
```


## Tissue-Specificity of Artemisinin-Related Genes
We determine which artemisinin genes are tissue-specific (τ ≥ 0.85).

```{r}
# Load gene classification
final_classified_genes <- read_csv(
  here("data/metadata/final_classified_genes.csv")
)

# Filter for tissue-specific artemisinin genes
specific_art_genes <- final_classified_genes %>%
  filter(Classification == "Tissue-specific") %>%
  inner_join(art_genes %>% rename(Gene = GENEID), by = "Gene")

# Save
write_csv(
  specific_art_genes,
  here("data/08_functional_annotation/artemisinin_related_genes/specific_artemisinin_genes.csv")
)
```

## Expression Classification of Artemisinin-Related Genes
We classify all artemisinin-related genes into expression categories.

```{r}
# Join all classifications
art_classified <- final_classified_genes %>%
  inner_join(art_genes %>% rename(Gene = GENEID), by = "Gene")

# Save
write_csv(
  art_classified,
  here("data/08_functional_annotation/artemisinin_related_genes/art_genes_category.csv")
)

# Summary statistics
art_summary <- art_classified %>%
  group_by(Classification) %>%
  summarise(Count = n(), .groups = "drop")

write_csv(
  art_summary,
  here("data/08_functional_annotation/artemisinin_related_genes/art_genes_category_stat.csv")
)

# Display
kable(art_summary, caption = "Expression Classification of Artemisinin-Related Genes")
```


# Final Wrap-Up
## Save Key Objects for Future Use
To enable **fast reloading** of results and plots without re-running the entire pipeline, we save essential objects in compressed `.rda` format.

```{r}
# Ensure directory exists
dir.create(here("data/figures"), showWarnings = FALSE, recursive = TRUE)

# Save key data objects
save(
  median_per_part,
  list = "median_per_part",
  file = here("data/figures/median_per_part.rda"),
  compress = "xz"
)

save(
  final_classified_genes,
  list = "final_classified_genes",
  file = here("data/figures/final_classified_genes.rda"),
  compress = "xz"
)

save(
  enrichment_df,
  list = "enrichment_df",
  file = here("data/figures/enrichment_df.rda"),
  compress = "xz"
)

# Save key plots (ggplotify::as.ggplot output from pheatmap/UpSet)
save(
  p_genes_per_group,
  list = "p_genes_per_group",
  file = here("data/figures/p_genes_per_group.rda"),
  compress = "xz"
)

save(
  p_upset,
  list = "p_upset",
  file = here("data/figures/p_upset.rda"),
  compress = "xz"
)

save(
  p_heatmap_median,
  list = "p_heatmap_median",
  file = here("data/figures/p_heatmap_median.rda"),
  compress = "xz"
)

save(
  p_heatmap_terms_pav,
  list = "p_heatmap_terms_pav",
  file = here("data/figures/p_heatmap_terms_pav.rda"),
  compress = "xz"
)

save(
  p_heatmap_specific_tfs,
  list = "p_heatmap_specific_tfs",
  file = here("data/figures/p_heatmap_specific_tfs.rda"),
  compress = "xz"
)
```
