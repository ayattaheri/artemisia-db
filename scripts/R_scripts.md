## Table of Contents
  - [Table of Contents](#table-of-contents)
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
- [UMAP — Nearest Neighbors Parameter Tuning](#umap-nearest-neighbors-parameter-tuning)
- [Visualize All UMAP Results](#visualize-all-umap-results)
- [Save UMAP Comparison Panel](#save-umap-comparison-panel)
- [Select Optimal `n_neighbors` and Finalize UMAP Plot](#select-optimal-n_neighbors-and-finalize-umap-plot)
- [Save Final UMAP Plot and Coordinates](#save-final-umap-plot-and-coordinates)
- [Update `atlas_counts_sce` with Optimal UMAP](#update-atlas_counts_sce-with-optimal-umap)
- [Export Gene-Level Expression Data to Partitioned Parquet](#export-gene-level-expression-data-to-partitioned-parquet)
- [Export Gene- and Transcript-Level Expression Data by Plant Part](#export-gene--and-transcript-level-expression-data-by-plant-part)
- [Calculate Gene Expression Specificity Using the Tau (τ) Index](#calculate-gene-expression-specificity-using-the-tau-τ-index)
- [Load Expression Data and Compute Median TPM per Tissue](#load-expression-data-and-compute-median-tpm-per-tissue)
- [Filter for Expressed Genes](#filter-for-expressed-genes)
- [Calculate Tau Index for Expressed Genes](#calculate-tau-index-for-expressed-genes)
- [Classify Genes by Expression Pattern](#classify-genes-by-expression-pattern)
- [Step 30: Identify Tissues of Specific Expression for Tissue-Specific Genes](#step-30-identify-tissues-of-specific-expression-for-tissue-specific-genes)
- [Visualize Gene Classification Distribution](#visualize-gene-classification-distribution)
  - [Bar chart](#bar-chart)
  - [Pie chart](#pie-chart)
  - [Donut chart](#donut-chart)
- [UpSet Plot — Co-Specificity Across Tissues](#upset-plot-co-specificity-across-tissues)
- [Heatmap — Expression Profiles of Tissue-Specific Genes](#heatmap-expression-profiles-of-tissue-specific-genes)
- [Integrate Multi-Source Functional Annotations](#integrate-multi-source-functional-annotations)
- [Clean and Filter Annotation for the database](#clean-and-filter-annotation-for-the-database)
- [Prepare Functional Ter](#prepare-functional-ter)

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



# UMAP — Nearest Neighbors Parameter Tuning
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

# Prepare Functional Ter