import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Set up scanpy settings for figure saving
sc.settings.figdir = "results/figures/"  # Directory where figures will be saved
sc.settings.set_figure_params(dpi=80, facecolor='white')

# Create figures directory if it doesn't exist
import os
os.makedirs("results/figures/", exist_ok=True)

print("Loading data...")

# Read the CORRECTED results (use filtered data for initial analysis)
# Modified to handle uncompressed files
adata = sc.read_10x_mtx(
    "results2/STARsolo/sample1_swapped_Solo.out/Gene/raw/",
    var_names="gene_symbols",
    make_unique=True,
    cache=False,  # Don't cache compressed files
    gex_only=False  # Handle different file formats
)

# Make variable names unique and set up
adata.var_names_make_unique()

adata.obs_names.is_unique


print(f"Initial data: {adata.n_obs} cells, {adata.n_vars} genes")

# Calculate QC metrics
print("Calculating QC metrics...")
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], inplace=True)

print("Basic statistics:")
print(adata.obs[['total_counts', 'n_genes_by_counts', 'pct_counts_mt']].describe())

# QC violin plots - SAVED
print("Generating QC violin plots...")
sc.pl.violin(
    adata, 
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter=0.4, 
    multi_panel=True,
    save="_qc_metrics.pdf"
)

# QC scatter plots - SAVED
print("Generating QC scatter plots...")
sc.pl.scatter(
    adata, 
    x='total_counts', 
    y='n_genes_by_counts',
    save="_counts_vs_genes.pdf"
)

sc.pl.scatter(
    adata, 
    x='total_counts', 
    y='pct_counts_mt',
    save="_counts_vs_mt.pdf"
)

# Show current cell numbers
print(f"Before filtering: {adata.n_obs} cells")

# Filter cells (adjust thresholds based on your QC plots)
print("Filtering cells and genes...")
sc.pp.filter_cells(adata, min_genes=200)  # Filter cells with too few genes
adata = adata[adata.obs.n_genes_by_counts < 5000, :]  # Filter cells with too many genes  
adata = adata[adata.obs.pct_counts_mt < 20, :]  # Filter high mitochondrial cells

# Filter genes
sc.pp.filter_genes(adata, min_cells=3)  # Filter genes expressed in few cells

print(f"After filtering: {adata.n_obs} cells, {adata.n_vars} genes")

# Post-filtering QC plots - SAVED
print("Generating post-filtering QC plots...")
sc.pl.violin(
    adata,
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    jitter=0.4,
    multi_panel=True,
    save="_qc_metrics_filtered.pdf"
)

# Save raw counts
adata.raw = adata

# Normalize to 10,000 reads per cell
print("Normalizing data...")
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Find highly variable genes
print("Finding highly variable genes...")
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

# Plot highly variable genes - SAVED
sc.pl.highly_variable_genes(adata, save="_highly_variable_genes.pdf")

print(f"Found {sum(adata.var.highly_variable)} highly variable genes")

# Keep only highly variable genes for downstream analysis
adata = adata[:, adata.var.highly_variable]

# Scale data
print("Scaling data...")
sc.pp.scale(adata, max_value=10)

# PCA
print("Running PCA...")
sc.tl.pca(adata, svd_solver='arpack')

# PCA variance plot - SAVED
sc.pl.pca_variance_ratio(
    adata, 
    log=True, 
    n_pcs=50,
    save="_pca_variance.pdf"
)


# PCA plot - SAVED  
sc.pl.pca(adata, save="_pca.pdf")

# Compute neighborhood graph
print("Computing neighborhood graph...")
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)

# UMAP
print("Running UMAP...")
sc.tl.umap(adata)

# UMAP plot - SAVED
sc.pl.umap(adata, save="_umap_basic.pdf")

# Leiden clustering
print("Performing Leiden clustering...")
sc.tl.leiden(adata, resolution=0.5)

# UMAP with clusters - SAVED
sc.pl.umap(
    adata, 
    color=['leiden'],
    legend_loc='on data',
    save="_umap_leiden_clusters.pdf"
)

# UMAP with QC metrics - SAVED
sc.pl.umap(
    adata,
    color=['total_counts', 'n_genes_by_counts', 'pct_counts_mt'],
    save="_umap_qc_metrics.pdf"
)

# Cluster statistics
print("\nCluster statistics:")
cluster_stats = adata.obs['leiden'].value_counts().sort_index()
print(cluster_stats)

# Save the processed data
print("Saving processed data...")
adata.write("results/processed_data.h5ad")

print("Analysis complete! All figures saved to results/figures/")
print("Processed data saved to results/processed_data.h5ad")

# Summary statistics
print(f"\nFinal results:")
print(f"Cells after filtering: {adata.n_obs}")
print(f"Genes after filtering: {adata.n_vars}")
print(f"Number of clusters: {len(adata.obs['leiden'].unique())}")
print(f"Median UMIs per cell: {np.median(adata.obs['total_counts'])}")
print(f"Median genes per cell: {np.median(adata.obs['n_genes_by_counts'])}")
