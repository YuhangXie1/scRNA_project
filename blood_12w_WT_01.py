# Code for scRNA-seq Hackathon project
# Quality control and clustering in sample "blood_12w_WT_01"
# 2024-12-17 Yueming Ren (and PyCharm AI assistant)

# Core scverse libraries
import scanpy as sc
import os

# Configure scanpy global parameters
sc.settings.set_figure_params(dpi=300, facecolor="white")

# Define the pull path to the main directory containing all sample folders
base_dir = "/Users/michaelren/Documents/2024_DPhil/RNAseq_DTC/Hackathon/mm_blood_10x" # Check the path!

# Use a function so I can process each dataset before integration
def scRNA_pp(
        adata,
        thresholds=None,
        plot_initial_qc=True,
        perform_scrublet=True,
        clustering_resolution=1.0,
):
    """
    Preprocess(pp) a single matrix from scRNA-seq results in AnnData
    Preprocess steps: QC, doublet detection, normalization, dimensionality reduction, and clustering.

    Parameters:
    -----------
    adata : AnnData
        The input AnnData object to preprocess.
    thresholds : dict, optional
        A dictionary of thresholds for filtering. Keys: 'pct_counts_mt', 'n_genes_by_counts', 'total_counts'.
        Example: {'pct_counts_mt': 10, 'n_genes_min': 200, 'n_genes_max': 6000, 'total_counts_min': 500}.
    plot_initial_qc : bool, optional
        If True, plots initial QC metrics for visual inspection.
    perform_scrublet : bool, optional
        If True, performs doublet detection using Scrublet.
    clustering_resolution : float, optional
        The resolution parameter for Leiden clustering.
    Plots will be saved in a subdirectory of the default directory /figures

    Returns:
    --------
    AnnData
        The processed AnnData object.
    """
    if thresholds is None:  # Default threshold
        thresholds = {
            'pct_counts_mt': 10,
            'n_genes_min': 200,
            'n_genes_max': 6000,
            'total_counts_min': 500
        }

    # Setup output directory for all plots
    output_dir = os.path.join("figures", dataset_name)  # Subdirectory in `figures/`
    os.makedirs(output_dir, exist_ok=True)  # Create directory if it doesn't exist
    sc.settings.figdir = output_dir  # Set Scanpy's global figure save location
    print(f"Saving plots in: {output_dir}")

    # Step 1: Calculate QC metrics
    adata.var["mt"] = adata.var_names.str.startswith("mt-")  # Identify mitochondrial genes (Mt- for mouse, MT- for human)
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, log1p=True)

    # Step 2: Save and Plot QC metrics before filtering
    if plot_initial_qc:
        print("Saving unfiltered QC metrics plots...")
        sc.pl.violin(
            adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True,
            save="_unfiltered_qc_violin.png"
        )
        sc.pl.scatter(
            adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
            save="_unfiltered_qc_scatterplot.png"
        )

    # Step 3: Apply filtering thresholds and filter the dataset
    print("Filtering cells based on thresholds...")
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    adata = adata[
        (adata.obs["pct_counts_mt"] < thresholds['pct_counts_mt']) &
        (adata.obs["n_genes_by_counts"] > thresholds['n_genes_min']) &
        (adata.obs["n_genes_by_counts"] < thresholds['n_genes_max']) &
        (adata.obs["total_counts"] > thresholds['total_counts_min'])
        ].copy()
    print(f"Number of cells remaining after filtering: {adata.n_obs}")

    # Step 4: Save and Plot QC metrics after filtering
    print("Saving filtered QC metrics plots...")
    sc.pl.violin(
        adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True,
        save="_filtered_qc_violin.png"  # Save in default directory './figures/'
    )
    sc.pl.scatter(
        adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt",
        save="_filtered_qc_total_vs_mt.png"  # Save with an informative filename
    )

    # Step 5: Doublet detection with Scrublet
    if perform_scrublet:
        print("Performing doublet detection using Scrublet...")
        sc.pp.scrublet(adata)  # Adds doublet scores and predictions to adata.obs

    # Step 6: Normalize and log-transform
    print("Normalizing data and performing log transformation...")
    # Save raw count data
    adata.layers["counts"] = adata.X.copy()  # Preserving raw count matrix in adata.layers["counts"]
    # Normalize to median total counts per cell
    sc.pp.normalize_total(adata)
    # Logarithmize the data
    sc.pp.log1p(adata)

    # Step 7: Identify highly variable genes
    print("Identifying highly variable genes...")
    sc.pp.highly_variable_genes(
        adata,
        # flavor="seurat_v3",
        n_top_genes=2000,
        # subset=True,  # Filter data object down to only HVGs
    )
    print("Saving HVG plot...")
    sc.pl.highly_variable_genes(adata, save="_hvg_plot.png")  # Save HVG plot in default figures directory

    # Step 8: PCA for dimensionality reduction
    print("Performing PCA...")
    sc.tl.pca(adata)
    print("Saving PCA results...")
    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, save="pca_variance_ratio.png")  # Save PCA elbow plot

    # Step 9: Compute neighborhood graph and UMAP embedding
    print("Computing UMAP embedding...")
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Step 10: Clustering (Leiden)
    print("Performing Leiden clustering on the neighborhood graph of cells...")
    sc.tl.leiden(adata, flavor="igraph", resolution=clustering_resolution, n_iterations=2)
    print("Saving Leiden clustering plots...")
    sc.pl.umap(adata, color=["leiden"], save="_leiden_clustering_umap.png")  # Save UMAP plot labeled by Leiden clusters

    return adata


# Choose dataset for testing
dataset_name = "blood_12w_WT_01"
dataset_path = os.path.join(base_dir, dataset_name, "mtx")

# Load the dataset
adata = sc.read_10x_mtx(dataset_path, var_names="gene_symbols", cache=True)
print(f"Dataset loaded: {adata.shape[0]} cells and {adata.shape[1]} genes")

# Run preprocessing with the scRNA_pp function
print(f"Preprocessing dataset: {dataset_name}")
processed_adata = scRNA_pp(
    adata,
    thresholds={
        'pct_counts_mt': 10,  # Threshold for mitochondrial content
        'n_genes_min': 200,  # Minimum genes per cell
        'n_genes_max': 6000,  # Maximum genes per cell
        'total_counts_min': 500  # Minimum total counts per cell
    },
    plot_initial_qc=True,  # Plot QC metrics
    perform_scrublet=True,  # Perform doublet detection
    clustering_resolution=1.0,  # Clustering resolution
)
