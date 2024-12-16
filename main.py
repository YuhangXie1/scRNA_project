import scanpy as sc
import anndata as ad
import pooch
import pandas as pd
import matplotlib
matplotlib.use("TkAgg")  # Or "Qt5Agg" depending on your setup
import matplotlib.pyplot as plt

sc.settings.set_figure_params(dpi=50, facecolor="white")

#loading in of files

samples = {
    "12w_KO_01":"blood_12w_KO_01",
    "12w_KO_02":"blood_12w_KO_02",
    "12w_WT_01":"blood_12w_WT_01",
    "12w_WT_02":"blood_12w_WT_01",
    "ss_KO":"blood_ss_KO_B6",
    "ss_WT":"blood_ss_KO_B6"
}

adata = sc.read_10x_mtx("mm_blood_10x/mm_blood_10x/blood_ss_KO_B6/mtx/", cache=True)

#QC
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("MT-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^HB[^(P)]")

sc.pp.calculate_qc_metrics(
    adata,  qc_vars=["mt", "ribo", "hb"], inplace=True, log1p=True
)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
    jitter=0.4,
    multi_panel=True,
)

sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt")

#removing cells with less than 100 genes. Removing genes with less than 3 cells
sc.pp.filter_cells(adata, min_genes=100)
sc.pp.filter_genes(adata, min_cells=3)

#doublet detection and removal
#sc.pp.scrublet(adata, batch_key="sample")

#normalisation
# Saving count data
adata.layers["counts"] = adata.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

#feature selection
#sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
#sc.pl.highly_variable_genes(adata)


#reduce dimensionality
#sc.tl.pca(adata)
#sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)