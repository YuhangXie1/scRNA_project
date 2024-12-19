import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
import scipy.sparse as sp
matplotlib.use("TkAgg")  # Or "Qt5Agg" depending on your setup
import matplotlib.pyplot as plt

sc.settings.set_figure_params(dpi=50, facecolor="white")

#loading in of files
file_paths = ["12w_KO_01_processed.h5ad",
            "12w_KO_02_processed.h5ad", 
            "12w_WT_01_processed.h5ad",
            "12w_WT_02_processed.h5ad",
            "ss_KO_processed.h5ad",             
            "ss_WT_processed.h5ad"
]

labels = ["12w_KO_01_processed.h5ad",
            "12w_KO_02_processed.h5ad", 
            "12w_WT_01_processed.h5ad",
            "12w_WT_02_processed.h5ad",
            "ss_KO_processed.h5ad",             
            "ss_WT_processed.h5ad"
]

adata_list = [sc.read_h5ad(f"mm_blood_10x/{file}") for file in file_paths]
adata = ad.concat(adata_list, axis=0, join="outer", label="source", keys=labels, index_unique="-")


# Check if the matrix is sparse
if sp.issparse(adata.X):
    # For sparse matrices, convert to dense format temporarily to check for NaN
    has_nan = np.isnan(adata.X.data).any()  # Check only the non-zero entries
else:
    # For dense matrices, use np.isnan directly
    has_nan = np.isnan(adata.X).any()

print(f"Does adata.X contain NaN values? {has_nan}")


sc.pp.neighbors(adata, n_neighbors=10)

#clustering
for res in [0.02, 0.35, 0.5, 2.0]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )

sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.35", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
    save = "cluster_different_leiden.png"
)

marker_genes = {
    "monocyte":["Itgam","Plac8","Csf1r","Ly6c2","Cx3cr1","Ace","Spn","Ccr2"],
    "T cells" : ["Cd3e","Cd3g","Cd3d","Ly6c2"],
    "NK cells": ["Nkg7", "Gzma"],
    "B cells" : ["Cd79a","Cd79b"],
    "Neutrophils": ["Ly6g","Csf3r"],
    "Densdritic": ["Clec10a","Clec9a","Xcr1","Flt3"]
}

sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var", save = "genes_leiden_res_0.02.png")
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.50", standard_scale="var", save = "genes_leiden_res_0.50.png")

"""
#differential gene expression
# Obtain cluster-specific differentially expressed genes
sc.tl.rank_genes_groups(adata, groupby="leiden_res_0.50", method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(
    adata, groupby="leiden_res_0.50", standard_scale="var", n_genes=5, save="cluster_genes_groups.png"
)
sc.get.rank_genes_groups_df(adata, group="7").head(5)

dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group="7").head(5)["names"]
sc.pl.umap(
    adata,
    color=[*dc_cluster_genes, "leiden_res_0.50"],
    legend_loc="on data",
    frameon=False,
    ncols=3,
    save="cluster_genes.png"
)
"""
#saving as data file
#adata.write("merged.h5ad")