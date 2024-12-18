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

adata = sc.read_h5ad("mm_blood_10x/monocytes_only.h5ad")

print(adata.obs.columns)

# Plot the UMAP colored by 'source' and Leiden clusters
sc.pl.umap(
    adata,
    color=["source"], 
    legend_loc="right margin",
    save="leiden_source_colored.png"
)

#differential gene expression
# Obtain cluster-specific differentially expressed genes
group = "leiden_res_0.02"
sc.tl.dendrogram(adata, groupby=group)

sc.tl.rank_genes_groups(adata, groupby=group, method="wilcoxon")
sc.pl.rank_genes_groups_dotplot(
    adata, groupby=group, standard_scale="var", n_genes=5, save=f"cluster_genes_groups_{group}.png"
)
sc.get.rank_genes_groups_df(adata, group="7").head(5)

dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group="7").head(5)["names"]
sc.pl.umap(
    adata,
    color=[*dc_cluster_genes, "leiden_res_0.50"],
    legend_loc="on data",
    frameon=False,
    ncols=3,
    save=f"cluster_genes_{group}.png"
)