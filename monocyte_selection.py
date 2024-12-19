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

marker_genes = {
    "monocyte":["Itgam","Plac8","Csf1r","Ly6c2","Cx3cr1","Ace","Spn","Ccr2"],
    "T cells" : ["Cd3e","Cd3g","Cd3d","Ly6c2"],
    "NK cells": ["Nkg7", "Gzma"],
    "B cells" : ["Cd79a","Cd79b"],
    "Neutrophils": ["Ly6g","Csf3r"],
    "Densdritic": ["Clec10a","Clec9a","Xcr1","Flt3"]
}

adata = sc.read_h5ad("mm_blood_10x/merged.h5ad")

def generate_cluster(adata, pass_key):
    sc.pp.pca(adata, n_comps=50)
    sc.pp.neighbors(adata, n_neighbors=10, use_rep="X_pca")

    for res in [0.02, 0.35, 0.50]:
        sc.tl.leiden(
            adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
        )

    #generating cluster map
    sc.pl.umap(
        adata,
        color=["leiden_res_0.02", "leiden_res_0.35", "leiden_res_0.50"],
        legend_loc="on data",
        save = f"cluster_different_leiden_pass_{pass_key}.png"
    )

    #generating gene plots
    sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var", save = f"genes_leiden_res_0.02_pass_{pass_key}.png")
    sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.35", standard_scale="var", save = f"genes_leiden_res_0.35_pass_{pass_key}.png")
    sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.50", standard_scale="var", save = f"genes_leiden_res_0.50_pass_{pass_key}.png")


#pass 1
generate_cluster(adata, pass_key = 1)

#pass 2
adata = adata[adata.obs["leiden_res_0.02"].isin(["0", "2"]), :]
generate_cluster(adata, pass_key = 2)

#pass 3
adata = adata[adata.obs["leiden_res_0.02"].isin(["0","1", "2"]), :]
generate_cluster(adata, pass_key = 3)




adata.write("monocyte_only.h5ad")

