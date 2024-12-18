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

adata = sc.read_h5ad("mm_blood_10x/merged.h5ad")

#clustering - pass 1
for res in [0.02]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )

adata = adata[adata.obs["leiden_res_0.02"].isin(["0", "2"]), :]
# Recompute PCA (optional, but recommended after subsetting)
sc.pp.pca(adata, n_comps=50)

# Recompute the neighbors graph
sc.pp.neighbors(adata, n_neighbors=10, use_rep="X_pca")  # Ensure correct representation is used

#clustering - pass 2
for res in [0.02]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )
    
adata = adata[adata.obs["leiden_res_0.02"].isin(["0","1", "2"]), :]
# Recompute PCA (optional, but recommended after subsetting)
sc.pp.pca(adata, n_comps=50)

# Recompute the neighbors graph
sc.pp.neighbors(adata, n_neighbors=10, use_rep="X_pca")  # Ensure correct representation is used


#clustering
for res in [0.02, 0.35, 0.5]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )

sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.35", "leiden_res_0.50"],
    legend_loc="on data",
    save = "cluster_different_leiden.png"
)

marker_genes = {
    "monocyte":["Itgam","Plac8","Csf1r","Ly6c2","Cx3cr1","Ace","Spn","Ccr2"],
    "T cells" : ["Cd3e","Cd3g","Cd3e","Ly6c2"],
    "NK cells": ["Nkg7", "Gzma"],
    "B cells" : ["Cd79a","Cd79b"],
    "Neutrophils": ["Ly6g","Csf3r"],
    "Densdritic": ["Clec10a","Clec9a","Xcr1","Flt3"]
}

missing_genes = [gene for gene in marker_genes["monocyte"] if gene not in adata.var_names]
print(f"Missing genes: {missing_genes}")


sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var", save = "genes_leiden_res_0.02.png")
adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.02"].map(
    {
        "0": "Other cells",
        "1": "Monocytes",
    }

)
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.50", standard_scale="var", save = "genes_leiden_res_0.50.png")
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.35", standard_scale="var", save = "genes_leiden_res_0.35.png")
#adata.write("monocyte_only.h5ad")

"""
#differential gene expression
# Obtain cluster-specific differentially expressed genes
group = "leiden_res_0.02"
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
    save="cluster_genes.png"
)
"""