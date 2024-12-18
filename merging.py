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
# Check for NaN values in the raw matrix
# print(np.isnan(adata.X.data).any())  # Check raw expression data

# # Check for NaN values in specific columns/rows of PCA
# print(np.isnan(adata.obsm['X_pca']).sum(axis=0))  # Per feature
# print(np.isnan(adata.obsm['X_pca']).sum(axis=1))  # Per cell




sc.pp.neighbors(adata, n_neighbors=10)  # Replace 'X_pca' with the desired representation

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
    "monocyte":["Csf1r","Csf1","Itgam","Ly6c1","Tlr4","Cd14","Ifngr1","Cx3cr1","Ccr2","Ccl2","Mafb","Spi1","Irf8","Batf3","Fcgr1","Mertk","Syk","Relb","Il1b","Nlrp3","Cd36","Hif1a","Stat3","Tnf"]
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