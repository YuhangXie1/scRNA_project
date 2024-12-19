import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
import scipy.sparse as sp
matplotlib.use("TkAgg")  # Or "Qt5Agg" depending on your setup
import matplotlib.pyplot as plt
import time

sc.settings.verbosity = 1  # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=80, frameon=False, figsize=(3, 3), facecolor="white")

adata = sc.read_h5ad("mm_blood_10x/monocytes_only.h5ad")

adata = adata[adata.obs['leiden_res_0.50'] != '3'].copy() #selected all data but cluster 3
#sc.pl.umap(adata, color=["leiden_res_0.50"], legend_loc="on data")

#print(adata.obs.columns)


sc.tl.pca(adata)
sc.external.pp.bbknn(adata, batch_key="source")  # running bbknn 1.3.6
sc.tl.umap(adata)
sc.pl.umap(adata, color=["source"], save = "source_after_bbknn.png")
sc.pl.umap(adata, color=["leiden_res_0.50"], legend_loc = "on data", save = "leiden_res_0.50_after_bbknn.png")

# sc.pp.pca(adata)
# sc.pp.neighbors(adata)
# sc.tl.umap(adata)
# sc.pl.umap(
#     adata, color=["source"], palette=sc.pl.palettes.vega_20_scanpy, save = "recluster_source_after_bbknn.png"
# )
# sc.pl.umap(
#     adata, color=["leiden_res_0.50"], legend_loc = "on data", palette=sc.pl.palettes.vega_20_scanpy, save = "recluster_leiden0.5_after_bbknn.png"
# )

sc.tl.embedding_density(adata, groupby="source")
sc.pl.embedding_density(adata, groupby="source", save = "density_source_bbknn.png")


file_paths = ["12w_KO_01_processed.h5ad",
            "12w_KO_02_processed.h5ad", 
            "12w_WT_01_processed.h5ad",
            "12w_WT_02_processed.h5ad",
            "ss_KO_processed.h5ad",             
            "ss_WT_processed.h5ad"
]

for batch in file_paths:
    sc.pl.umap(adata, color="source", groups=[batch], save = f"by_source_{batch}.png")

for cluster in range(0,14,1):
    sc.pl.umap(adata, color="leiden_res_0.50", groups=[str(cluster)], save = f"by_leiden0.5_{cluster}.png")

#adata.write("integrated_data.h5ad")