import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("TkAgg")  # Or "Qt5Agg" depending on your setup
import matplotlib.pyplot as plt

sc.settings.set_figure_params(dpi=50, facecolor="white")

#loading in of files
adata = sc.read("mm_blood_10x/ss_KO_processed.h5ad")

#clustering
for res in [0.02, 0.5, 2.0]:
    sc.tl.leiden(
        adata, key_added=f"leiden_res_{res:4.2f}", resolution=res, flavor="igraph"
    )

sc.pl.umap(
    adata,
    color=["leiden_res_0.02", "leiden_res_0.50", "leiden_res_2.00"],
    legend_loc="on data",
)

marker_genes = {
    "monocyte":["Csf1r","Csf1","Itgam","Ly6c1","Tlr4","Cd14","Ifngr1","Cx3cr1","Ccr2","Ccl2","Mafb","Spi1","Irf8","Batf3","Fcgr1","Mertk","Syk","Relb","Il1b","Nlrp3","Cd36","Hif1a","Stat3","Tnf"]
}

missing_genes = [gene for gene in marker_genes["monocyte"] if gene not in adata.var_names]
print(f"Missing genes: {missing_genes}")


sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.02", standard_scale="var")
adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.02"].map(
    {
        "0": "Other cells",
        "1": "Monocytes",
    }
)
sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_0.50", standard_scale="var")