import scanpy as sc
import anndata as ad
import pandas as pd
import numpy as np
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
# print(adata.shape)
# for i in range(0,5):
#     print(adata[i])

#QC
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("mt-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
# hemoglobin genes
adata.var["hb"] = adata.var_names.str.contains("^Hb[^(P)]")

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


# Step 2: Calculate the percentage of mitochondrial genes per cell
# Sum the counts for mitochondrial genes
mito_genes = adata.var_names.str.startswith('mt-')
adata.obs['mito_percentage'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1 * 100

# Step 3: Filter cells with high mitochondrial content
# You can adjust the threshold (e.g., 10% in this case)
mito_threshold = 10
adata = adata[adata.obs['mito_percentage'] < mito_threshold, :]

#removing cells with less than 200 genes. Removing genes with less than 3 cells
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

#doublet detection and removal
sc.pp.scrublet(adata)

#normalisation
# Saving count data
adata.layers["counts"] = adata.X.copy()
# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data
sc.pp.log1p(adata)

#feature selection
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pl.highly_variable_genes(adata)


#reduce dimensionality
sc.tl.pca(adata)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True)

# sc.pl.pca(
#     adata,
#     #color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
#     dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
#     ncols=2,
#     size=2,
# )

#nearest neighbour
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    #color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2,
)

#clustering
# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
sc.pl.umap(adata, color=["leiden"])

#refiltering
sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
)

sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
)

#saving as data file
#adata.write("ss_KO_processed.h5ad")

marker_genes = {
    "monocyte":["Csf1r","Csf1","Itgam","Ly6c1","Tlr4","Cd14","Ifngr1","Nos2","Cx3cr1","Ccr2","Ccl2","Mafb","Spi1","Irf8","Batf3","Fcgr1","Mertk","Syk","Relb","Il1b","Nlrp3","Cd36","Hif1a","Stat3","Tnf"]
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