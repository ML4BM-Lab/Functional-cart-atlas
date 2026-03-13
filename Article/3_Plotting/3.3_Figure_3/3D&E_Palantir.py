###############################################################################
###############################################################################

# Program: 3D&E_Palantir.py
# Author: Sergio Cámara Peña
# Date: 16/01/2026
# Version: FINAL

###############################################################################
###############################################################################

# %% Import libraries
import palantir
import scanpy as sc
import pandas as pd
import os
import numpy as np
import random
import matplotlib
import matplotlib.pyplot as plt
import scipy.stats as st

## Warnings
import warnings
from numba.core.errors import NumbaDeprecationWarning

warnings.filterwarnings(action="ignore", category=NumbaDeprecationWarning)
warnings.filterwarnings(
    action="ignore", module="scanpy", message="No data for colormapping"
)

## Custom functions
def get_expr_vector(adata, gene):
    """Devuelve expresion 1D (cells,) para un gen en ad.X."""
    if gene not in adata.var_names:
        return None
    gidx = np.where(adata.var_names == gene)[0][0]
    v = adata.X[:, gidx]
    if hasattr(v, "toarray"):  # sparse
        v = v.toarray().ravel()
    else:
        v = np.array(v).ravel()
    return v

def binned_means_zscore(adata, genes, pseudotime, n_bins=25):
    # quantile-based bins --> bins with a similar number of cells
    bins = pd.qcut(pseudotime, q=n_bins, duplicates="drop")
    codes = np.asarray(bins.codes)
    n_eff = codes.max() + 1

    M = np.zeros((len(genes), n_eff), dtype=float)
    for i, g in enumerate(genes):
        v = get_expr_vector(adata, g)
        for b in range(n_eff):
            idx = (codes == b)
            M[i, b] = np.mean(v[idx]) if np.any(idx) else np.nan

    # z-score per gen (row)
    mu = np.nanmean(M, axis=1, keepdims=True)
    sd = np.nanstd(M, axis=1, keepdims=True) + 1e-9
    Mz = (M - mu) / sd

    # mean pseudotime value per bin (to label the X-axis)
    bin_mids = np.array([np.mean(pseudotime[codes == b]) for b in range(n_eff)])
    return Mz, bin_mids

# %% Set random seed
random.seed(2504)

# %% Load data
os.chdir("/home/scamara/data_a/scamara/Atlas/Input")
ad_orig = sc.read("Python_scVI_adata_big_V4_state4_Normalized.h5ad")

# %% Set PATH to save figs
os.chdir("/home/scamara/data_a/scamara/Atlas/Figuras/Figura_3")

# %% Filter object
# Antigen == "Blood"
adata_filtered = ad_orig[ad_orig.obs['Antigen'] == "Blood"].copy()
print(adata_filtered.shape[0])

# manual_celltype_annotation_high == "CD8 cytotoxic"
adata_filtered = adata_filtered[adata_filtered.obs['manual_celltype_annotation_high'] == "CD8 cytotoxic"].copy()
print(adata_filtered.shape[0])

# STATUS is not "Healthy"
adata_filtered = adata_filtered[adata_filtered.obs['STATUS'] == "DISEASE"].copy()
print(adata_filtered.shape[0])

# Remove rows where Max_Response is NA
adata_filtered = adata_filtered[~adata_filtered.obs['Max_Response'].isna()].copy()

# Remove cells where Time_Point_Ranges == "Infusion_Product" and Stimulated == "YES"
mask = ~((adata_filtered.obs['Time_Point_Ranges'] == "Infusion_Product") & (adata_filtered.obs['Stimulated'] == "YES"))
adata_filtered = adata_filtered[mask].copy()
print(adata_filtered.shape[0])

# Remove cells stimulated in vitro
mask = ~((adata_filtered.obs['Stimulated'] == "YES") & (adata_filtered.obs['Stimulation_Location'] == "In_vitro"))
adata_filtered = adata_filtered[mask].copy()
print(adata_filtered.shape[0])

# Max_Response == "CR"
adata_filtered = adata_filtered[adata_filtered.obs['Max_Response'] == "CR"].copy()
print(adata_filtered.shape[0])

ad = adata_filtered.copy()

# %% Run diffusion maps
dm_res = palantir.utils.run_diffusion_maps(ad, n_components=5, pca_key="X_scVI")
ms_data = palantir.utils.determine_multiscale_space(ad)

# %% View UMAP representation
sc.pl.embedding(
    ad,
    basis="umap",
    frameon=False,
)

# %% Calculate start cell
X = ad.obsm["X_scVI"]
ip_idx = np.where(ad.obs["Time_Point_Ranges"].values == "Infusion_Product")[0]

ip_center = X[ip_idx].mean(axis=0)
start_cell = ad.obs_names[ip_idx[np.argmin(np.linalg.norm(X[ip_idx] - ip_center, axis=1))]]

print("start_cell:", start_cell)

# Show start cell on UMAP
palantir.plot.highlight_cells_on_umap(ad, [start_cell])
plt.show()

# %% Run Palantir
pr_res = palantir.core.run_palantir(
    ad,
    early_cell=start_cell,
    num_waypoints=1200, # Default
    terminal_states=None # Default
)

# %% 3D - Plot Palantir pseudotime
sc.pl.embedding(ad, basis="umap", color="palantir_pseudotime", frameon=False, save="_3D_UMAP_palantir_pseudotime.pdf")

# %% 3E - Signatures
# Signatures to keep:
# "Cytotoxicity",
# "ER1439_Selenocysteine_Synthesis_R-HSA-2408557",
# "ER231_T_Cell_Receptor_Signaling_Pathway_(GO0050852)",
# "ER781_Cholesterol_Biosynthetic_Process_(GO0006695)"

Firma_Cytotoxicity = ["PRF1", "IFNG", "NKG7", "GNLY", "GZMA", "GZMK", "GZMB", "GZMH", "CD7", "CCL5", "CD8A", "CCL4"]

Firma_Selenocysteine = [
    "RPS20", "SARS1", "RPL26L1", "RPL18", "RPL31", "RPS5", "RPL6", "RPLP0",
    "RPL3", "RPS16", "RPS19", "RPL18A", "RPL28", "RPL19", "RPL34", "SEPSECS",
    "RPS13", "RPS12", "RPL24", "RPS15", "RPL22", "RPS25", "RPL21", "RPL5",
    "RPS10", "RPL23", "RPS4Y1", "RPL36", "RPL27", "EEFSEC", "RPS15A", "RPL35",
    "RPS6", "RPLP1", "RPS24", "RPL3L", "RPS2", "RPS11", "RPL13A", "RPL11",
    "RPS8", "RPS27A", "RPL32", "RPS3A", "RPL37", "RPL10", "RPL7", "RPL7A",
    "RPS3", "FAU", "RPL30", "RPL8", "RPL26", "RPL29", "RPL22L1", "RPL9",
    "RPL39L", "RPS14", "RPL10L", "RPL36AL", "RPL27A", "RPL13", "RPSA", "RPS9",
    "RPS21", "RPS7", "RPL38", "RPL4", "RPL15", "RPLP2", "RPS27", "SEPHS2",
    "PSTK", "RPS17", "RPL35A", "RPS27L", "RPS23", "SECISBP2", "RPL14", "RPS26",
    "RPL37A", "RPL12", "RPS4X", "RPL23A", "RPL10A", "RPL39", "RPS29", "UBA52",
    "RPL41", "RPS18", "RPS28", "RPL36A", "RPL17", "RPS4Y2"
]

Firma_TCR_Signaling = [
    "ADA", "BTNL10P", "BTN3A3", "BTN2A2", "SIVA1", "RBCK1", "KHDRBS1", "CD226",
    "MALT1", "BTNL3", "LILRB4", "BTN3A2", "BTN3A1", "BTN2A1", "CD160", "HHLA2",
    "CD300A", "ERMAP", "SPPL3", "CCR7", "CRKL", "CSK", "RC3H1", "CTLA4",
    "ITPRIPL1", "BTNL9", "CYLD", "DENND1B", "DUSP3", "EIF2B1", "ELF1", "FYB2",
    "FCHO1", "RFTN1", "UBR2", "ICOSLG", "FOSL2", "CD2AP", "ABL1", "FYB1",
    "FYN", "ZNF683", "SLC39A6", "PRKD2", "PTPN22", "GATA3", "GBP1", "LAT",
    "TNFRSF21", "TRDC", "TRBC2", "TRBC1", "TRAC", "PHPT1", "STOML2", "HLA-A",
    "HLA-DPB1", "HLA-DQB1", "HLA-DRB1", "HLA-DRB3", "HRAS", "IKBKB", "INPP5D",
    "ITK", "KCNN4", "THEMIS", "LCK", "LCP2", "LGALS3", "LIPA", "SH2D1A",
    "MOG", "NCK1", "PAWR", "TRAT1", "FOXP3", "PDE4B", "PDE4D", "PIK3CA",
    "PIK3CD", "UBASH3A", "PLCG1", "PLCG2", "RC3H2", "BTN2A3P", "LIME1", "RNF31",
    "MAPK1", "PRNP", "BTNL2", "PSEN1", "DUSP22", "PTPN2", "PTPN6", "PTPRC",
    "PTPRJ", "NECTIN2", "RELA", "RPS3", "CEACAM1", "NFKBIZ", "SHB", "WNK1",
    "BRAF", "STK11", "MAP3K7", "BTK", "BTN1A1", "TRGC1", "TRGC2", "TEC",
    "THY1", "TRAF6", "TXK", "UBE2N", "EZR", "ZAP70", "LAPTM5", "CACNB3",
    "CACNB4", "PVRIG", "VTCN1", "BTNL8", "ZC3H12A", "CD276", "PRAM1", "SLA2",
    "CARD11", "NFKBID", "IKBKG", "DGKZ", "SKAP1", "CBLB", "RIPK2", "EIF2B4",
    "EIF2B3", "EIF2B2", "EIF2B5", "BCL10", "RAB29", "BTRC", "CD3D", "CD3E",
    "CD3G", "CD247", "CD8A", "CD8B", "CD28", "THEMIS2", "BCAR1", "CD81",
    "TESPA1"
]

Firma_Cholesterol_Biosynthesis = [
    "DHCR7", "ACAA2", "HSD17B7", "MSMO1", "LBR", "APOE", "INSIG2", "NPC1L1",
    "IDI2", "APOA5", "CYP51A1", "NSDHL", "DHCR24", "PMVK", "EBP", "IDI1",
    "PRKAA1", "MVK", "HMGCS1", "CNBP", "HMGCS2", "PRKAA2", "MVD", "ACLY",
    "LSS", "FDFT1", "CES1", "FDPS", "CFTR", "G6PD", "APOA4", "HMGCR",
    "APOA1", "CYB5R3", "TM7SF2", "INSIG1", "PRKAG2"
]

# Merge signatures
pathways = {
    "Cholesterol biosynthesis": Firma_Cholesterol_Biosynthesis,
    "TCR signaling": Firma_TCR_Signaling,
    "Cytotoxicity": Firma_Cytotoxicity,
    "Selenocysteine synthesis": Firma_Selenocysteine,
}

# %% 3E - Define variables
PSEUDOTIME_COL = "palantir_pseudotime"

# %% 3E - MANUAL gene selection
top_genes = {
    "Cholesterol biosynthesis": ["FDPS", "PMVK", "FDFT1", "DHCR24", "INSIG1", "DHCR7", "MSMO1", "ACAA2"],
    "TCR signaling": ["PTPRC", "CD3D", "CD28", "CD247", "LCK", "LCP2", "ZAP70", "LAT"],
    "Cytotoxicity": ["NKG7", "CCL5", "GZMH", "GZMA", "IFNG", "GZMB", "CD8A", "CCL4"],
    "Selenocysteine synthesis": ["SARS1", "SEPHS2", "EEFSEC", "SEPSECS", "PSTK", "SECISBP2", "RPL30", "RPL41"],
}

# %% 3E - Heatmaps per signature (top genes) along pseudotime
N_BINS = 25
VMIN, VMAX = -2, 2  # Color range for Z-score

pt = ad.obs[PSEUDOTIME_COL].astype(float).values

# Sort panels to match the story
panel_order = [
    "Cholesterol biosynthesis",
    "TCR signaling",
    "Cytotoxicity",
    "Selenocysteine synthesis",
]

fig, axes = plt.subplots(4, 1, figsize=(7.5, 9.5), constrained_layout=True)

last_im = None
for ax, pname in zip(axes, panel_order):
    genes = top_genes[pname]
    Mz, bin_mids = binned_means_zscore(ad, genes, pt, n_bins=N_BINS)

    last_im = ax.imshow(Mz, aspect="auto", interpolation="nearest", vmin=VMIN, vmax=VMAX)
    ax.set_title(pname)
    ax.set_yticks(np.arange(len(genes)))
    ax.set_yticklabels(genes, fontsize=8)

    n_eff = Mz.shape[1]
    ax.set_xlim(-0.5, n_eff - 0.5)
    ax.set_xticks([-0.5, n_eff - 0.5])
    ax.set_xticklabels(["0", "1"])
    ax.set_xlabel("Pseudotime")


cbar = fig.colorbar(last_im, ax=axes, fraction=0.02, pad=0.01)
cbar.set_label("Expression (z-score per gene)")

# Save
out_pdf = "_3E_Palantir_TopGenes_Heatmaps_4Pathways.pdf"
fig.savefig(out_pdf, bbox_inches="tight")
plt.show()

print("Saved:", out_pdf)

# %% End of script
