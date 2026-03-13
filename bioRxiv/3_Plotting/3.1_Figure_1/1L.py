###############################################################################
###############################################################################

# Program: 1L.py
# Author: Sergio Cámara Peña
# Date: 04/08/2025
# Version: FINAL

###############################################################################
###############################################################################

# %% Import libraries
import os
import scanpy as sc
import anndata
import numpy as np
import random
import pandas as pd

# %% Set random seed
random.seed(2504)

# %% Load data
path = "/home/scamara/data_a/scamara/Atlas/Input"
file = f"{path}/Python_scVI_adata_big_V4_state4.h5ad"
adata = anndata.read_h5ad(file)
print(adata.shape[0])

# Normalize
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Optionally, set raw to keep a copy of the normalized data
adata.raw = adata

# %% Filter object
# Antigen == "Blood"
adata_filtered = adata[adata.obs['Antigen'] == "Blood"].copy()
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

mask = ~((adata_filtered.obs['Stimulated'] == "YES") & (adata_filtered.obs['Stimulation_Location'] == "In_vitro"))
adata_filtered = adata_filtered[mask].copy()
print(adata_filtered.shape[0])

adata_filtered_bis = adata_filtered.copy()

# Max_Response == "CR"
adata_filtered = adata_filtered[adata_filtered.obs['Max_Response'] == "CR"].copy()
print(adata_filtered.shape[0])

# %% Signatures
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

# Merge
gene_signatures = {
    "Cytotoxicity": Firma_Cytotoxicity,
    "Selenocysteine_Synthesis": Firma_Selenocysteine,
    "T_Cell_Receptor_Signaling": Firma_TCR_Signaling,
    "Cholesterol_Biosynthesis": Firma_Cholesterol_Biosynthesis
}

# %% Calculate score for each signature and save in adata.obs
for name, genes in gene_signatures.items():
    genes_upper = [g.upper() for g in genes]
    var_names_upper = [g.upper() for g in adata_filtered.var_names]
    genes_present = [genes[i] for i, g in enumerate(genes_upper) if g in var_names_upper]

    if len(genes_present) == 0:
        print(f"No genes of signature {name} found in dataset — skipping.")
        continue

    print(f"Scoring signature {name} with {len(genes_present)} genes (out of {len(genes)} total)")

    sc.tl.score_genes(
        adata_filtered,
        gene_list=genes_present,
        score_name=name + "_score",
        use_raw=False
    )

# %% Plot
os.chdir("/home/scamara/data_a/scamara/Atlas/Figuras/Figura_1")
for name in gene_signatures.keys():
    score_col = name + "_score"
    if score_col in adata_filtered.obs.columns:
        sc.pl.umap(
            adata_filtered,
            color=score_col,
            title=f"Signature: {name}",
            vmin=0,
            show=False,
            save=f"_1L_{name}.pdf"
        )
    else:
        print(f"{score_col} not found in adata.obs, skipping plot.")

# %% End of script
