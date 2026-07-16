###############################################################################
###############################################################################

# Program: 4.5_Raw_data_obtainer_V5.py
# Author: Sergio Cámara Peña
# Date: 21/11/2024
# Version: FINAL

###############################################################################
###############################################################################

# %% Import libraries
import os
import scanpy as sc

# %% Switches
Guardar = True

# %% Read files - V5
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
adata_V5 = sc.read_h5ad("Atlas_integ_scArches_FINAL_V5.h5ad")

# %% Save files normalized - V5
os.chdir("/home/scamara/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Raw_Atlas")
if Guardar:
    sc.pp.normalize_total(adata_V5, target_sum=1e4)
    sc.pp.log1p(adata_V5)
    adata_V5.raw = adata_V5.copy()
    adata_V5.write("Atlas_integ_scArches_FINAL_V5_Normalized.h5ad")
else:
    adata_V5 = sc.read_h5ad("Atlas_integ_scArches_FINAL_V5_Normalized.h5ad")

# %% End of script
