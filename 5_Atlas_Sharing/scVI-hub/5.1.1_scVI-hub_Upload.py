###############################################################################
###############################################################################

# Program: 5.1.1_scVI-hub_Upload.py
# Author: Sergio Cámara Peña
# Date: 27/09/2025
# Version: FINAL

###############################################################################
###############################################################################

# %% Load needed libraries
import scvi
import scanpy as sc
from scvi.hub import HubMetadata, HubModel, HubModelCardHelper
import os

# %% 0. Load adata used to train model
First_Time = False
os.chdir("/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/Python-Celltypist/V4/")
if First_Time:
    adata = sc.read(
        "Seurat_merged_With_Celltypist.h5ad"
    )

    adata.layers["counts"] = adata.X.copy()

    ## Normalize and log scale
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata  # Freeze the state in ".raw"

    ## Remove samples with fewer than 45 CD3⁺ and CAR⁺ cells — since they cause issues in the next step.
    sample_column = 'Product_norm'

    # Count the number of cells per sample
    cell_counts = adata.obs[sample_column].value_counts(sort=False)

    # Filter samples with less than 45 cells
    filtered_samples = cell_counts[cell_counts >= 45].index

    # Subset the AnnData object to keep only samples with at least 45 cells
    adata_filtered = adata.copy()
    adata_filtered = adata_filtered[adata_filtered.obs[sample_column].isin(filtered_samples)]

    del adata
    adata = adata_filtered.copy()
    del adata_filtered

    # Identify and only keep the 2000 most variable genes
    sc.pp.filter_genes(adata, min_cells=10)

    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        n_top_genes=2000,
        layer="counts",
        batch_key="Product_norm",
        span=0.6,
        subset=True,
    )

    adata.obs.rename(columns={"orig.ident": "orig_ident"}, inplace=True)

    adata.write("scVI_hub_adata.h5ad")

else:
    adata = sc.read_h5ad("scVI_hub_adata.h5ad")

# %% 1. Load pretrained model
os.chdir(
    "/mnt/md0/data/scamara/Atlas_Mieloma_Multiple/Resultados/Joined_datasets/Integration/scVI/V4"
)

local_dir = "Atlas_integ_scVI_V4"
model = scvi.model.SCVI.load(local_dir, adata=adata)

# %% 2. Create metadata for the model
hm = HubMetadata.from_dir(
    local_dir,
    anndata_version="0.9.2",
    training_data_url="https://zenodo.org/records/17213453/files/scVI_hub_adata.h5ad?download=1"  # link to Zenodo with the full AnnData
)

# %% 3. Create the Model Card from scvi-tools template
hmch = HubModelCardHelper.from_dir(
    local_dir,
    license_info="cc-by-4.0",
    anndata_version="0.9.2",
    data_modalities=["rna"],
    data_is_annotated=False,
    description="Pretrained scVI integration model for the Functional CAR-T Cell Atlas",
    references="Reference pending - manuscript in preparation.",
    training_data_url="https://doi.org/10.5281/zenodo.17213452",  # link to Zenodo with the full AnnData
    training_code_url="https://github.com/ML4BM-Lab/Functional-cart-atlas",
    data_is_minified=True
)

print(hmch.model_card.content)

# (optional) Save the model card to disk if you want to manually edit/extend it
hmch.model_card.save("Functional_CART_Atlas_model_card.md")

# %% 4. Create the HubModel object - AnnData NOT included (too large)
hmo = HubModel(local_dir, metadata=hm, model_card=hmch)
hmo

# %% 5. Push to Hugging Face Hub
# Load Hugging Face token
with open(os.path.expanduser("/home/scamara/data/scamara/.huggingface/token.txt"), "r") as f:
    repo_token = f.read().strip()

# Push
hmo.push_to_huggingface_hub(
    repo_name="sergiocamarap/Functional-cart-atlas-model",
    repo_create=False,
    repo_token=repo_token
)

# %% End of script
