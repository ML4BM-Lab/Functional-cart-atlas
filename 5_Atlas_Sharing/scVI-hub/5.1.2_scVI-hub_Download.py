###############################################################################
###############################################################################

# Program: 5.1.2_scVI-hub_Download.py
# Author: Sergio Cámara Peña
# Date: 27/09/2025
# Version: FINAL

###############################################################################
###############################################################################

# %% Load needed libraries
import os
import shutil
from scvi.hub import HubModel

# %% Define cache directory
cache_dir = os.path.join(os.getcwd(), "cache_scvi")

# If cache already exists, remove it
if os.path.exists(cache_dir):
    shutil.rmtree(cache_dir)

# Recreate empty cache directory
os.makedirs(cache_dir, exist_ok=True)

# %% Download model
hmo = HubModel.pull_from_huggingface_hub(
    repo_name="sergiocamarap/Functional-cart-atlas-model",
    revision="main",
    cache_dir=cache_dir
)

hmo.model

# %% End of script
