###############################################################################
###############################################################################

# Program: scVI-hub_Download.py
# Author: Sergio Cámara Peña
# Date: 27/09/2025
# Version: FINAL
# Machine: Margaret

###############################################################################
###############################################################################

# %% Command-line and environment path configuration

import argparse
import os
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
_ARTICLE_DIR = next(
    (candidate for candidate in (_SCRIPT_DIR, *_SCRIPT_DIR.parents) if candidate.name == "Article"),
    _SCRIPT_DIR / "Article" if (_SCRIPT_DIR / "Article").is_dir() else _SCRIPT_DIR,
)
_path_parser = argparse.ArgumentParser()
_path_parser.add_argument(
    "--cache-dir",
    default=os.environ.get("CART_ATLAS_CACHE_DIR", str(_ARTICLE_DIR / "cache_scvi")),
    help="Directory used for the downloaded scVI hub cache.",
)
_path_args, _unknown_path_args = _path_parser.parse_known_args()
cache_dir = Path(_path_args.cache_dir).expanduser().resolve()

# %% Load needed libraries
import shutil
from scvi.hub import HubModel

# %% Define cache directory
# If cache already exists, remove it
if cache_dir.exists():
    shutil.rmtree(cache_dir)

# Recreate empty cache directory
cache_dir.mkdir(parents=True, exist_ok=True)

# %% Download model
hmo = HubModel.pull_from_huggingface_hub(
    repo_name="sergiocamarap/Functional-cart-atlas-model",
    revision="main",
    cache_dir=cache_dir
)

hmo.model

# %% End of script
