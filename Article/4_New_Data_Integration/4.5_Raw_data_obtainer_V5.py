###############################################################################
###############################################################################

# Program: 4.5_Raw_data_obtainer_V5.py
# Author: Sergio Cámara Peña
# Date: 21/11/2024
# Version: V5
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
    "--project-dir",
    default=os.environ.get("CART_ATLAS_PROJECT_DIR", str(_ARTICLE_DIR)),
    help="Root containing the repository-relative data, results, figure, model, and metadata subdirectories.",
)
_path_args, _unknown_path_args = _path_parser.parse_known_args()
project_dir = Path(_path_args.project_dir).expanduser().resolve()
if not project_dir.is_dir():
    raise FileNotFoundError(
        f"Project directory does not exist: {project_dir}. "
        "Set --project-dir or CART_ATLAS_PROJECT_DIR."
    )

_generated_input_paths = {
    project_dir / "Resultados" / "Joined_datasets" / "Raw_Atlas" / "Atlas_integ_scArches_FINAL_V5.h5ad",
}

def _require_path(path):
    if not path.exists():
        raise FileNotFoundError(f"Required input path does not exist: {path}")
    return path


def _input_path(directory, filename):
    path = directory / filename
    if path.exists():
        return path
    if path in _generated_input_paths:
        input_path = project_dir / "Input" / path.name
        if input_path.exists():
            return input_path
        raise FileNotFoundError(
            "Required generated input file does not exist. Checked: "
            + ", ".join(str(candidate) for candidate in (path, input_path))
        )
    raise FileNotFoundError(f"Required input path does not exist: {path}")


def _output_path(directory, filename):
    directory.mkdir(parents=True, exist_ok=True)
    return directory / filename

# %% Import libraries
import os
import scanpy as sc

# %% Switches
Guardar = True

# %% Read files - V5
_current_dir = Path(project_dir / 'Resultados' / 'Joined_datasets' / 'Raw_Atlas')
adata_V5 = sc.read_h5ad(_input_path(_current_dir, "Atlas_integ_scArches_FINAL_V5.h5ad"))

# %% Save files normalized - V5
_current_dir = Path(project_dir / 'Resultados' / 'Joined_datasets' / 'Raw_Atlas')
if Guardar:
    sc.pp.normalize_total(adata_V5, target_sum=1e4)
    sc.pp.log1p(adata_V5)
    adata_V5.raw = adata_V5.copy()
    adata_V5.write(_output_path(_current_dir, "Atlas_integ_scArches_FINAL_V5_Normalized.h5ad"))
else:
    adata_V5 = sc.read_h5ad(_input_path(_current_dir, "Atlas_integ_scArches_FINAL_V5_Normalized.h5ad"))

# %% End of script
