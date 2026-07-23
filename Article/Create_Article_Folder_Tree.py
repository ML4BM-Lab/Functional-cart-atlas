###############################################################################
###############################################################################

# Program: Create_Article_Folder_Tree.py
# Author: Sergio Cámara Peña
# Date: 23/07/2026
# Version: FINAL
# Purpose: Create the directory tree used by the Article analysis scripts

###############################################################################
###############################################################################

# %% Command-line and environment path configuration

import argparse
import os
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent if "__file__" in globals() else Path.cwd()
_path_parser = argparse.ArgumentParser(
    description=(
        "Create missing directories used by the Article analysis scripts. "
        "Existing directories and files are never overwritten or deleted."
    )
)
_path_parser.add_argument(
    "--project-dir",
    default=os.environ.get("CART_ATLAS_PROJECT_DIR", str(_SCRIPT_DIR)),
    help=(
        "Root where the directory tree will be created "
        "(env: CART_ATLAS_PROJECT_DIR; default: directory containing this script)."
    ),
)
if "__file__" in globals():
    _path_args = _path_parser.parse_args()
else:
    _path_args, _unknown_path_args = _path_parser.parse_known_args()

project_dir = Path(_path_args.project_dir).expanduser().resolve()

# %% Fixed directory tree derived from the Article scripts

DIRECTORY_TREE = (
    ".huggingface",
    "Codigo/Codigo_datasets_atlas/Datasets_Integration/Initial_Version_Atlas",
    "Codigo/Codigo_datasets_atlas/Datasets_Integration/Tests_labs",
    "Codigo/Gene_Markers_Info",
    "Codigo/Rocinante_DEA/Funciones",
    "Datasets/Bai_et_al/Aligned_data_of_JITC_paper",
    "Datasets/Good_et_al/Downloaded_count_matrices",
    "Datasets/Haradhvala_et_al/Downloaded_data/Authors_Metadata_Reduced",
    "Datasets/Li_X_Cancer_Cell_letter_et_al",
    "Datasets/Li_X_Cancer_Cell_letter_et_al/Sample_files",
    "Datasets_Metadata",
    "Figuras/Figura_2",
    "Figuras/Figura_3",
    "Figuras/Figura_4",
    "Figuras/Figura_5",
    "Figuras/Suplementarias",
    "Gene_Markers_Info",
    "Input",
    "Models/CellTypist",
    "New_Datasets/Jordana_et_al/Datos_lorea",
    "Para_Cluster_Profiler",
    "Resultados/AUCell_CD8_IL10_CR",
    "Resultados/AUCell_CD8_IL10_NR",
    "Resultados/Bai_et_al/Plots/Exploratory_analysis",
    "Resultados/Bai_et_al/Plots/Worst_case_scenario",
    "Resultados/Bai_et_al/RDS",
    "Resultados/Boroughs_et_al/Count_Matrices/Cell_Ranger",
    "Resultados/Boroughs_et_al/Plots/Exploratory_analysis",
    "Resultados/Boroughs_et_al/Plots/Worst_case_scenario",
    "Resultados/Boroughs_et_al/RDS",
    "Resultados/CD8_IP_comparison",
    "Resultados/CD8_Short_Anytime_CR",
    "Resultados/CD8_clusters",
    "Resultados/Deng_et_al/Count_Matrices/Cell_Ranger",
    "Resultados/Deng_et_al/Plots/Exploratory_analysis",
    "Resultados/Deng_et_al/Plots/Worst_case_scenario",
    "Resultados/Deng_et_al/RDS",
    "Resultados/Dreamlet_Vs_PyDESeq2",
    "Resultados/Good_et_al/Plots/Exploratory_analysis",
    "Resultados/Good_et_al/Plots/Worst_case_scenario",
    "Resultados/Good_et_al/RDS",
    "Resultados/Haradvala_et_al",
    "Resultados/Haradvala_et_al/Count_Matrices/Cell_Ranger",
    "Resultados/Haradvala_et_al/Plots/Exploratory_analysis",
    "Resultados/Haradvala_et_al/Plots/Worst_case_scenario",
    "Resultados/Haradvala_et_al/RDS",
    "Resultados/IACs_vs_Rest",
    "Resultados/IP_comparison",
    "Resultados/Joined_datasets/Dreamlet_to_ClusterProfiler",
    "Resultados/Joined_datasets/Integration/Plots/V4/Annotation/Low_Res",
    "Resultados/Joined_datasets/Integration/Plots/V4/Colorized_Atlas_Metadata",
    "Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/CD4",
    "Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/CD4/GSA",
    "Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/CD8",
    "Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/CD8/GSA",
    "Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/Final",
    "Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/GATA3",
    "Resultados/Joined_datasets/Integration/Plots/V4/Subclustering/GATA3/GSA",
    "Resultados/Joined_datasets/Integration/Plots/V4/WO_integ",
    "Resultados/Joined_datasets/Integration/Plots/V4/scVI_integ",
    "Resultados/Joined_datasets/Integration/Plots/V4/scVI_integ/figures",
    "Resultados/Joined_datasets/Integration/Plots/V4/scVI_markers_exploration",
    "Resultados/Joined_datasets/Integration/Plots/V4/scVI_markers_exploration_2",
    "Resultados/Joined_datasets/Integration/Plots/V4/scVI_markers_exploration_3",
    "Resultados/Joined_datasets/Integration/Plots/V5/WO_integ",
    "Resultados/Joined_datasets/Integration/Plots/WO_integ",
    "Resultados/Joined_datasets/Integration/Python-Celltypist/V4",
    "Resultados/Joined_datasets/Integration/Python-Celltypist/V5",
    "Resultados/Joined_datasets/Integration/RDS/WO_integ/V4",
    "Resultados/Joined_datasets/Integration/RDS/WO_integ/V5",
    "Resultados/Joined_datasets/Integration/Subclustering_V4",
    "Resultados/Joined_datasets/Integration/Subclustering_V4/Clustree_res",
    "Resultados/Joined_datasets/Integration/Subclustering_V4/Data/CD4",
    "Resultados/Joined_datasets/Integration/Subclustering_V4/Data/CD8",
    "Resultados/Joined_datasets/Integration/Subclustering_V4/Data/GATA3",
    "Resultados/Joined_datasets/Integration/scArches/V5",
    "Resultados/Joined_datasets/Integration/scVI/V4",
    "Resultados/Joined_datasets/Integration/scVI/V4/Clustree_resolutions",
    "Resultados/Joined_datasets/Integration/scVI/V4/Data",
    "Resultados/Joined_datasets/Integration/scVI/V5",
    "Resultados/Joined_datasets/Integration_methods_lab",
    "Resultados/Joined_datasets/Integration_methods_lab/Harmony/Sin_GT",
    "Resultados/Joined_datasets/Integration_methods_lab/LIGER/Sin_GT",
    "Resultados/Joined_datasets/Integration_methods_lab/Merged_WO_integration/Sin_GT",
    "Resultados/Joined_datasets/Integration_methods_lab/STACAS/Sin_GT",
    "Resultados/Joined_datasets/Integration_methods_lab/Seurat_RPCA/Sin_GT",
    "Resultados/Joined_datasets/Integration_methods_lab/fastMNN/Sin_GT",
    "Resultados/Joined_datasets/Integration_methods_lab/scGraph",
    "Resultados/Joined_datasets/Integration_methods_lab/scIB_Merged",
    "Resultados/Joined_datasets/Integration_methods_lab/scVI/Sin_GT_With_Python",
    "Resultados/Joined_datasets/Post_Integration/Data",
    "Resultados/Joined_datasets/Raw_Atlas",
    "Resultados/Jordana_et_al/Plots/Worst_case_scenario",
    "Resultados/Jordana_et_al/RDS",
    "Resultados/Li_X_Cancer_Cell_letter_et_al/Plots/Exploratory_analysis",
    "Resultados/Li_X_Cancer_Cell_letter_et_al/Plots/Worst_case_scenario",
    "Resultados/Li_X_Cancer_Cell_letter_et_al/RDS",
    "Resultados/Li_X_et_al/Count_Matrices/Cell_Ranger",
    "Resultados/Li_X_et_al/Plots/Exploratory_analysis",
    "Resultados/Li_X_et_al/Plots/Worst_case_scenario",
    "Resultados/Li_X_et_al/RDS",
    "Resultados/Lynn_et_al/Count_Matrices/Cell_Ranger/CD19_vs_GD2_exp",
    "Resultados/Lynn_et_al/Count_Matrices/Cell_Ranger/Control_vs_JUN",
    "Resultados/Lynn_et_al/Plots/Exploratory_analysis",
    "Resultados/Lynn_et_al/Plots/Worst_case_scenario",
    "Resultados/Lynn_et_al/RDS",
    "Resultados/Melenhorst_et_al/Count_Matrices/Cell_Ranger",
    "Resultados/Melenhorst_et_al/Plots/Exploratory_analysis",
    "Resultados/Melenhorst_et_al/Plots/Worst_case_scenario",
    "Resultados/Melenhorst_et_al/RDS",
    "Resultados/Rodriguez-Marquez_et_al/Count_Matrices/Cell_Ranger",
    "Resultados/Rodriguez-Marquez_et_al/Plots/Exploratory_analysis",
    "Resultados/Rodriguez-Marquez_et_al/Plots/Worst_case_scenario",
    "Resultados/Rodriguez-Marquez_et_al/RDS",
    "Resultados/Sheih_et_al/Count_Matrices/Cell_Ranger",
    "Resultados/Sheih_et_al/Plots/Exploratory_analysis",
    "Resultados/Sheih_et_al/Plots/Worst_case_scenario",
    "Resultados/Sheih_et_al/RDS",
    "Resultados/Wang_et_al/Count_Matrices/Cell_Ranger",
    "Resultados/Wang_et_al/Plots/Exploratory_analysis",
    "Resultados/Wang_et_al/Plots/Worst_case_scenario",
    "Resultados/Wang_et_al/RDS",
    "Resultados/Xhangolli_et_al/Count_Matrices/Dropseq",
    "Resultados/Xhangolli_et_al/Plots/Cancer_cells_Removal",
    "Resultados/Xhangolli_et_al/Plots/Exploratory_analysis",
    "Resultados/Xhangolli_et_al/Plots/Worst_case_scenario",
    "Resultados/Xhangolli_et_al/RDS",
    "Resultados_Figuras/Data",
    "Resultados_Figuras/Figura_1",
    "Resultados_Figuras/Figura_2",
    "Resultados_Figuras/Figura_3",
    "Resultados_Figuras/Figura_4",
    "Resultados_Figuras/Figura_5",
    "Resultados_Figuras/Suplementarias",
    "Resultados_V5/BCMA_vs_CD19_IP_All",
    "Resultados_V5/BCMA_vs_CD19_MID",
    "Resultados_V5/Supplementary_Table_13",
    "Shiny_app",
    "Shiny_app/shinyAtlas",
    "UCSC_Browser",
    "cache_scvi",
)

# %% Create missing directories without modifying existing content

def create_directory_tree(root, relative_directories):
    created = 0
    existing = 0
    conflicts = []
    if root.exists() and not root.is_dir():
        raise NotADirectoryError(f"Project root exists but is not a directory: {root}")
    root.mkdir(parents=True, exist_ok=True)
    for relative_directory in relative_directories:
        directory = root / relative_directory
        if directory.is_dir():
            existing += 1
            continue
        if directory.exists():
            conflicts.append(f"{directory} (a non-directory entry already exists)")
            continue
        try:
            directory.mkdir(parents=True, exist_ok=True)
            created += 1
        except OSError as error:
            conflicts.append(f"{directory} ({error})")
    return created, existing, conflicts


if __name__ == "__main__":
    created_count, existing_count, directory_conflicts = create_directory_tree(
        project_dir,
        DIRECTORY_TREE,
    )
    print(f"Project directory: {project_dir}")
    print(f"Directory-tree entries created: {created_count}")
    print(f"Directory-tree entries already present: {existing_count}")
    if directory_conflicts:
        print("Directory-tree entries not created because of path conflicts:", file=sys.stderr)
        for conflict in directory_conflicts:
            print(f"  {conflict}", file=sys.stderr)
        raise SystemExit(1)
