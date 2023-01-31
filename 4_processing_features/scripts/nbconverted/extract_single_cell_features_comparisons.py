# ---
# jupyter:
#   jupytext:
#     formats: ipynb,scripts/nbconverted//py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.4
#   kernelspec:
#     display_name: Python [conda env:4.processing_features_comparisons] *
#     language: python
#     name: conda-env-4.processing_features_comparisons-py
# ---

# ## Process single cell morphology features for CellProfiler readouts - Comparisons
#
# Compare the output of `pycytominer` and `pycytominer-transform` using NF1_SchwannCells_data output.

# +
import os
import pathlib
import warnings

import pandas as pd
from pycytominer import normalize
from pycytominer.cyto_utils import cells, output
from pycytominer_transform import config, convert

# ignore warnings
warnings.filterwarnings("ignore")

# +
# Set file and directory constants
cp_dir = "../CellProfiler_pipelines"
output_dir = "data"

# gather sorted list of source data filepaths
single_cell_filepaths = sorted(
    list(
        pathlib.Path(f"{cp_dir}/Analysis_Output/").glob(
            "**/NF1_data_allcp_plate*.sqlite"
        )
    )
)

# prepare platemap filepath
platemap_file = f"{cp_dir}/Metadata/platemap_NF1_CP.csv"

# prepare output data filepaths
sc_output_file = pathlib.Path(f"{output_dir}/nf1_sc_cellprofiler.csv.gz")
sc_norm_output_file = pathlib.Path(f"{output_dir}/nf1_sc_norm_cellprofiler.csv.gz")

print("single_cell_filepaths:")
single_cell_filepaths
# -

# Define custom linking columns between compartments
linking_cols = {
    "Per_Cytoplasm": {
        "Per_Cells": "Cytoplasm_Parent_Cells",
        "Per_Nuclei": "Cytoplasm_Parent_Nuclei",
    },
    "Per_Cells": {"Per_Cytoplasm": "Cells_Number_Object_Number"},
    "Per_Nuclei": {"Per_Cytoplasm": "Nuclei_Number_Object_Number"},
}

# Load platemap file
platemap_df = pd.read_csv(platemap_file)
platemap_df

# show file to be processed below
single_cell_filepath_one = single_cell_filepaths[0]
single_cell_filepath_one

# %%time
# pycytominer
# perform merge single cells without annotation
# and export to parquet format, re-reading the result
# from the parquet file for precision in comparison
pycytominer_sc_df_without_annotation = pd.read_parquet(
    # path is the output filename of SingleCells.merge_single_cells
    path=cells.SingleCells(
        sql_file=f"sqlite:///{single_cell_filepath_one}",
        compartments=["Per_Cells", "Per_Cytoplasm", "Per_Nuclei"],
        compartment_linking_cols=linking_cols,
        image_table_name="Per_Image",
        strata=["Image_Metadata_Well", "Image_Metadata_Plate"],
        merge_cols=["ImageNumber"],
        image_cols="ImageNumber",
        load_image_data=True,
        # perform merge_single_cells without annotation
        # and receive parquet filepath
    ).merge_single_cells(
        sc_output_file="pycytominer_singlecells_merge.parquet",
        output_type="parquet",
    )
)

pycytominer_sc_df_without_annotation.info()
pycytominer_sc_df_without_annotation.head()

# show configuration preset SQL for pycytominer-transform
print(config["cellprofiler_sqlite"]["CONFIG_JOINS"])

# make a copy of default configuration presets from pycytominer-transform
modified_config = config.copy()
# replace Cytoplasm_Parent_OrigNuclei value in modified config
modified_config["cellprofiler_sqlite"]["CONFIG_JOINS"] = modified_config[
    "cellprofiler_sqlite"
]["CONFIG_JOINS"].replace(
    "per_cytoplasm.Cytoplasm_Parent_OrigNuclei", "per_cytoplasm.Cytoplasm_Parent_Nuclei"
)
# show the modified configuration
print(modified_config["cellprofiler_sqlite"]["CONFIG_JOINS"])

# %%time
# pycytominer-transform
# perform merge without annotation and export
# to parquet format, reading the result
# from the parquet file for comparison
pycytominer_transform_sc_df_without_annotation = pd.read_parquet(
    path=convert(
        source_path=str(single_cell_filepath_one),
        dest_path="./pycytominer-transform_singlecells_merge.parquet",
        dest_datatype="parquet",
        merge=True,
        merge_chunk_size=100,
        preset="cellprofiler_sqlite",
        joins=modified_config["cellprofiler_sqlite"]["CONFIG_JOINS"],
    )
)

pycytominer_transform_sc_df_without_annotation.info()
pycytominer_transform_sc_df_without_annotation.head()

# check for missing cols from pycytominer to pycytominer-transform
pycytominer_missing_cols = [
    col
    for col in pycytominer_sc_df_without_annotation
    if col not in pycytominer_transform_sc_df_without_annotation.columns
]
pycytominer_missing_cols

# check for missing cols from pycytominer-transform to pycytominer
pt_missing_cols = [
    col
    for col in pycytominer_transform_sc_df_without_annotation.columns
    if col not in pycytominer_sc_df_without_annotation.columns
]
pt_missing_cols

# form renaming dictionary
rename_for_metadata_prefix = {
    orig: new
    for orig, new in zip(pt_missing_cols, pycytominer_missing_cols)
    if orig in new
}
rename_for_metadata_prefix

# minor rename for existing data
pycytominer_transform_sc_df_without_annotation = (
    pycytominer_transform_sc_df_without_annotation.rename(
        columns=rename_for_metadata_prefix
    )
)

# check for missing cols from pycytominer to pycytominer-transform
pycytominer_missing_cols = [
    col
    for col in pycytominer_sc_df_without_annotation
    if col not in pycytominer_transform_sc_df_without_annotation.columns
]
pycytominer_missing_cols

# check the shape of both dataframes
print(pycytominer_sc_df_without_annotation.shape)
print(pycytominer_transform_sc_df_without_annotation.shape)

# test dataframe equality
pd.testing.assert_frame_equal(
    left=pycytominer_sc_df_without_annotation,
    right=pycytominer_transform_sc_df_without_annotation[
        # use the pycytominer column order as a reference for pycytominer-transform output
        pycytominer_sc_df_without_annotation.columns
    ],
)

# +
# Merge single cells across compartments
anno_kwargs = {"join_on": ["Metadata_well_position", "Image_Metadata_Well"]}

sc_df = sc.merge_single_cells(
    platemap=platemap_df,
    **anno_kwargs,
)

# Save level 2 data as a csv
output(sc_df, sc_output_file)

print(sc_df.shape)
sc_df.head()

# +
# Normalize single cell data and write to file
normalize_sc_df = normalize(sc_df, method="standardize")

output(normalize_sc_df, sc_norm_output_file)

print(normalize_sc_df.shape)
normalize_sc_df.head()
# -

# ### Visualize basic count statistics

sc_df.Metadata_genotype.value_counts()

pd.crosstab(sc_df.Metadata_genotype, sc_df.Metadata_Well)
