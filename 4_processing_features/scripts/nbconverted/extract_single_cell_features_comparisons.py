# ---
# jupyter:
#   jupytext:
#     formats: ipynb,scripts/nbconverted//py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python [conda env:4.process-nf1-features_dev] *
#     language: python
#     name: conda-env-4.process-nf1-features_dev-py
# ---

# ## Process single cell morphology features for CellProfiler readouts - Comparisons
#
# Compare the output of `pycytominer` and `pycytominer-transform` using NF1_SchwannCells_data output.

# +
import os
import pathlib
import warnings

import duckdb
import pandas as pd
import pyarrow.parquet as parquet
from pycytominer import normalize
from pycytominer.cyto_utils import cells, output
from pycytominer_transform import convert

# ignore warnings
warnings.filterwarnings("ignore")

# +
# Set file and directory constants
cp_dir = "../CellProfiler_pipelines"
output_dir = "data"

sql_file = "NF1_data.sqlite"
single_cell_filepath = f"{cp_dir}/Analysis_Output/{sql_file}"
single_cell_file = f"sqlite:///{cp_dir}/Analysis_Output/{sql_file}"
platemap_file = f"{cp_dir}/Metadata/platemap_NF1_CP.csv"

sc_output_file = pathlib.Path(f"{output_dir}/nf1_sc_cellprofiler.csv.gz")
sc_norm_output_file = pathlib.Path(f"{output_dir}/nf1_sc_norm_cellprofiler.csv.gz")
# -

# Define custom linking columns between compartments
linking_cols = {
    "Per_Cytoplasm": {
        "Per_Cells": "Cytoplasm_Parent_Cells",
        "Per_Nuclei": "Cytoplasm_Parent_OrigNuclei",
    },
    "Per_Cells": {"Per_Cytoplasm": "Cells_Number_Object_Number"},
    "Per_Nuclei": {"Per_Cytoplasm": "Nuclei_Number_Object_Number"},
}

# Load platemap file
platemap_df = pd.read_csv(platemap_file)
platemap_df

# pycytominer
# perform merge single cells without annotation
# and export to parquet format, re-reading the result
# from the parquet file for precision in comparison
pycytominer_sc_df_without_annotation = pd.read_parquet(
    path=cells.SingleCells(
        sql_file=single_cell_file,
        compartments=["Per_Cells", "Per_Cytoplasm", "Per_Nuclei"],
        compartment_linking_cols=linking_cols,
        image_table_name="Per_Image",
        strata=["Image_Metadata_Well", "Image_Metadata_Plate"],
        merge_cols=["ImageNumber"],
        image_cols="ImageNumber",
        load_image_data=True,
        # perform merge_single_cells without annotation
        # and send ask for parquet based output, returning a filepath
    ).merge_single_cells(
        sc_output_file="pycytominer_singlecells_merge.parquet",
        output_type="parquet",
    )
)
pycytominer_sc_df_without_annotation.info()
pycytominer_sc_df_without_annotation.head()

# pycytominer-transform
# perform merge without annotation and export
# to parquet format, reading the result
# from the parquet file for comparison
pycytominer_transform_sc_df_without_annotation = pd.read_parquet(
    path=convert(
        source_path=single_cell_filepath,
        dest_path="./pycytominer-transform_singlecells_merge.parquet",
        dest_datatype="parquet",
        merge=True,
        merge_chunk_size=100,
        preset="cellprofiler_sqlite",
    )
)
pycytominer_transform_sc_df_without_annotation.info()
pycytominer_transform_sc_df_without_annotation.head()

# check for missing cols from pycytominer to pycytominer-transform
[
    col
    for col in pycytominer_sc_df_without_annotation
    if col not in pycytominer_transform_sc_df_without_annotation.columns
]

# check for missing cols from pycytominer-transform to pycytominer
[
    col
    for col in pycytominer_transform_sc_df_without_annotation.columns
    if col not in pycytominer_sc_df_without_annotation.columns
]

# minor rename for existing data
pycytominer_transform_sc_df_without_annotation = (
    pycytominer_transform_sc_df_without_annotation.rename(
        columns={
            "Cytoplasm_Parent_Cells": "Metadata_Cytoplasm_Parent_Cells",
            "Cytoplasm_Parent_OrigNuclei": "Metadata_Cytoplasm_Parent_OrigNuclei",
        }
    )
)

# append columns which already exist but are differently named
pycytominer_transform_sc_df_without_annotation[
    "Metadata_Cells_Number_Object_Number"
] = pycytominer_transform_sc_df_without_annotation["Metadata_Cytoplasm_Parent_Cells"]
pycytominer_transform_sc_df_without_annotation[
    "Metadata_Nuclei_Number_Object_Number"
] = pycytominer_transform_sc_df_without_annotation[
    "Metadata_Cytoplasm_Parent_OrigNuclei"
]

# check for missing cols from pycytominer to pycytominer-transform (after column changes)
[
    col
    for col in pycytominer_sc_df_without_annotation
    if col not in pycytominer_transform_sc_df_without_annotation.columns
]

# check for missing cols from pycytominer-transform to pycytominer (after column changes)
[
    col
    for col in pycytominer_transform_sc_df_without_annotation.columns
    if col not in pycytominer_sc_df_without_annotation.columns
]

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
# create a table from duckdb read of sqlite tables
table = (
    duckdb.connect()
    .execute(
        f"""
            /* install and load sqlite plugin for duckdb */
            INSTALL sqlite_scanner;
            LOAD sqlite_scanner;
            
            /* attach sqlite db to duckdb for full table awareness */
            CALL sqlite_attach('{single_cell_filepath}');
            
            /* perform query on sqlite tables through duckdb */
            with Per_Image_Filtered as 
                (select ImageNumber, Image_Metadata_Well, Image_Metadata_Plate from Per_Image)
            SELECT * from Per_Cytoplasm cytoplasm
            left join Per_Cells cells on
                cells.ImageNumber = cytoplasm.ImageNumber
                and cells.Cells_Number_Object_Number = cytoplasm.Cytoplasm_Parent_Cells
            left join Per_Nuclei nuclei on
                nuclei.ImageNumber = cytoplasm.ImageNumber
                and nuclei.Nuclei_Number_Object_Number = cytoplasm.Cytoplasm_Parent_OrigNuclei
            left join Per_Image_Filtered image on
                image.ImageNumber = cytoplasm.ImageNumber
            """
    )
    .arrow()
    .drop_null()
)

# account for duplicate column names from joins
cols = []
# reversed order column check as col removals will change index order
for i, colname in reversed(list(enumerate(table.column_names))):
    if colname not in cols:
        cols.append(colname)
    else:
        table = table.remove_column(i)

# write the arrow table to parquet
parquet.write_table(
    table=table,
    where="./duckdb_singlecells_merge.parquet",
)

# read the parquet file
duckdb_sc_df_without_annotation = pd.read_parquet(
    path="./duckdb_singlecells_merge.parquet"
)
duckdb_sc_df_without_annotation.info()
duckdb_sc_df_without_annotation.head()

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
