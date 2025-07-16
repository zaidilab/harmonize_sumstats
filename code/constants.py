import os
import polars as pl

# File path to the reference dataset
REF_SNP_LIST_PATH = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", "data", "1kg_hg38_hm3.snplist")
)

# Schema Definition for Reference SNP list with column names and their datatypes
REF_DF_SCHEMA = {
    "chromosome": pl.Utf8,
    "position": pl.Int64,
    "id": pl.Utf8,
    "ref": pl.Utf8,
    "alt": pl.Utf8
}

# Name of the file which would have the new, processed GWAS file
OUTPUT_FILENAME = "GWAS_new_1M.tsv"