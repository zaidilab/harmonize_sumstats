from io_utils import load_ref_data, load_gwas_data, write_output
from liftover import perform_liftover
from join_utils import join_and_filter
from transformation import apply_transformations
from constants import OUTPUT_FILENAME
import argparse, os
import polars as pl

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--genome-build",
        type=str,
        required=True,
        choices=["hg19", "hg38"],
        help="Genome build of the GWAS file (e.g: hg19, hg38)"
    )
    parser.add_argument(
        "--ip-path",
        type=str,
        required=True,
        help="File path to input GWAS dataset"
    )
    parser.add_argument(
        "--op-path",
        type=str,
        required=True,
        help="File path to save output file"
    )

    args = parser.parse_args()

    # Load the reference dataset
    ref_df = load_ref_data()
    print(f"Reference SNP list preview: {ref_df.head()}")

    # Load the GWAS dataset
    gwas_df = load_gwas_data(args.ip_path)
    print(f"GWAS summary statistics preview: {gwas_df.head()}")

    # If the input genome build is not hg38, perform the liftover
    if args.genome_build != "hg38":
        gwas_df = perform_liftover(gwas_df, args.genome_build)
    else:
        # If no liftover, ensure '_hg38' column is present and same as 'base_pair_location'
        gwas_df = gwas_df.with_columns(pl.col("base_pair_location").alias("_hg38"))
    print(f"GWAS DataFrame after liftover {gwas_df.head()}")

    # Prepare for Join: Cast 'chromosome' column in 'gwas_df' to Utf8 
    # This ensures matching data types for join keys between gwas_df and ref_df.
    gwas_df = gwas_df.with_columns(pl.col("chromosome").cast(pl.Utf8))

    # Perform inner join and filtering of the DataFrames
    joined_df = join_and_filter(gwas_df, ref_df)

    # Apply the transformations to beta and frequency when necessary
    joined_df = apply_transformations(joined_df)

    # Drop unnecessary columns
    joined_df = joined_df.drop(["base_pair_location", "ref", "alt", "id"])
    print(f"Updated GWAS DataFrame preview: {joined_df.head()}")

    # Write the Output DataFrame to a .tsv file 
    # Ensure the output directory exists
    os.makedirs(args.op_path, exist_ok=True)
    write_output(joined_df, os.path.join(args.op_path, OUTPUT_FILENAME))

if __name__ == "__main__":
    main()