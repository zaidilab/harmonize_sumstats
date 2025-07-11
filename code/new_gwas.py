import polars as pl
from pyliftover import LiftOver
import time
import os
import argparse 

print("Script started.")

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

genome_build = args.genome_build
gwas_summary_path = args.ip_path
op_path = args.op_path

ref_snp_list_path = "../dataset/1kg_hg38_hm3.snplist"
output_filename = "GWAS_new_1M.tsv"

# Schema Definition for Reference SNP list
ref_df_schema = {
    "chromosome": pl.Utf8,
    "position": pl.Int64,
    "id": pl.Utf8,
    "ref": pl.Utf8,
    "alt": pl.Utf8
}

# Read Reference SNP list File
print(f"Loading reference SNP list")
ref_df = pl.read_csv(
    ref_snp_list_path,
    separator="\t",
    has_header=False,
    new_columns=["chromosome", "position", "id", "ref", "alt"],
    schema=ref_df_schema
)
print(f"Reference SNP list preview: {ref_df.head()}")

# Read GWAS Summary Statistics File 
print(f"Loading GWAS summary statistics")
gwas_df = pl.read_csv(
    gwas_summary_path,
    separator="\t"
)
print(f"GWAS summary statistics preview: {gwas_df.head()}")

# Conditional Liftover Operation 
# Polars' 'map_batches' function will call this for each chunk of the DataFrame.
def liftover_batch(df):
    """
    Python function to perform the liftover operation on a batch (chunk) of data.
    'df': a Polars DataFrame chunk.
    """
    # Initialize an empty list to store the lifted coordinates.
    lifted = []

    # Initialize Liftover object for the current batch.
    lo = LiftOver(genome_build, "hg38")

    # Iterate over pairs of 'chromosome' and 'base_pair_location' from the input DataFrame chunk.
    for chrom, pos in zip(df["chromosome"], df["base_pair_location"]):
        # Ensure chromosome is a string and add 'chr' prefix for pyliftover
        chrom_str = str(chrom) if not isinstance(chrom, str) else chrom
        result = lo.convert_coordinate(f"chr{chrom_str}", pos)

        # Append the converted location to the 'lifted' list.
        # 'result[0][1]' extracts the base pair location from the first successful conversion.
        # If 'result' is empty (meaning conversion failed), append None.
        lifted.append(result[0][1] if result else None)

    # Add the newly computed '_hg38' column to the original DataFrame chunk.
    return df.with_columns(pl.Series("_hg38", lifted, dtype=pl.Int64))

# Define the expected schema after liftover_batch, including the new '_hg38' column.
output_schema_after_liftover = {
    "chromosome": pl.Int64,
    "base_pair_location": pl.Int64,
    "effect_allele": pl.String,
    "other_allele": pl.String,
    "beta": pl.Float64,
    "standard_error": pl.Float64,
    "effect_allele_frequency": pl.Float64,
    "p_value": pl.Float64,
    "variant_id": pl.String,
    "rs_id": pl.String,
    "n": pl.Int64,
    "CHISQ": pl.Float64,
    "_hg38": pl.Int64
}

# Only perform liftover if genome_build is not hg38
if genome_build != "hg38":
    print(f"Performing liftover from {genome_build} to hg38...")
    start_time = time.time()
    # Convert the 'gwas_df' into a LazyFrame (operations are not executed immediately but are instead recorded as a plan) 
    # because map_batches() method are for lazyFrames only
    gwas_df = (gwas_df.lazy()
                .map_batches(liftover_batch, schema=output_schema_after_liftover)       # Applies the 'liftover_batch' Python function to chunks (batches) of the LazyFrame.
                .collect())                                                             # Executes the LazyFrame's plan and returns the result as an eager DataFrame.
    end_time = time.time()
    print(f"Liftover execution time: {end_time - start_time:.4f} seconds")
else:
    print("Input genome build is hg38, skipping liftover.")
    # If no liftover, ensure '_hg38' column is present and same as 'base_pair_location'
    gwas_df = gwas_df.with_columns(
        pl.col("base_pair_location").alias("_hg38")
    )
print(f"GWAS DataFrame after liftover {gwas_df.head()}")


# Prepare for Join: Cast 'chromosome' column in 'gwas_df' to Utf8 
# This ensures matching data types for join keys between gwas_df and ref_df.
gwas_df = gwas_df.with_columns(
    pl.col("chromosome").cast(pl.Utf8).alias("chromosome")
)

# Perform Inner Join 
joined_df = gwas_df.join(
    ref_df,                                         # The right-hand side DataFrame for the join.
    left_on=["chromosome", "_hg38"],                # Columns from 'gwas_df' to use as join keys.
    right_on=["chromosome", "position"],            # Columns from 'ref_df' to use as join keys.
    how="inner"                                     # 'inner' join to only include rows where keys exist in BOTH DataFrames.

# After joining, filter the joined DataFrame where,
# effect allele and other allele from the GWAS summary statistics
# match the reference and alternate alleles from the reference SNP list, respectively.
).filter(
    # Filter for rows where GWAS alleles match reference alleles directly OR are flipped
    ((pl.col("effect_allele") == pl.col("ref")) &
        (pl.col("other_allele") == pl.col("alt"))) |
    ((pl.col("effect_allele") == pl.col("alt")) &
        (pl.col("other_allele") == pl.col("ref")))
)
print(f"Joined and filtered DataFrame preview: {joined_df.head()}")

# Apply transformations (frequency, beta, allele swapping) only if the alleles need flipping.
joined_df = joined_df.with_columns([
    # Conditional transformation for 'effect_allele_frequency'
    pl.when(
        (pl.col("effect_allele") == pl.col("ref")) & # If GWAS effect_allele matches ref_df alt allele
        (pl.col("other_allele") == pl.col("alt"))    # AND GWAS other allele matches ref_df ref allele
    )
    .then(1 - pl.col("effect_allele_frequency")) # THEN, flip the frequency
    .otherwise(pl.col("effect_allele_frequency")) # ELSE, keep original frequency
    .alias("effect_allele_frequency"),

    # Conditional transformation for 'beta'
    pl.when(
        (pl.col("effect_allele") == pl.col("ref")) &
        (pl.col("other_allele") == pl.col("alt"))
    )
    .then(-pl.col("beta")) # THEN, flip the sign of beta
    .otherwise(pl.col("beta")) # ELSE, keep original beta
    .alias("beta"),

    # Conditional transformation for 'effect_allele'
    pl.when(
        (pl.col("effect_allele") == pl.col("ref")) &
        (pl.col("other_allele") == pl.col("alt"))
    )
    .then(pl.col("ref")) # THEN, set effect_allele to ref_df's 'ref'
    .otherwise(pl.col("effect_allele")) # ELSE, keep original effect_allele
    .alias("effect_allele"),

    # Conditional transformation for 'other_allele'
    pl.when(
        (pl.col("effect_allele") == pl.col("ref")) &
        (pl.col("other_allele") == pl.col("alt"))
    )
    .then(pl.col("alt")) # THEN, set other_allele to ref_df's 'alt'
    .otherwise(pl.col("other_allele")) # ELSE, keep original other_allele
    .alias("other_allele"),

    # Assign 'variant_id' directly from 'id' column from ref_df
    pl.col("id").alias("variant_id")
])

# Drop unnecessary columns
joined_df = joined_df.drop(["base_pair_location", "ref", "alt", "id"])

print(f"Final processed DataFrame preview: {joined_df.head()}")

# Write the Output DataFrame to a .tsv file 
# Ensure the output directory exists
os.makedirs(op_path, exist_ok=True)
final_output_path = os.path.join(op_path, output_filename)

print(f"Writing processed data to: {final_output_path}")
joined_df.write_csv(final_output_path, separator='\t')
print("Process complete.")