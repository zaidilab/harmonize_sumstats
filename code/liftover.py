import polars as pl
from pyliftover import LiftOver
import time

def liftover_batch(df, genome_build):
    """
    Python function to perform the liftover operation on a batch (chunk) of data.
    'df': a Polars DataFrame chunk.
    'genome_build': genome build of the GWAS file
    """
    # Initialize an empty list to store the lifted coordinates.
    lifted = []

    # Initialize Liftover object for the current batch.
    lo = LiftOver(genome_build, "hg38")

    # Iterate over pairs of 'chromosome' and 'base_pair_location' from the input DataFrame chunk.
    for chrom, pos in zip(df["chromosome"], df["base_pair_location"]):
        # Add 'chr' prefix for pyliftover
        result = lo.convert_coordinate(f"chr{chrom}", pos)

        # Append the converted location to the 'lifted' list.
        # 'result[0][1]' extracts the base pair location from the first successful conversion.
        # If 'result' is empty (meaning conversion failed), append None.
        lifted.append(result[0][1] if result else None)

    # Add the newly computed '_hg38' column to the original DataFrame chunk.
    return df.with_columns(pl.Series("_hg38", lifted, dtype=pl.Int64))

def perform_liftover(gwas_df, genome_build):
    # Define the expected schema after liftover_batch, including the new '_hg38' column.
    # Modify according to how your input GWAS file is structured
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

    print(f"Performing liftover from {genome_build} to hg38...")
    start = time.time()
    # Convert the 'gwas_df' into a LazyFrame (operations are not executed immediately but are instead recorded as a plan) 
    # because map_batches() method are for lazyFrames only
    gwas_df = (gwas_df.lazy()
                .map_batches(                                       # Applies the 'liftover_batch' Python function to chunks (batches) of the LazyFrame.
                    lambda df: liftover_batch(df, genome_build),      
                    schema=output_schema_after_liftover)            
                .collect())                                         # Executes the LazyFrame's plan and returns the result as an eager DataFrame                                                                
    print(f"Liftover completed in {time.time() - start:.2f} seconds")
    return gwas_df
