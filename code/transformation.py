import polars as pl

def apply_transformations(df):
    '''
    Apply transformations (frequency, beta, allele swapping) only if the alleles need flipping.
    df: DataFrame to perform the transformations on
    '''
    return df.with_columns([
        # Changing 'effect_allele_frequency' based on the specified condition
        pl.when(
        (pl.col("effect_allele") == pl.col("ref")) & # If GWAS effect_allele matches ref_df's ref allele
        (pl.col("other_allele") == pl.col("alt"))    # AND GWAS other allele matches ref_df's alt allele
    )
    .then(1 - pl.col("effect_allele_frequency"))  # THEN, change the frequency
    .otherwise(pl.col("effect_allele_frequency")) # ELSE, keep original frequency
    .alias("effect_allele_frequency"),

    # Changing 'beta' based on the specified condition
    pl.when(
        (pl.col("effect_allele") == pl.col("ref")) &
        (pl.col("other_allele") == pl.col("alt"))
    )
    .then(-pl.col("beta"))      # THEN, flip the sign of beta
    .otherwise(pl.col("beta"))  # ELSE, keep original beta
    .alias("beta"),

    # Changing 'effect_allele' based on the specified condition
    pl.when(
        (pl.col("effect_allele") == pl.col("ref")) &
        (pl.col("other_allele") == pl.col("alt"))
    )
    .then(pl.col("alt"))                # THEN, set effect_allele to ref_df's alt
    .otherwise(pl.col("effect_allele")) # ELSE, keep original effect_allele
    .alias("effect_allele"),

    # Changing 'other_allele' based on the specified condition
    pl.when(
        (pl.col("effect_allele") == pl.col("ref")) &
        (pl.col("other_allele") == pl.col("alt"))
    )
    .then(pl.col("ref"))                # THEN, set other_allele to ref_df's ref
    .otherwise(pl.col("other_allele"))  # ELSE, keep original other_allele
    .alias("other_allele"),

    # Assign 'variant_id' directly from 'id' column from ref_df
    pl.col("id").alias("variant_id")
    ])