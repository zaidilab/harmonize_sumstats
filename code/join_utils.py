import polars as pl

def join_and_filter(gwas_df, ref_df):
    '''
    Performs the inner join and filtering on the reference and GWAS DataFrames.
    gwas_df: DataFrame having the GWAS Summary Statistics file
    ref_df: DataFrame having the reference data (.snplist)
    '''
    # Perform Inner Join 
    return gwas_df.join(
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
