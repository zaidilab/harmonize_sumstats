import polars as pl
from constants import REF_SNP_LIST_PATH, REF_DF_SCHEMA

def load_ref_data():
    '''
    Loads the reference .snplist data from the file path stored in 'REF_SNP_LIST_PATH'
    '''
    print("Loading reference SNP list...")

    # Read Reference SNP list file as a Polars DataFrame and return it
    return pl.read_csv(
        REF_SNP_LIST_PATH,                                            # Specifies the path to the file.
        separator="\t",                                               # Defines the character used to separate columns in the file (here, it's a tab character).
        has_header=False,                                             # Indicates that the file does NOT have a header row.
        new_columns=["chromosome", "position", "id", "ref", "alt"],   # Renames the columns as the default would be "column_1", "column2" and so on
        schema=REF_DF_SCHEMA                                          # Applies the schema 
    )

def load_gwas_data(gwas_summary_path):
    '''
    Loads the GWAS summary data from the file path specified in 'gwas_summary_path' 
    gwas_summary_path: File path given by the user that refers to the GWAS summary data
    '''
    print(f"Loading GWAS summary statistics")
    # Read GWAS Summary Statistics file as a Polars DataFrame and return it
    return pl.read_csv(
        gwas_summary_path,  # Specifies the path to the file
        separator="\t"      # Defines the column separator as a tab
    )

def write_output(df, out_path):
    '''
    Writes the final, processed DataFrame to a .tsv file called GWAS_new_1M
    df: DataFrame to write to the output file
    out_path: File path to GWAS_new_1M.tsv 
    '''
    print(f"Writing processed data to: {out_path}")
    df.write_csv(out_path, separator="\t")
    print("Process complete.")
