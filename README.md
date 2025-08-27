# Standardizing GWAS Summary Statistics

## Description
This project provides a Python-based pipeline to standardize GWAS (Genome-Wide Association Study) summary statistics using a reference SNP dataset aligned to the human genome build hg38. It:
- Changes the position of genetic variants if the GWAS data is in hg19, so it matches the hg38 version using a tool called [PyLiftover](https://github.com/konstantint/pyliftover)
- Makes sure that the effect allele in the GWAS file match the alternate allele in the reference file
- Adjusts beta (effect size) and frequency if the alleles are flipped
- Saves the new, updated GWAS file as a .tsv named `GWAS_new_1M.tsv`

The containerization is done using Conda for reproducibility and can be run from the command line with user-defined input paths and genome build.

## Project File Descriptions

The code is modularized in the `code/` directory as follows:

| File               | Description |
|--------------------|-------------|
| `main.py`          | Entry point of the pipeline. Parses CLI arguments, performs data loading, liftover (if needed), joining, allele transformation, and output writing. |
| `liftover.py`      | Contains functions for performing coordinate liftover from other genome builds (e.g., `hg19`) to `hg38` using `pyliftover`. |
| `io_utils.py`      | Handles file I/O operations. Provides functions to load the reference SNP list and GWAS summary statistics, and to write the processed results to disk. |
| `join_utils.py`    | Joins GWAS data with the reference data and filters for matching or flipped allele orientations. |
| `transformation.py`| Applies conditional transformations to modify effect sizes, allele frequencies, and allele labels. |
| `constants.py`     | Stores shared constants like file paths, column schemas, and default output filenames. |

## Getting Started

### 1. Clone the Repository
Clone the GitHub repository and navigate to the directory where the Python script is present.

```
git clone https://github.com/zaidilab/gwas-practice-project.git
cd gwas-practice-project/code
```

### 2. Create and activate conda environment

These commands creates a new Conda environment named `gwas_env` using the dependencies listed in environment.yml and activates the `gwas_env`, so you're ready to run the project.

Make sure to have [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed before executing the following commands. 

```
conda env create --file environment.yml
conda activate gwas_env
```

### 3. Run the Python script

```
python main.py \
  --genome-build hg19 \
  --ip-path "path/to/your/gwas_summary.tsv.gz" \
  --op-path "path/to/save/output/"
```

#### Arguments

- `--genome-build` : Genome build of the input GWAS dataset (hg19 or hg38)
- `--ip-path` : Path to the input GWAS summary statistics file
- `--op-path` : Directory path where the processed and updated GWAS output file should be saved

### Notes
- The input GWAS file should atleast have these columns, named in this format: `chromosome`, `base_pair_location`, `effect_allele`, `other_allele`, `beta`, `effect_allele_frequency`, `variant_id`.
- Ensure your input file is tab-delimited (.tsv) and follows a consistent format.
- The script expects a reference SNP list `1kg_hg38_hm3.snplist` located in a `data/` directory one level above the script. This file is included in the repository. 
- The processed GWAS data will be saved as `GWAS_new_1M.tsv` within the directory specified by `--op-path`.
