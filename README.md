# PlasmoGenomics: A Command-Line Toolkit for Plasmodium Genomic Analysis

**PlasmoGenomics** is a command-line toolkit designed for processing and analyzing *Plasmodium* sequencing data. It includes features for identifying genetic variants, exploring population genetics, and checking for potential drug resistance.

## Key Features

- **Variant Calling**: Processes sequencing data to identify genetic variants using freebayes.
- **Drug Resistance Profiling**: Investigates mutations in key genes associated with antimalarial drug resistance, including:
  - `crt`
  - `mdr1`
  - `dhfr`
  - `dhps`
  - `kelch13`
- **Mutational Burden Analysis**: Assesses the accumulation and distribution of mutations across the genome.
- **Population Genetics Analyses**:
  - **PCA**: Principal Component Analysis for visualizing population structure.
  - **Admixture**: Analysis of genetic ancestry proportions in populations.
  - **F<sub>ST</sub> (Fixation Index)**: Measures genetic differentiation between populations.
  - **IBD (Identity by Descent)**: Detects genomic regions shared by descent.
  - **LD (Linkage Disequilibrium)**: Examines non-random associations of alleles.
  - **Tajima's D**: Tests for selection and demographic events.
  - **MAF (Minor Allele Frequency)**: Profiles the distribution of low-frequency alleles.

## Installation

### Prerequisites

- Python 3.10 or 3.11 (as the project is not compatible with Python versions outside the range `>=3.10,<3.11`).
- Poetry for dependency management.

### Setting up the Environment

1. **Clone the repository**:
    ```bash
    git clone https://github.com/halimasalam/PlasmoGenomics.git
    cd PlasmoGenomics
    ```

2. **Install dependencies with Poetry**:
   - If you don't have Poetry installed, you can install it by following the instructions [here](https://python-poetry.org/docs/#installation).
   - Once Poetry is installed, run the following command to install project dependencies:
     ```bash
     poetry install
     ```
   This will create a virtual environment and install the necessary dependencies.

3. **Activate the virtual environment**:
   - To activate the Poetry-managed virtual environment, use:
     ```bash
     poetry shell
     ```

### Python Version Setup

Ensure that your Python version is within the supported range. If you need to switch to Python `3.10` or `3.11`, use the following command after installing the compatible version of Python:

```bash
poetry env use python3.10  # Or python3.11
```

You can verify that Poetry is using the correct Python version with:

```bash
poetry run python --version
```

This should display a Python version in the range >=3.10,<3.11.


## Dependencies

The following Python packages and command-line tools are required:

### Python Packages:
- pandas
- numpy
- scipy
- matplotlib
- seaborn
- biopython
- pyvcf
- scikit-allel
- plinkio

### Command-Line Tools:
- vcftools
- bcftools
- samtools
- plink
- GATK

Poetry will handle these dependencies automatically when you run poetry install.

## Usage

### Preprocessing Workflow

After setting up your environment and installing dependencies, you can begin preprocessing your genomic data using the provided scripts. Here's an example workflow:

1. **Prepare Your Data**: Ensure your raw sequencing data is organized appropriately, typically in FASTQ format.

2. **Run Preprocessing Script**: Execute the preprocessing script to clean and format your data for analysis. Replace `your_data.fastq` with your actual data file.

```bash
python scripts/preprocess_data.py --input your_data.fastq --output cleaned_data.fastq
```
This script performs quality control, trimming, and formatting to prepare your data for downstream analyses.

3. **Proceed to Analysis**: Once preprocessing is complete, you can use the cleaned data for various analyses such as PCA, drug resistance profiling, and population genetics studies.

## Customizing Preprocessing

The preprocessing script accepts various parameters to customize its behavior. Use the --help flag to see all available options:

```bash
python scripts/preprocess_data.py --help
```

This will display information on how to adjust settings like quality thresholds, adapter trimming, and output formats to suit your specific needs.

Feel free to modify the script or add new ones to accommodate different preprocessing requirements or data types.