# BIOSTAT 666 Final Project: Comparing Polygenic Risk Scores Calculated from GWAS Summary Statistics

**Authors:** Gabrielle Dotson, Mengtong Hu, Jeremy Kaplan, and Soumik Purkayastha

## Setup

### Conda Environment
All of the library requirements for this project can be installed in a conda environment. From a command line you can build the environment as
    
    conda env create -f environment.yml

Once the environment builds, you can activate the environment by running:

    conda activate b666prj

### Data
Reference genotypes were downloaded from the [1000 Genomes Project FTP Server](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) 
Download files to `input/vcf`

GWAS summary statistics were downloaded from the [UNC Psychiatric Genomics Consortium](https://www.med.unc.edu/pgc/). Download and unzip files into `input/gwas`

## Running the workflow
Currently the workflow will process VCFs and generate clumped SNPs on the merged denotypes. You can run the entire workflow by running:

    snakemake output/output/gwas_snps.clumped 
This will create a number of outputs in `output/` containing information about the clumped SNPs.