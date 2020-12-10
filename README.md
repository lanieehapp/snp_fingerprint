# SNP Fingerprinting

This tool generates a SNP fingerprint for a BAM.

## Contents:
* fp.R: R script meant to be run from the command line.
* fingerprint_snps.Rdata: Rdata file that contains a dataframe of the SNP positions identified. Used by fp.R.
* fingerprinting_locations.txt: locations used by mpileup to generate variant calls for 

## Usage:

`Rscript fp.R <BAM file> <Reference FASTA> <Reference FASTA index> <Sample ID> <Rdata output> <Vector output> <Count output>`

