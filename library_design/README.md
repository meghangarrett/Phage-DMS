# Phage-DMS Library Design
This repository contains scripts and sequences for making Phage-DMS libraries for the gp41, gp120, and V3 regions of HIV Env. Scripts were written by Kate H.D. Crawford and Jesse D. Bloom and experiment design - including sequence selection and optimization - was conducted by Meghan Garrett, Hannah Itell, and Julie Overbaugh.

## Contents

This directory contains the scripts and HIV Env sequences necessary for designing Phage-DMS libraries for gp41 and V3 from the BF520, BG505, and ZAA1197 strains and for gp120 from the BG505, B41, and DU422 strains. 
The wildtype sequences for these proteins are in the `sequences` directory. 
Note that the scripts require both the DNA and protein sequences as `.fasta` files.
The oligos are designed based on the DNA sequences, but the protein sequences are used to double check no unintended changes to the protein sequence are introduced when designing libraries.

## Making Phage-DMS libraries

Phage-DMS libraries consist of oligos encoding overlapping peptides spanning the protein of interest where the middle residue for each oligo has been randomized to encode any of the 20 amino acids.
A complete library thereby consists of 20 oligos per site in the protein of interest. 

The scripts contained in this repository create oligo libraries for the gp41, gp120, and V3 regioons of Env from several HIV strains.
The oligos created by these scripts are 93 nucleotides long and each oligo advances by 3 nucleotides compared to the previous oligo. (These parameters are specified by the `oligo_length` and `tile` variables defined in the beginning of each script.)
These scripts also add on the necessary adaptor sequences for cloning these oligos into the [T7 Select phage system](https://www.emdmillipore.com/US/en/product/T7Select10-3-Cloning-Kit,EMD_BIO-70550#anchor_USP).
Additionally, user-specified linker sequences are added to the first and last oligos ordered to ensure the first and last sites to be mutated are in the middle of those oligos. 
The scripts also use synonymous substitution to remove any sequences matching the restriction sites used for cloning. 
For the T7 Select system, those restriction sites are `GAATTC` and `AAGCTT` matching the motifs for the EcoRI and HindIII restriction enzymes, respectively.
The sequences to remove from the designed oligos are referred to as `avoid_motifs` and can be specified in the beginning of the scripts.
Finally, since there is a lot of sequence similarity between the different HIV proteins for which we were designing oligos, the script removes all duplicate sequences, so each oligo ordered is unique.

There are different scripts for each V3 library as each virus in this library had different linker sequences for the first and last oligos. The `gp41_oligos.py` script makes libraries for gp41 from the BG505, BF520, and ZAA1197 strains. The `gp120_oligos_181219.py` script makes Phage-DMS oligo libraries for gp120 from BG505, B41, and DU422. Of note, only the `gp120_oligos_181219.py` script allows for the input wildtype sequences to be different lengths.

The output for each script is a `.csv` file with the columns `Virus`, `Rand_Loc`, `Rand_AA`, and `Oligo`. 
The `Virus` column specifies the virus strain for which the oligos are designed; `Rand_Loc` specifies what residue is in the middle of the oligo and, thus, being randomized; `Rand_AA` specifies what amino acid is encoded in this middle position; and `Oligo` is the final oligo sequence including necessary adaptors.
These `.csv` files are output to an `oligos` directory.

To run these scripts and re-make the oligos we designed, clone this repo and run each script from the command line (for example, use the command: `python gp41_oligos.py` to re-make the gp41 oligos).
These scripts require Python version 3.6 and pandas.
To facilitate running these scripts, we have also included an `envrionment.yml` file.
This can be used to create a conda environment compatible with running these scripts for library design using the command: `conda env create -f environment.yml`.
