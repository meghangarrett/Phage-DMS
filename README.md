# Phage-DMS
<<<<<<< Updated upstream
analysis of Phage-DMS data to determine antibody epitopes
=======
Analysis of Phage-DMS data to determine antibody epitopes

## Contents

This directory contains the scripts and data associated with the manuscript "Phage-DMS: a comprehensive method for fine mapping of antibody epitopes" from Garrett, et al. 

Code used to generate the sequences used to create each Phage-DMS library is described within the directory `library_design`. 

Code used to analyze the results of Phage-DMS experiments is described below.

## Analyzing Phage-DMS experiments

A Phage-DMS experiment is performed by incubating an antibody of interest with the phage display library (for these experiments, either a gp41/V3 Phage-DMS library or a gp120 Phage-DMS library). Antibody will bind to specific peptides displayed by phage, and these antibody-phage complexes are immunoprecipitated, sequences PCR amplified with barcodes, then pooled and deeply sequenced. All data described here was sequenced by single end 125 bp reads on an Illumnia MiSeq. Resulting sequences are then demultiplexed and aligned to the computationally designed sequences (specified by the `metadata` files in the `data` directory), and "annotated counts" sheets are generated with raw count data.

The `data` folder contains files named as `[experiment_date].annotatedCounts.csv` from several Phage-DMS experiments taking place during 2018 and 2019. These sheets have columns specifying the identity of each unique sequence in the Phage-DMS libraries, in addition to a column for the raw "counts" of each sequence within a sample, which were generated by demultiplexing and aligning the reads from deep sequencing experiments. 

The `data` folder also contains files called `gp41V3_PhageDMS_library_key.csv` and `gp120_PhageDMS_library_key.csv` that provide the appropriate amino acid position numbering for for both the gp41/V3 libraries and gp120 libraries. Typically, HIV Envelope positions are numbered according to their relative HXB2 Env position. These "key" tables map the correct HXB2 numbering relative to the native numbering for each Env strain in the library. The "keys" are generated with either the script `make_gp120_library_key.RMD` or `make_gp120_library_key.RMD`. Both require individual fasta files containing the alignment between the Virus/Protein (ex: BG505 V3) and HXB2. This can be done by the LANL tool [HIVAlign](https://www.hiv.lanl.gov/content/sequence/VIRALIGN/viralign.html), aligning against the HIV Env protein sequence with HXB2 with HMM-align.

**Analysis of Phage-DMS libraries**
The scripts called `gp41V3_library_analysis.Rmd` and `gp120_library_analysis.Rmd` are used to examine the composition, diversity, correlation, and coverage of Phage-DMS libraries used in these experiments. These files read in `[experiment_date].annotatedCounts.csv` files containing deeply sequenced input libraries from the `data` folder, in addition to reading in the appropriate library key file from the `data` folder to facilitate labeling with the corresponding HXB2 site number. The plots and results from analysis with these scripts is deposited to the `results/library_analysis` directory.

**Analysis of results with monoclonal antibodies (mAbs)**
The scripts called `[library]_[mAb]_analysis.Rmd` are used to analyze results of Phage-DMS experiments with various mAbs of interest. The scripts will read in in `[experiment_date].annotatedCounts.csv` files from biological replicate experiments, as well as the appropriate library key file. Each script will examine the correlation between technical and biological replicates, as well as calculate the fold enrichment of WT peptides and the scaled differential selection of each peptide. Heatmaps are made in Prism, so a `.csv` file with a matrix of the scaled differential selection values for each amino acid mutation within a region of interest are exported for convenience. Plots and results from analysis with these scripts is deposited to a unique folder for each library/mAb combination within the `results/mAb_analysis` directory.# Phage-DMS
Analysis of Phage-DMS data to determine antibody epitopes

## Contents

This directory contains the scripts and data associated with the manuscript "Phage-DMS: a comprehensive method for fine mapping of antibody epitopes" from Garrett, et al. 

Code used to generate the sequences used to create each Phage-DMS library is described within the directory `library_design`. 

Code used to analyze the results of Phage-DMS experiments is described below.

## Analyzing Phage-DMS experiments

A Phage-DMS experiment is performed by incubating an antibody of interest with the phage display library (for these experiments, either a gp41/V3 Phage-DMS library or a gp120 Phage-DMS library). Antibody will bind to specific peptides displayed by phage, and these antibody-phage complexes are immunoprecipitated, sequences PCR amplified with barcodes, then pooled and deeply sequenced. All data described here was sequenced by single end 125 bp reads on an Illumnia MiSeq. Resulting sequences are then demultiplexed and aligned to the computationally designed sequences (specified by the `metadata` files in the `data` directory), and "annotated counts" sheets are generated with raw count data.

The `data` folder contains files named as `[experiment_date].annotatedCounts.csv` from several Phage-DMS experiments taking place during 2018 and 2019. These sheets have columns specifying the identity of each unique sequence in the Phage-DMS libraries, in addition to a column for the raw "counts" of each sequence within a sample, which were generated by demultiplexing and aligning the reads from deep sequencing experiments. 

The `data` folder also contains files called `gp41V3_PhageDMS_library_key.csv` and `gp120_PhageDMS_library_key.csv` that provide the appropriate amino acid position numbering for for both the gp41/V3 libraries and gp120 libraries. Typically, HIV Envelope positions are numbered according to their relative HXB2 Env position. These "key" tables map the correct HXB2 numbering relative to the native numbering for each Env strain in the library. The "keys" are generated with either the script `make_gp120_library_key.RMD` or `make_gp120_library_key.RMD`. Both require individual fasta files containing the alignment between the Virus/Protein (ex: BG505 V3) and HXB2. This can be done by the LANL tool [HIVAlign](https://www.hiv.lanl.gov/content/sequence/VIRALIGN/viralign.html), aligning against the HIV Env protein sequence with HXB2 with HMM-align.

**Analysis of Phage-DMS libraries**
The scripts called `gp41V3_library_analysis.Rmd` and `gp120_library_analysis.Rmd` are used to examine the composition, diversity, correlation, and coverage of Phage-DMS libraries used in these experiments. These files read in `[experiment_date].annotatedCounts.csv` files containing deeply sequenced input libraries from the `data` folder, in addition to reading in the appropriate library key file from the `data` folder to facilitate labeling with the corresponding HXB2 site number. The plots and results from analysis with these scripts is deposited to the `results/library_analysis` directory.

**Analysis of results with monoclonal antibodies (mAbs)**
The scripts called `[library]_[mAb]_analysis.Rmd` are used to analyze results of Phage-DMS experiments with various mAbs of interest. The scripts will read in in `[experiment_date].annotatedCounts.csv` files from biological replicate experiments, as well as the appropriate library key file. Each script will examine the correlation between technical and biological replicates, as well as calculate the fold enrichment of WT peptides and the scaled differential selection of each peptide. Heatmaps are made in Prism, so a `.csv` file with a matrix of the scaled differential selection values for each amino acid mutation within a region of interest are exported for convenience. Plots and results from analysis with these scripts is deposited to a unique folder for each library/mAb combination within the `results/mAb_analysis` directory.
>>>>>>> Stashed changes
