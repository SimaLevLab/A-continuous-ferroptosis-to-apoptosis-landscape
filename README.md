# A-continuous-ferroptosis-to-apoptosis-landscape
Here we provide the code used for recreating the computational results for our publication "Programming a Ferroptosis-to-Apoptosis Transition Landscape Revealed Ferroptosis Biomarkers and Repressors for Cancer Therapy" published in Advanced Science (DOI: 10.1002/advs.202307263)

The code was written in R 4.2.1

Notes:
1. The first script, "1 - RNAseq analysis - master script.R", should be run before any other script. It contains the libraries needed for the entire code, definition of some helper functions, the code for analyzing the differentially expressed genes in our RNAseq, and the code used to create the "Gradient Gene Signature" (GGS). This script provides the data and libraries to all other scripts which analyze the RNAseq and GGS, and therefore should run first before running any of the other scripts.
2. All the data used for these scripts are supplied in the "data" folder.
3. The figures created by the codes will be saved in the "results" folder.
4. The apoptosis and ferroptosis datasets used in scripts #6 and #7 are in the "data/Datasets" folder.
