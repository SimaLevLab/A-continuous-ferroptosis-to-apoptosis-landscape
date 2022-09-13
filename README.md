# A-continuous-ferroptosis-to-apoptosis-landscape
Here we provide the code used for recreating the computational results for our publication "A continuous ferroptosis-to-apoptosis landscape identifies  ferroptosis biomarkers and repressors for cancer therapy" submitted to "Cancer Discovery".

The code was written in R 3.6, and verified in R 4.2.

Notes:
1. the first script, "1 - RNAseq analysis - master script.R", should be run before any other script. It contains the libraries needed for the code, definition of some helper functions, the code for analyzing the differentially expressed genes in our RNAseq, and the code used to create the "Gradient Gene Signature" (GGS). These provides the data to all other scripts (2 to 8) which analyze this RNAseq and GGS.
