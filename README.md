# A-continuous-ferroptosis-to-apoptosis-landscape
Here we provide the code used for recreating the computational results for our publication "A continuous ferroptosis-to-apoptosis landscape identifies  ferroptosis biomarkers and repressors for cancer therapy" submitted to "Cancer Discovery".

The code was written in R 3.6, and verified in R 4.2.

Notes:
1. The first script, "1 - RNAseq analysis - master script.R", should be run before any other script. It contains the libraries needed for the code, definition of some helper functions, the code for analyzing the differentially expressed genes in our RNAseq, and the code used to create the "Gradient Gene Signature" (GGS). These provides the data to all other scripts (2 to 8) which analyze this RNAseq and GGS.
2. All the data for these script is supplied in the "data" folder.
3. All figures in the paper should be reproduced exactly as they are using this code, except the UMAP shown in figure 6C and the classification metrics shown in figure 6D (and their supplementary figures S6B and S6C) - those are based on a random seed which was changed from the R version used to create the figures in the paper (V 3.5) and the current one (V. 4.2). The results, however, are very similar.
