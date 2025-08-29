# Sz-vertebrae-SIA
Code and data for reproducing results and figures of the manuscript "Pronounced ontogenetic shifts of Atlantic smooth hammerhead shark revealed by CNS stable isotopes in vertebrae"

1_Main code.R 
Contains the main functions for the data transformation, calculation of standard ellipse areas, respective overlap plotting.
Within this code, either the RAW isotope data "VertebraeSIA_RAW.csv" and corresponding age readings "Vertebrae_age_reading.csv" can be imported to reproduce the entirety of the calculations, OR pre-processed data "VertebraeSIA_Aged.csv", "SEVb.csv", and "Overlap_summary.csv" can be used to run only selected calculations or plotting without running the more time-consuming functions. "Skinner et al. 2019 sppl._ellipsoidcode_functions_final_EDIT.r" contains necessary functions used in the code.

2_C GAM model selection.R
3_N GAM model selection.R
4_S GAM model selection.R
Contain the code to reproduce the generalized mixed effect models (GAMMs) of the manuscript including model selection, and assessment using the pre-processed "VertebraeSIA_Aged.csv".
