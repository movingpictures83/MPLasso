# MPLasso
# Language: R
# Input: TXT
# Output: CSV
# Tested with: PluMA 1.1, R 4.0.0
# Dependency: MASS_7.3.51.6, Matrix_1.2.18, matrixcalc_1.0.3, glasso_1.11, ROCR_1.0.11, huge_1.3.4.1, MCMCpack_1.4.8, glmnet_4.0.2, parallel_4.0.0, matrixStats_0.56.0

PluMA plugin to compute association networks using the MPLasso algorithm  (Lo and Marculescu, 2017).

The plugin takes as input a TXT file of keyword-value pairs, tab-delimited:
interactionFile: initial interactions
occurrenceFile: co-occurrence prior information
countFile: abundances of each taxon, tab delimited (rows are taxa, columns are samples)

It then computes new association networks using graph learning.  These are output as a CSV file,
an NxN matrix where entry (i, j) represents the association between taxon i and taxon j.
