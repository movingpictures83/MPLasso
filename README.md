# Language: R
# Input: TXT
# Output: CSV
# Tested with: PluMA 1.0, R 4.0

PluMA plugin to compute association networks using the MPLasso algorithm  (Lo and Marculescu, 2017).

The plugin takes as input a TXT file of keyword-value pairs, tab-delimited:
interactionFile: initial interactions
occurrenceFile: co-occurrence prior information
countFile: abundances of each taxon, tab delimited (rows are taxa, columns are samples)

It then computes new association networks using graph learning.  These are output as a CSV file,
an NxN matrix where entry (i, j) represents the association between taxon i and taxon j.
