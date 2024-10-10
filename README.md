# HD_sample_enrichment

This repository contains the analysis and statistical simulation code used for my Master's Paper, "Sample Enrichment Strategies for Prodromal Huntington's Disease Trials:Power and Generalizability".

A PDF of the final version of the paper can be found in this repository [*(Sample Enrichment Strategies for Prodromal Huntingtons Disease.pdf)*]("https://github.com/bbodek/HD_sample_enrichment/blob/main/Sample%20Enrichment%20Strategies%20for%20Prodromal%20Huntingtons%20Disease.pdf")

This code contains scripts which accomplish 2 main tasks:

1) Run a statistical analysis on a dataset from the ENROLL-HD observational study to obtain parameter estimates for use in the simulation. This dataset is not publically available, so this analysis cannot be reproduced using only this repository.

2) Use the parameter estimates obtained from part 1 to run a sequence of Power Simulations. Parameter estimates are saved in "data/parameters_estimates.csv", so the simulation portion of this code can be reproduced using this repository.
