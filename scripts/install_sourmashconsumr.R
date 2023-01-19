library(devtools)
install_github("Arcadia-Science/sourmashconsumr")
library(sourmashconsumr)
file.create(snakemake@output[['install']])
