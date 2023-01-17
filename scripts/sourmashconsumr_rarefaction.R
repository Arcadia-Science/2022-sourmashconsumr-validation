library(readr)
library(sourmashconsumr)

signature_df <- read_signature(snakemake@input[['sig']]) %>%
  dplyr::filter(ksize == 31)

rarefaction_df <- from_signatures_to_rarefaction_df(signatures_df = signature_df)

write_tsv(rarefaction_df, snakemake@output[['tsv']])
