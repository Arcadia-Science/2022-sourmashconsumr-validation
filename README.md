# Validating new functionality in the sourmashconsumr R package

This repository validates new functions in the [sourmashconsumr](https://github.com/Arcadia-Science/sourmashconsumr/) R package.
It currently tests the validity of the function `from_taxonomy_annotate_to_multi_strains()` to detect the presence of multiple strain of the same species in a single metagenome.

## Getting started

This repository uses snakemake to run the pipeline and conda to manage software environments and installations.
You can find operating system-specific instructions for installing miniconda [here](https://docs.conda.io/en/latest/miniconda.html).
After installing conda and [mamba](https://mamba.readthedocs.io/en/latest/), run the following command to create the pipeline run environment.

```
mamba env create -n sourmashconsumrpub -f environment.yml
conda activate sourmashconsumrpub
```

After installing the environment, you can execute any of the snakefiles in the repository.
For example, to run the CAMI workflow, run:

```
snakemake -s CAMI.snakefile -j 1 --use-conda
```
