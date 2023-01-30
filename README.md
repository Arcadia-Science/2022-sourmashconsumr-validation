# Validating new functionality in the sourmashconsumr R package

This repository validates new functions in the [sourmashconsumr](https://github.com/Arcadia-Science/sourmashconsumr/) R package.

It currently tests the validity of:
1. The function `from_taxonomy_annotate_to_multi_strains()` to detect the presence of multiple strain of the same species in a single metagenome.
2. The function `from_signatures_to_rarefaction_df()` to assess sequencing saturation.

It also contains notebooks for the figures presented in the sourmashconsumr pub ([DOI: 10.57844/arcadia-1896-ke33](https://arcadia-research.pubpub.org/pub/resource-sourmashconsumr)).

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

The snakefiles `CAMI.snakefile` and `PRJNA60717.snakefile` assess `from_taxonomy_annotate_to_multi_strains()`. 
The snakefile `accumulation.snakefile` assesses `from_signatures_to_rarefaction_df()`.

After running the snakemake workflows, the notebooks in the `notebooks/` directory can be executed using the workflow outputs and the built-in data sets in the sourmashconsumr package.
To run the notebooks, first create and activate the notebook/plotting environment:

```
mamba env create -n pltenv -f pltenv.yml
conda activate pltenv
```

Then, start a jupyter notebook session and run the notebooks

```
jupyter notebook
```
