import pandas as pd

taxonomy = pd.read_csv("inputs/CAMI_low_taxonomy.csv", header = 0)
SOURCE_GENOMES = taxonomy['ident'].unique().tolist()

rule all:
    input: "outputs/sourmash_taxonomy/CAMI_low_vs_source_genomes.with-lineages.csv"

rule download_CAMI:
    output: "inputs/CAMI_low.tar"
    conda: "envs/wget.yml"
    shell:'''
    wget -O {output} https://ftp.cngb.org/pub/gigadb/pub/10.5524/100001_101000/100344/ChallengeDatasets.dir/CAMI_low.tar
    '''

rule decompress_CAMI:
    input: "inputs/CAMI_low.tar"
    output:
        "inputs/CAMI_low/RL_S001__insert_270.fq.gz",
        "inputs/CAMI_low/source_genomes_low.tar.gz",
        "inputs/CAMI_low/gsa_mapping.binning",
        "inputs/CAMI_low/gs_read_mapping.binning.gz",
    shell:'''
    tar xvf {input} -C inputs/
    '''

rule decompress_source_genomes:
    input: "inputs/CAMI_low/source_genomes_low.tar.gz"
    output: expand("inputs/CAMI_low/source_genomes/{source_genome}.fna", source_genome = SOURCE_GENOMES)
    shell:'''
    tar xvf {input} -C inputs/CAMI_low
    # these next few lines are not great, but the source genomes
    # 1) have two separate file endings (fasta, fna);
    # 2) file prefixes do not match the genome identifiers used by the rest of the documents in CAMI low.
    # To fix these problems, the following lines
    # 1) give all source genomes the same file ending (.fna)
    # 2) unify the source genome file prefixes with those used in the rest of the CAMI low documents.
    
    cd inputs/CAMI_low/source_genomes/
    for infile in *fasta 
    do
      bn=$(basename $infile .fasta)
      mv $infile ${{bn}}.fna
    done
    
    for infile in *.gt1kb.fna
    do
      bn=$(basename $infile .gt1kb.fna)
      mv $infile ${{bn}}.fna
    done

    mv 1139_AG_run158_run162.final.scaffolds.fna 1139_AG.fna
    mv 1220_AD_run172_run176.final.scaffolds.fna 1220_AD.fna
    mv 1220_AJ_run172_run176_run188.final.scaffolds.fna 1220_AJ.fna
    mv 1285_BH_run189.final.scaffolds.fna 1285_BH.fna
    mv 1286_AP_run191_run197.final.scaffolds.fna 1286_AP.fna
    mv 1365_A_run201.final.scaffolds.fna 1365_A.fna
    
    cd ../..
    '''

rule sketch_CAMI_source_genomes:
    input: "inputs/CAMI_low/source_genomes/{source_genome}.fna"
    output: "outputs/sourmash_sketch/{source_genome}.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name {wildcards.source_genome} -o {output} {input}
    '''

rule sig_combine_source_genomes:
    input: expand("outputs/sourmash_sketch/{source_genome}.sig", source_genome = SOURCE_GENOMES)
    output: "outputs/sourmash_database/CAMI_low_source_genomes.zip"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig cat -o {output} {input}
    '''

rule sketch_CAMI_mgx:
    input: "inputs/CAMI_low/RL_S001__insert_270.fq.gz"
    output: "outputs/sourmash_sketch/CAMI_low.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name CAMI_low -o {output} {input}
    '''

rule gather_source_genomes_CAMI:
    input:
        sig = "outputs/sourmash_sketch/CAMI_low.sig",
        db = "outputs/sourmash_database/CAMI_low_source_genomes.zip"
    output: csv = "outputs/sourmash_gather/CAMI_low_vs_source_genomes.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -k 31 --scaled 1000 --threshold-bp 0 -o {output.csv} {input.sig} {input.db}
    '''

rule prepare_taxonomy_source_genomes_CAMI:
    input: "inputs/CAMI_low_taxonomy.csv"
    output: "outputs/sourmash_taxonomy/CAMI_low_prepared_lineages.sqldb"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash tax prepare --taxonomy-csv {input} -o {output}
    '''

rule taxonomy_annotate_source_genomes_CAMI:
    input:
       lin_prepared="outputs/sourmash_taxonomy/CAMI_low_prepared_lineages.sqldb",
       gather="outputs/sourmash_gather/CAMI_low_vs_source_genomes.csv"
    output: "outputs/sourmash_taxonomy/CAMI_low_vs_source_genomes.with-lineages.csv"
    params: outdir = "outputs/sourmash_taxonomy/"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash tax annotate -g {input.gather} -t {input.lin_prepared} -o {params.outdir}
    '''
