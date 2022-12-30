import pandas as pd


m = pd.read_csv("inputs/PRJNA60717.csv", header = 0)
SAMPLES = m['sample'].unique().tolist()

# parse out all run accessions (some samples have two accessions, some only have one)
def flatten(l):
    return [item for sublist in l for item in sublist]

RUNSTMP = m['run_accessions'].unique().tolist()
RUNS = list()
for run in RUNSTMP:
    tmp = run.split(";")
    RUNS.append(tmp)

RUNS = flatten(RUNS)

# set sourmash runtime parameters
KSIZES = [31]
LINEAGES=['bacteria']

rule all:
    input: expand("outputs/sourmash_taxonomy/{sample}-vs-genbank-2022.03-k{ksize}.with-lineages.csv", sample = SAMPLES, ksize = KSIZES)

rule download_runs:
    output:
        r1 = "inputs/raw/{runs}_pass_1.fastq.gz",
        r2 = "inputs/raw/{runs}_pass_2.fastq.gz"
    conda: 'envs/sratools.yml'
    shell:'''
    fastq-dump --gzip --outdir inputs/raw --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {wildcards.runs}
    '''

rule sketch_samples:
    input:
        m = "inputs/PRJNA60717.csv",
        r1 = expand("inputs/raw/{runs}_pass_1.fastq.gz", runs = RUNS),
        r2 = expand("inputs/raw/{runs}_pass_2.fastq.gz", runs = RUNS)
    output: expand("outputs/sourmash_sketch/{sample}.sig", sample = SAMPLES)
    run:
        shell("mkdir -p outputs/sourmash_sketch")
        for sample in SAMPLES:
            row = m.loc[m['sample'] == sample]
            sig_path = os.path.join("outputs", "sourmash_sketch", sample + ".sig")
            runs = row['run_accessions'].values[0]
            runs = runs.split(";")
            if len(runs) == 1:
                # if there was only one run accession for the sample, sketch the run accession with the sample name as the sketch name
                r1_path = os.path.join("inputs", "raw", runs + "_pass_1.fastq.gz")
                r2_path = os.path.join("inputs", "raw", runs + "_pass_2.fastq.gz")
                shell("sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name {sample} -o {sig_path} {r1_path} {r2_path}")
            else:
                # else there are two accessions and they need to be sketched into a single signature
                runs0 = runs[0]
                runs1 = runs[1]
                runs0_r1_path = os.path.join("inputs", "raw", runs0 + "_pass_1.fastq.gz")
                runs0_r2_path = os.path.join("inputs", "raw", runs0 + "_pass_2.fastq.gz")
                runs1_r1_path = os.path.join("inputs", "raw", runs1 + "_pass_1.fastq.gz")
                runs1_r2_path = os.path.join("inputs", "raw", runs1 + "_pass_2.fastq.gz")
                shell("sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --name {sample} -o {sig_path} {runs0_r1_path} {runs0_r2_path} {runs1_r1_path} {runs1_r2_path}")

#############################################################
## Download and prep sourmash databases
#############################################################

rule download_sourmash_databases_genbank:
    input: "inputs/sourmash_databases/sourmash-database-info.csv"
    output: "inputs/sourmash_databases/genbank-2022.03-{lineage}-k{ksize}.zip"
    run:
        sourmash_database_info = pd.read_csv(str(input[0]))
        ksize = int(wildcards.ksize)
        lineage_df = sourmash_database_info.loc[(sourmash_database_info['lineage'] == wildcards.lineage) & (sourmash_database_info['ksize'] == ksize)]
        if lineage_df is None:
            raise TypeError("'None' value provided for lineage_df. Are you sure the sourmash database info csv was not empty?")

        osf_hash = lineage_df['osf_hash'].values[0] 
        shell("curl -JLo {output} https://osf.io/{osf_hash}/download")


rule download_sourmash_lineages_genbank:
    input: "inputs/sourmash_databases/sourmash-lineage-info.csv"
    output: "inputs/sourmash_databases/genbank-2022.03-{lineage}.lineages.csv.gz"
    run:
        sourmash_lineage_info = pd.read_csv(str(input[0]))
        lineage_df = sourmash_lineage_info.loc[sourmash_lineage_info['lineage'] == wildcards.lineage]
        if lineage_df is None:
            raise TypeError("'None' value provided for lineage_df. Are you sure the sourmash database info csv was not empty?")

        osf_hash = lineage_df['osf_hash'].values[0] 
        shell("curl -JLo {output} https://osf.io/{osf_hash}/download")

rule gunzip_lineage_csvs:
    input: "inputs/sourmash_databases/genbank-2022.03-{lineage}.lineages.csv.gz"
    output: "inputs/sourmash_databases/genbank-2022.03-{lineage}.lineages.csv"
    shell:'''
    gunzip -c {input} > {output}
    '''

rule sourmash_taxonomy_prepare:
    input: expand("inputs/sourmash_databases/genbank-2022.03-{lineage}.lineages.csv", lineage = LINEAGES),
    output: "outputs/sourmash_taxonomy/genbank-2022.03-prepared-lineages.sqldb"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash tax prepare --taxonomy-csv {input} -o {output}
    '''

###############################################################
## Run taxonomic classification
###############################################################

rule sourmash_gather:
    input:
        sig="outputs/sourmash_sketch/{sample}.sig",
        databases=expand("inputs/sourmash_databases/genbank-2022.03-{lineage}-k{{ksize}}.zip", lineage = LINEAGES)
    output: csv="outputs/sourmash_gather/{sample}-vs-genbank-2022.03-k{ksize}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash gather -k {wildcards.ksize} --scaled 1000 --threshold-bp 0 -o {output.csv} {input.sig} {input.databases}
    '''

rule sourmash_taxonomy_annotate:
   input:
       lin_prepared="outputs/sourmash_taxonomy/genbank-2022.03-prepared-lineages.sqldb",
       gather="outputs/sourmash_gather/{sample}-vs-genbank-2022.03-k{ksize}.csv"
   output: "outputs/sourmash_taxonomy/{sample}-vs-genbank-2022.03-k{ksize}.with-lineages.csv"
   params: outdir = "outputs/sourmash_taxonomy/"
   conda: "envs/sourmash.yml"
   shell:'''
   sourmash tax annotate -g {input.gather} -t {input.lin_prepared} -o {params.outdir}
   '''
