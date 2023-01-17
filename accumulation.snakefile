import pandas as pd

#ACC = ["SRR5936131", "SRR5947006", "SRR5935765", "SRR5936197", "SRR5946923", "SRR5946920"]

m = pd.read_csv("inputs/nonpareil_fig2_samples.csv", header = 0)
# filter to paired end libraries
m = m[(m["paired_info"] == "pe")] # filter to PE reads
m = m[(m["data_availability"] == "available")] # filter to only data available
ACC = m['Run'].unique().tolist()

rule all:
    input: 
        expand("outputs/nonpareil/{acc}.npo", acc = ACC),
        expand("outputs/sourmash_sketch/{acc}_{types}_100k.sig", acc = ACC, types = ["fastp", "raw"])

rule download_runs:
    output:
        r1 = temp("inputs/raw/{acc}_pass_1.fastq.gz"),
        r2 = temp("inputs/raw/{acc}_pass_2.fastq.gz")
    conda: 'envs/sratools.yml'
    shell:'''
    fastq-dump --gzip --outdir inputs/raw --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {wildcards.acc}
    '''

rule fastp:
    """
    The nonpareil documentation states, "Nonpareil expects that the sequencing error is always well below 5%, 
    so we suggest using an expected error cutoff of 1% (i.e., Q>20, or 1 error in 100 nucleotides)."
    This rule trims the raw reads to a Q >= 20 using fastp.
    It also sets the minimum read length to the nonpareil k-mer length of 24 (sourmash auto skips these sequences).
    """
    input:
        r1 = "inputs/raw/{acc}_pass_1.fastq.gz",
        r2 = "inputs/raw/{acc}_pass_2.fastq.gz"
    output:
        r1 = temp("outputs/fastp/{acc}_1.fq.gz"),
        r2 = temp("outputs/fastp/{acc}_2.fq.gz"),
        html = "outputs/fastp/{acc}.html",
        json = "outputs/fastp/{acc}.json"
    conda: "envs/fastp.yml"
    shell:'''
    fastp -q 20 --length_required 24 -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -j {output.json} -h {output.html} -R {wildcards.acc}
    '''

rule tmp_gunzip_for_nonpareil:
    input: "outputs/fastp/{acc}_1.fq.gz"
    output: temp("outputs/temp/{acc}_1.fq")
    shell:'''
    gunzip -c {input} > {output}
    '''

rule nonpareil:
    """
    The nonpareil documentation states, "if you have paired-end reads, you should use only one sister read per pair in Nonpareil." 
    """
    input: "outputs/temp/{acc}_1.fq"
    output: "outputs/nonpareil/{acc}.npo"
    params: prefix = lambda wildcards: "outputs/nonpareil/" + wildcards.acc 
    conda: "envs/nonpareil.yml"
    shell:'''
    nonpareil -s {input} -T kmer -f fastq -b {params.prefix}
    '''

rule sourmash_sketch_raw:
    input: 
        r1 = "inputs/raw/{acc}_pass_1.fastq.gz",
        r2 = "inputs/raw/{acc}_pass_2.fastq.gz"
    output: "outputs/sourmash_sketch/{acc}_raw.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=10000,abund --name {wildcards.acc} -o {output} {input.r1} {input.r2}
    '''

rule sourmash_sketch_fastp:
    input: 
        r1 = "outputs/fastp/{acc}_1.fq.gz",
        r2 = "outputs/fastp/{acc}_2.fq.gz",
    output: "outputs/sourmash_sketch/{acc}_fastp.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sketch dna -p k=21,k=31,k=51,scaled=10000,abund --name {wildcards.acc} -o {output} {input.r1} {input.r2}
    '''

rule sourmash_sketch_downsample:
    input: "outputs/sourmash_sketch/{acc}_{types}.sig"
    output: "outputs/sourmash_sketch/{acc}_{types}_100k.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig downsample --scaled 100000 -o {output} {input}
    '''
    
