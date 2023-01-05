ACC = ["SRR5936131", "SRR5947006", "SRR5935765", "SRR5936197", "SRR5946923", "SRR5946920"]

rule download_runs:
    output:
        r1 = "inputs/raw/{acc}_pass_1.fastq.gz",
        r2 = "inputs/raw/{acc}_pass_2.fastq.gz"
    conda: 'envs/sratools.yml'
    shell:'''
    fastq-dump --gzip --outdir inputs/raw --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip {wildcards.acc}
    '''

rule fastp:
    """
    The nonpareil documentation states, "Nonpareil expects that the sequencing error is always well below 5%, 
    so we suggest using an expected error cutoff of 1% (i.e., Q>20, or 1 error in 100 nucleotides)."
    This rule trims the raw reads to a Q >= 20 using fastp
    """
    input:
        r1 = "inputs/raw/{acc}_pass_1.fastq.gz",
        r2 = "inputs/raw/{acc}_pass_2.fastq.gz"
    output:
        r1 = "outputs/fastp/{acc}_1.fq.gz",
        r2 = "outputs/fastp/{acc}_2.fq.gz",
        html = "outputs/fastp/{acc}.html",
        json = "outputs/fastp/{acc}.json"
    conda: "envs/fastp.yml"
    shell:'''
    fastp -q 20 -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} -j {output.json} -h {output.html} -R {wildcards.acc}
    '''

rule nonpareil:
    """
    The nonpareil documentation states, "if you have paired-end reads, you should use only one sister read per pair in Nonpareil." 
    """
    input: "outputs/fastp/{acc}_1.fq.gz"
    output: "outputs/nonpareil/{acc}.npo"
    params: prefix = lambda wildcards: "outputs/nonpareil/" + {wildcards.acc} 
    conda: "envs/nonpareil.html"
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