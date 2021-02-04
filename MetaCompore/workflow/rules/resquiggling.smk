# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name="resquiggling"

if config["gpu_acceleration"]:
    f5c_container="library://aleg/default/f5c:gpu-0.6"
else:
    f5c_container="library://aleg/default/f5c:cpu-0.6"

rule_name="f5c_index"
rule f5c_index:
    input:
        fast5_dir=get_fast5,
        fastq=rules.merge_fastq.output.fastq
    output:
        index=rules.merge_fastq.output.fastq+".index",
        index_fai=rules.merge_fastq.output.fastq+".index.fai",
        index_gzi=rules.merge_fastq.output.fastq+".index.gzi",
        index_readdb=rules.merge_fastq.output.fastq+".index.readdb"
    log: join("logs", module_name, rule_name, "{cond}_{rep}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: f5c_container
    shell: "f5c index {params.opt} -t {threads} -d {input.fast5_dir} {input.fastq} 2> {log}"

rule_name="f5c_eventalign"
rule f5c_eventalign:
    input:
        fastq=rules.merge_fastq.output.fastq,
        index=rules.f5c_index.output.index,
        bam=rules.alignmemt_postfilter.output.bam,
        fasta=rules.get_transcriptome.output.fasta,
        kmer_model="resources/f5c/r9.4_70bps.u_to_t_rna.5mer.template.model"
    output:
        tsv=join("results", module_name, rule_name, "{cond}_{rep}_data.tsv"),
        summary=join("results", module_name, rule_name, "{cond}_{rep}_summary.tsv")
    log: join("logs", module_name, rule_name, "{cond}_{rep}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: f5c_container
    shell: "f5c eventalign {params.opt} -t {threads} --kmer-model {input.kmer_model} -r {input.fastq} -b {input.bam} -g {input.fasta} --summary {output.summary}  > {output.tsv} 2> {log}"
