# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name = "resquiggling"

# if config["gpu_acceleration"]:
#     f5c_container = "library://aleg/default/f5c:gpu-0.5"
# else:
#     f5c_container = "library://aleg/default/f5c:cpu-0.5"
#
# rule_name="f5c_index"
# rule f5c_index:
#     input:
#         fast5_dir = get_fast5,
#         fastq = rules.merge_fastq.output.fastq
#     output:
#         index = rules.merge_fastq.output.fastq+".index",
#         index_fai = rules.merge_fastq.output.fastq+".index.fai",
#         index_gzi = rules.merge_fastq.output.fastq+".index.gzi",
#         index_readdb = rules.merge_fastq.output.fastq+".index.readdb"
#     log: join("logs", module_name, rule_name, "{sample}.log")
#     threads: get_threads(config, rule_name)
#     params: opt=get_opt(config, rule_name)
#     resources: mem_mb=get_mem(config, rule_name)
#     container: f5c_container
#     shell: "f5c index {params.opt} -t {threads} -d {input.fast5_dir} {input.fastq} 2> {log}"
#
# rule_name="f5c_eventalign"
# rule f5c_eventalign:
#     input:
#         fastq = rules.merge_fastq.output.fastq,
#         index = rules.f5c_index.output.index,
#         bam = rules.pbt_alignmemt_filter.output.bam,
#         fasta = rules.get_transcriptome.output.fasta
#     output: tsv=join("results", module_name, rule_name, "{sample}.tsv")
#     log: join("logs", module_name, rule_name, "{sample}.log")
#     threads: get_threads(config, rule_name)
#     params: opt=get_opt(config, rule_name)
#     resources: mem_mb=get_mem(config, rule_name)
#     container: f5c_container
#     shell: "f5c eventalign {params.opt} -t {threads} -r {input.fastq} -b {input.bam} -g {input.fasta} > {output.tsv} 2> {log}"

rule_name="nanopolish_index"
rule nanopolish_index:
    input:
        fast5_dir = get_fast5,
        fastq = rules.merge_fastq.output.fastq,
        seqsum = rules.ont_guppy.output.seqsum
    output:
        index = rules.merge_fastq.output.fastq+".index",
        index_fai = rules.merge_fastq.output.fastq+".index.fai",
        index_gzi = rules.merge_fastq.output.fastq+".index.gzi",
        index_readdb = rules.merge_fastq.output.fastq+".index.readdb"
    log: join("logs", module_name, rule_name, "{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/nanopolish:0.13.2"
    shell: "nanopolish index {params.opt} -v -d {input.fast5_dir} -s {input.seqsum} {input.fastq} 2> {log}"

# Chunked processing to speed up nanopolish
try:
    nchunk=int(config["pbt_alignment_split"]["n_chunks"])
except:
    nchunk=4
chunk_list=list(range(nchunk))

rule_name="pbt_alignment_split"
rule pbt_alignment_split:
    input: bam = rules.pbt_alignmemt_filter.output.bam
    output:
        bam=temp(expand(join("results", module_name, rule_name, "{{sample}}", "{chunk}.bam"), chunk=chunk_list)),
        bam_index=temp(expand(join("results", module_name, rule_name, "{{sample}}", "{chunk}.bam.bai"), chunk=chunk_list))
    log: join("logs", module_name, rule_name, "{sample}.log")
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/pybiotools:0.2.4"
    shell: "pyBioTools Alignment Split {params.opt} -i {input.bam} -l {output.bam} --verbose &>> {log}"

rule_name="nanopolish_eventalign"
rule nanopolish_eventalign:
    input:
        fastq = rules.merge_fastq.output.fastq,
        index = rules.nanopolish_index.output.index,
        fasta = rules.get_transcriptome.output.fasta,
        bam = join("results", module_name, "pbt_alignment_split", "{sample}", "{chunk}.bam"),
        bam_index = join("results", module_name, "pbt_alignment_split", "{sample}", "{chunk}.bam.bai")
    output: tsv=temp(join("results", module_name, rule_name, "{sample}","{chunk}.tsv"))
    log: join("logs", module_name, rule_name, "{sample}", "{chunk}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/nanopolish:0.13.2"
    shell: "nanopolish eventalign {params.opt} -v -t {threads} -r {input.fastq} -b {input.bam} -g {input.fasta} > {output.tsv} 2> {log}"

rule_name="nanopolish_eventalign_gather"
rule nanopolish_eventalign_gather:
    input: tsv_list=expand(join("results", module_name, "nanopolish_eventalign","{{sample}}", "{chunk}.tsv"), chunk=chunk_list)
    output: tsv=join("results", module_name, rule_name, "{sample}.tsv.gz")
    log: join("logs", module_name, rule_name, "{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    script: "../scripts/nanopolish_eventalign_gather.py"
