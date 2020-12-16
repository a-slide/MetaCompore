# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name = "alignment"

rule_name = "minimap2_index"
rule minimap2_index:
    input: fasta = rules.get_transcriptome.output.fasta
    output: idx = join("results", module_name, rule_name, "transcriptome_reference.mmi")
    log: join("logs",module_name, rule_name, "name.log")
    threads: get_threads(config, rule_name)
    params: opt = get_opt(config, rule_name)
    resources: mem_mb = get_mem(config, rule_name)
    container: "library://aleg/default/minimap2:2.17"
    shell: "minimap2 -t {threads} {params.opt} -d {output.idx} {input.fasta} &> {log}"

rule_name = "minimap2_align"
rule minimap2_align:
    input:
        idx = rules.minimap2_index.output.idx,
        fastq = rules.merge_fastq.output.fastq
    output:
        bam = join("results", module_name, rule_name, "{sample}.bam"),
        bam_index = join("results", module_name, rule_name, "{sample}.bam.bai")
    log: join("logs",module_name, rule_name, "{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt = get_opt(config, rule_name)
    resources: mem_mb = get_mem(config, rule_name)
    container: "library://aleg/default/minimap2:2.17"
    shell: "minimap2 -t {threads} {params.opt} {input.idx} {input.fastq} 2> {log} | \
            samtools view -bh 2>> {log} | \
            samtools sort > {output.bam} 2>> {log} \
            && samtools index {output.bam} 2>> {log}"

rule_name = "pbt_alignmemt_filter"
rule pbt_alignmemt_filter:
    input:
        bam = rules.minimap2_align.output.bam,
        bam_index = rules.minimap2_align.output.bam_index
    output:
        bam = join("results", module_name, rule_name, "{sample}.bam"),
        bam_index = join("results", module_name, rule_name, "{sample}.bam.bai")
    log: join("logs",module_name, rule_name, "{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt = get_opt(config, rule_name)
    resources: mem_mb = get_mem(config, rule_name)
    container: "library://aleg/default/pybiotools:0.2.4"
    shell: "pyBioTools Alignment Filter {params.opt} -i {input.bam} -o {output.bam} --verbose &> {log}"
