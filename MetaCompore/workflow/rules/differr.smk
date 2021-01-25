# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name = "differr"

rule_name="differr_compare"
rule differr_compare:
    input:
        control_bam=join("results", "alignment", "alignmemt_merge", "control.bam"),
        test_bam=join("results", "alignment", "alignmemt_merge", "test.bam"),
        fasta = rules.get_transcriptome.output.fasta
    output:
        bed=join("results", module_name, rule_name, "differr_results.bed"),
    log: join("logs",module_name, rule_name, "out.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/differr_nanopore_drs:latest"
    shell: "differr -p {threads} {params.opt} -a {input.control_bam} -b {input.test_bam} -r {input.fasta} -o {output.bed} &> {log}"

# rule_name="differr_postprocess"
# rule differr_postprocess:
