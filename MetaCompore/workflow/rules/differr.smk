# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name = "differr"

rule_name="differr_compare"
rule differr_compare:
    input:
        control_bam=expand(join("results", "alignment", "alignmemt_postfilter", "control_{rep}.bam"), rep=replicates_list),
        test_bam=expand(join("results", "alignment", "alignmemt_postfilter", "test_{rep}.bam"), rep=replicates_list),
        fasta = rules.get_transcriptome.output.fasta
    output:
        res_bed=join("results", module_name, rule_name, "differr_results.bed"),
    log: join("logs",module_name, rule_name, "out.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/differr_nanopore_drs:latest"
    script: "../scripts/differr_compare.py"

rule_name="differr_postprocess"
rule differr_postprocess:
    input:
        res_bed=rules.differr_compare.output.res_bed
    output:
        res_tsv=join("results", "final", "{}_results.tsv".format(module_name))
    log: join("logs",module_name, rule_name, "out.log"),
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    script: "../scripts/differr_postprocess.py"
