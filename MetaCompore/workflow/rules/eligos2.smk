# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name = "eligos2"

rule_name="eligos2_fasta_to_bed"
rule eligos2_fasta_to_bed:
    input:
        fasta = rules.get_transcriptome.output.fasta
    output:
        bed = temp(join("results", module_name, rule_name, "transcriptome_reference.bed"))
    log: join("logs",module_name, rule_name, "out.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    shell: "faidx --transform bed {input.fasta} > {output.bed}"

rule_name="eligos2_pair_diff_mod"
rule eligos2_pair_diff_mod:
    input:
        control_bam=join("results", "alignment", "alignmemt_merge", "control.bam"),
        test_bam=join("results", "alignment", "alignmemt_merge", "test.bam"),
        fasta=rules.get_transcriptome.output.fasta,
        bed=rules.eligos2_fasta_to_bed.output.bed
    output:
        res_tsv=join("results", module_name, rule_name, "test_vs_control_on_transcriptome_reference_combine.txt")
    log: join("logs",module_name, rule_name, "out.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/eligos2:2.0.0"
    script: "../scripts/eligos2_pair_diff_mod.py"

rule_name="eligos2_postprocess"
rule eligos2_postprocess:
    input:
        res_tsv=rules.eligos2_pair_diff_mod.output.res_tsv
    output:
        res_tsv=join("results", "final", "{}_results.tsv".format(module_name))
    log: join("logs",module_name, rule_name, "out.log"),
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    script: "../scripts/eligos2_postprocess.py"
