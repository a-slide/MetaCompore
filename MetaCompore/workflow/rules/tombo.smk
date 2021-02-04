# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name="tombo"

rule_name="tombo_preprocess"
rule tombo_preprocess:
    input:
        fast5_dir=get_fast5,
        fastq=rules.merge_fastq.output.fastq,
        fasta=rules.get_transcriptome.output.fasta,
        selected_reads_fn=rules.alignmemt_postfilter.output.selected_reads_fn
    output:
        fast5_dir=directory(join("results", module_name, rule_name, "{cond}_{rep}"))
    log: join("logs",module_name, rule_name, "{cond}_{rep}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/tombo:1.5.1"
    script: "../scripts/tombo_preprocess.py"

rule_name="tombo_level_sample_compare"
rule tombo_level_sample_compare:
    input:
        control_fast5=expand(join("results", module_name, "tombo_preprocess", "control_{rep}"), rep=replicates_list),
        test_fast5=expand(join("results", module_name, "tombo_preprocess", "test_{rep}"), rep=replicates_list)
    output:
        res_h5=join("results", module_name, rule_name, "results.tombo.stats")
    log: join("logs",module_name, rule_name, "out.log")
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
        bn=join("results", module_name, rule_name, "results")
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/tombo:1.5.1"
    shell: "tombo detect_modifications level_sample_compare {params.opt} --processes {threads} --fast5-basedirs {input.test_fast5} --alternate-fast5-basedirs {input.control_fast5} --statistics-file-basename {params.bn} &> {log}"

rule_name="tombo_postprocess"
rule tombo_postprocess:
    input:
        res_h5=rules.tombo_level_sample_compare.output.res_h5,
        fasta=rules.get_transcriptome.output.fasta
    output:
        res_tsv=join("results", "final", "{}_results.tsv".format(module_name))
    log: join("logs", module_name, rule_name, "out.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    script: "../scripts/tombo_postprocess.py"
