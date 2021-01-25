# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name="tombo"

#################################################################################### Add first step to select reads based on alignment sampling results

rule_name="tombo_preprocess"
rule tombo_preprocess:
    input:
        fast5_dir=get_fast5,
        fastq=rules.merge_fastq.output.fastq,
        fasta=rules.get_transcriptome.output.fasta,
        selected_reads_fn=rules.alignmemt_postfilter.output.selected_reads_fn
    output:
        fast5_dir=temp(directory(join("results", module_name, rule_name, "{cond}_{rep}")))
    log: join("logs",module_name, rule_name, "{cond}_{rep}.log"),
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/tombo:1.5.1"
    script: f"../scripts/tombo_preprocess.py"

rule_name="tombo_detect_modifications"
rule tombo_detect_modifications:
    input:
        control_fast5=expand(join("results", module_name, "tombo_preprocess", "control_{rep}"), rep=replicates_list),
        test_fast5=expand(join("results", module_name, "tombo_preprocess", "test_{rep}"), rep=replicates_list)
    output:
        stats_h5=join("results", module_name, rule_name, "results.tombo.stats")
    log: join("logs",module_name, rule_name, "out.log"),
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/tombo:1.5.1"
    script: f"../scripts/tombo_detect_modifications.py"

rule_name="tombo_postprocess"
rule tombo_postprocess:
    input:
        stats_h5=rules.tombo_detect_modifications.output.stats_h5,
        fasta=rules.get_transcriptome.output.fasta
    output:
        stats_tsv=join("results", module_name, rule_name, "tombo_results.tsv"),
    log: join("logs",module_name, rule_name, "out.log"),
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    script: f"../scripts/tombo_postprocess.py"
