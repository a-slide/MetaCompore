# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name="mines"

rule_name="tombo_de_novo"
rule tombo_de_novo:
    input:
        fast5=expand(join("results", "tombo", "tombo_preprocess", "test_{rep}"), rep=replicates_list),
    output:
        stats_h5=join("results", module_name, rule_name, "test.tombo.stats")
    log: join("logs",module_name, rule_name, "test.log")
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
        bn=join("results", module_name, rule_name, "test")
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/tombo:1.5.1"
    shell: "tombo detect_modifications de_novo {params.opt} --processes {threads} --fast5-basedirs {input.fast5} --statistics-file-basename {params.bn} &> {log}"

rule_name="tombo_de_novo_text_output"
rule tombo_de_novo_text_output:
    input:
        fast5=expand(join("results", "tombo", "tombo_preprocess", "test_{rep}"), rep=replicates_list),
        stats_h5=rules.tombo_de_novo.output.stats_h5
    output:
        coverage=join("results", module_name, rule_name, "test.coverage.plus.bedgraph"),
        fraction=temp(join("results", module_name, rule_name, "test.fraction_modified_reads.plus.wig")),
        coverage_minus=temp(join("results", module_name, rule_name, "test.coverage.minus.bedgraph")),
        fraction_minus=temp(join("results", module_name, rule_name, "test.fraction_modified_reads.minus.wig"))
    log: join("logs",module_name, rule_name, "test.log"),
    threads: get_threads(config, rule_name)
    params:
        opt=get_opt(config, rule_name),
        bn=join("results", module_name, rule_name, "test")
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/tombo:1.5.1"
    shell: "tombo text_output browser_files --fast5-basedirs {input.fast5} --statistics-filename {input.stats_h5} --browser-file-basename {params.bn} --file-types coverage fraction &> {log}"

rule_name="mines_wig2bed"
rule mines_wig2bed:
    input:
        fraction=rules.tombo_de_novo_text_output.output.fraction
    output:
        fraction=join("results", module_name, rule_name, "test.fraction_modified_reads.plus.bed")
    log: join("logs",module_name, rule_name, "out.log"),
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/mines:latest"
    shell: "wig2bed < {input.fraction} > {output.fraction} 2> {log}"

rule_name="mines_cdna"
rule mines_cdna:
    input:
        coverage=rules.tombo_de_novo_text_output.output.coverage,
        fraction=rules.mines_wig2bed.output.fraction,
        fasta=rules.get_transcriptome.output.fasta,
        kmer_models="resources/mines/names.txt"
    output:
        res_bed=join("results", module_name, rule_name, "mines_results.bed")
    log: join("logs",module_name, rule_name, "out.log"),
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/mines:latest"
    shell: "cDNA_MINES --kmer_models {input.kmer_models} --fraction_modified {input.fraction} --coverage {input.coverage} --ref {input.fasta} --output {output.res_bed} &> {log}"

rule_name="mines_postprocess"
rule mines_postprocess:
    input:
        res_bed=rules.mines_cdna.output.res_bed
    output:
        res_tsv=join("results", "final", "{}_results.tsv".format(module_name))
    log: join("logs",module_name, rule_name, "out.log"),
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    script: "../scripts/mines_postprocess.py"
