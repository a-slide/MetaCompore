# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name="nanocompore"

rule_name="nanocompore_eventalign_collapse"
rule nanocompore_eventalign_collapse:
    input:
        tsv=rules.f5c_eventalign.output.tsv
    output:
        outdir=directory(join("results", module_name, rule_name, "{cond}_{rep}")),
        tsv=join("results", module_name, rule_name, "{cond}_{rep}", "out_eventalign_collapse.tsv"),
        idx=join("results", module_name, rule_name, "{cond}_{rep}", "out_eventalign_collapse.tsv.idx")
    log: join("logs", module_name, rule_name, "{cond}_{rep}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/nanocompore:1.0.3"
    shell: "nanocompore eventalign_collapse -t {threads} {params.opt} --overwrite -i {input.tsv} -o {output.outdir} &> {log}"

rule_name="nanocompore_sampcomp"
rule nanocompore_sampcomp:
    input:
        control_tsv=expand(join("results", module_name, "nanocompore_eventalign_collapse", "control_{rep}", "out_eventalign_collapse.tsv"), rep=replicates_list),
        test_tsv=expand(join("results", module_name, "nanocompore_eventalign_collapse", "test_{rep}", "out_eventalign_collapse.tsv"), rep=replicates_list),
        fasta=rules.get_transcriptome.output.fasta
    output:
        res_tsv=join("results", module_name, rule_name, "outnanocompore_results.tsv"),
        shift_tsv=join("results", module_name, rule_name, "outnanocompore_shift_stats.tsv"),
        res_db=join("results", module_name, rule_name, "outSampComp.db"),
    log: join("logs",module_name, rule_name, "sampcomp.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/nanocompore:1.0.3"
    script: f"../scripts/nanocompore_sampcomp.py"

rule_name="nanocompore_postprocess"
rule nanocompore_postprocess:
    input:
        res_tsv=rules.nanocompore_sampcomp.output.res_tsv,
        fasta=rules.get_transcriptome.output.fasta
    output:
        join("results", "final", "nanocompore_results_GMM_context_0.tsv"),
        join("results", "final", "nanocompore_results_GMM_context_2.tsv"),
        join("results", "final", "nanocompore_results_KS_dwell_context_0.tsv"),
        join("results", "final", "nanocompore_results_KS_dwell_context_2.tsv"),
        join("results", "final", "nanocompore_results_KS_intensity_context_0.tsv"),
        join("results", "final", "nanocompore_results_KS_intensity_context_2.tsv")
    log: join("logs", module_name, rule_name, "out.log"),
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    script: "../scripts/nanocompore_postprocess.py"
