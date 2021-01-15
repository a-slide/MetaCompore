# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name = "nanocompore"

rule_name = "nanocompore_eventalign_collapse"
rule nanocompore_eventalign_collapse:
    input: tsv = rules.nanopolish_eventalign_gather.output.tsv
    output:
        outdir = directory(join("results", module_name, rule_name, "{cond}_{rep}")),
        tsv = join("results", module_name, rule_name, "{cond}_{rep}", "out_eventalign_collapse.tsv"),
        idx = join("results", module_name, rule_name, "{cond}_{rep}", "out_eventalign_collapse.tsv.idx")
    log: join("logs", module_name, rule_name, "{cond}_{rep}.log")
    threads: get_threads(config, rule_name)
    params: opt = get_opt(config, rule_name)
    resources: mem_mb = get_mem(config, rule_name)
    container: "library://aleg/default/nanocompore:1.0.2"
    shell: "nanocompore eventalign_collapse -t {threads} {params.opt} --overwrite -i {input.tsv} -o {output.outdir} &> {log}"

rule_name = "nanocompore_sampcomp"
rule nanocompore_sampcomp:
    input:
        control_tsv=expand(join("results", module_name, "eventalign_collapse", "control_{rep}", "out_eventalign_collapse.tsv"), rep=replicates_list),
        test_tsv=expand(join("results", module_name, "eventalign_collapse", "test_{rep}", "out_eventalign_collapse.tsv"), rep=replicates_list),
        fasta = rules.get_transcriptome.output.fasta
    output:
        outdir = directory(join("results", module_name, rule_name))
    log: join("logs",module_name, rule_name, "sampcomp.log")
    threads: get_threads(config, rule_name)
    params: opt = get_opt(config, rule_name)
    resources: mem_mb = get_mem(config, rule_name)
    container: "library://aleg/default/nanocompore:1.0.2"
    shell: "nanocompore sampcomp -t {threads} {params.opt} --overwrite --file_list1 {input.control_tsv} --file_list2 {input.test_tsv} --label1 control --label2 test --fasta {input.fasta} --outpath {output.outdir}"

# rule_name = "nanocompore_postprocess"
# rule nanocompore_postprocess:
