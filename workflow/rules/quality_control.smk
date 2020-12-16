# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name = "quality_control"

rule_name = "pycoQC"
rule pycoQC:
    input:
        seqsum = rules.ont_guppy.output.seqsum,
        bam = rules.pbt_alignmemt_filter.output.bam
    output:
        json = join("results", module_name, rule_name, "pycoQC_{sample}.json"),
        html = join("results", module_name, rule_name, "pycoQC_{sample}.html")
    log: join("logs", module_name, rule_name, "{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt = get_opt(config, rule_name)
    resources: mem_mb = get_mem(config, rule_name)
    container: "library://aleg/default/pycoqc:2.5.2"
    shell: "pycoQC {params.opt} -f {input.seqsum} -a {input.bam} -o {output.html} -j {output.json} &> {log}"
