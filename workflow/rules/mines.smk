# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name = "input"

rule_name="rule1"
rule rule1:
    input: fa=transcriptome_fa
    output: fa=join("results", module_name, rule_name, "transcriptome_reference.fa")
    log: join("logs",module_name, rule_name, "name.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/CONTAINER:TAG"
    script: "scripts/SCRIPT.py" 