# -*- coding: utf-8 -*-

##### Imports #####

from os.path import join
from snakemake.remote.FTP import RemoteProvider as FTP
from snakemake.remote.HTTP import RemoteProvider as HTTP

##### Define RemoteProvider if needed #####

transcriptome_fa = config["transcriptome_ref"]

if transcriptome_fa.startswith("ftp"):
    transcriptome_fa = FTP().remote(transcriptome_fa)

elif transcriptome_fa.startswith("http"):
    transcriptome_fa = HTTP().remote(transcriptome_fa)


##### Rules #####

module_name = "input"

rule_name="get_transcriptome"
rule get_transcriptome:
    input: fa=transcriptome_fa
    output:
        fa = join("results", module_name, rule_name, "transcriptome_reference.fa"),
        fai = join("results", module_name, rule_name, "transcriptome_reference.fa.fai")
    log: join("logs",module_name, rule_name, "transcriptome_reference.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    script: f"../scripts/{rule_name}.py"
