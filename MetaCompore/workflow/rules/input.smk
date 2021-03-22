# -*- coding: utf-8 -*-

##### Imports #####

from os.path import join
from snakemake.remote.FTP import RemoteProvider as FTP
from snakemake.remote.HTTP import RemoteProvider as HTTP

##### Define RemoteProvider if needed #####

transcriptome_ref = config["transcriptome_ref"]

if transcriptome_ref.startswith("ftp"):
    transcriptome_ref = FTP().remote(transcriptome_ref)

elif transcriptome_ref.startswith("http"):
    transcriptome_ref = HTTP().remote(transcriptome_ref)


##### Rules #####

module_name = "input"

rule_name="get_transcriptome"
rule get_transcriptome:
    input: fasta=transcriptome_ref
    output:
        fasta = join("results", module_name, rule_name, "transcriptome_reference.fa"),
        fai = join("results", module_name, rule_name, "transcriptome_reference.fa.fai")
    log: join("logs",module_name, rule_name, "transcriptome_reference.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    script: f"../scripts/get_transcriptome.py"

rule_name="generate_transcriptome_picard_index"
rule generate_transcriptome_picard_index:
    input: fasta=rules.get_transcriptome.input.fasta
    output:
        picard_dict = join("results", module_name, "get_transcriptome", "transcriptome_reference.fa.dict"),
    log: join("logs",module_name, rule_name, "transcriptome_reference_picard.log")
    threads: 1
    resources: mem_mb=get_mem(config, rule_name)
    container: "docker://broadinstitute/picard"
    shell: 
        """
        perl -pe "s/'$/primesymbol/" {input.fasta} > {input.fasta}.sanitized
        java -jar /usr/picard/picard.jar CreateSequenceDictionary R={input.fasta}.sanitized O={output.picard_dict}.sanitized
        perl -pe "s/primesymbol/'/" {output.picard_dict}.sanitized > {output.picard_dict}
        """
