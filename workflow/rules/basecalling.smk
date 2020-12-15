# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name = "basecalling"

rule_name = "ont_guppy"
rule ont_guppy:
    input: fast5_dir = get_fast5
    output:
        seq_sum = join("results", module_name, rule_name, "{sample}","sequencing_summary.txt"),
        fastq_dir = directory(join("results", module_name, rule_name, "{sample}"))
    log: join("logs", module_name, rule_name, "{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt = get_opt(config, rule_name)
    resources: mem_mb = get_mem(config, rule_name)
    container: "library://aleg/default/ont_guppy_gpu:4.2.2" if config["gpu_basecalling"] else "library://aleg/default/ont_guppy_cpu:4.2.2"
    shell: "guppy_basecaller {params.opt} -i {input.fast5_dir} -s {output.fastq_dir} &> {log}"

rule_name = "merge_fastq"
rule merge_fastq:
    input: fastq_dir = rules.ont_guppy.output.fastq_dir
    output: fastq = join("results", module_name, rule_name, "{sample}.fastq")
    log: join("logs",module_name, rule_name, "{sample}.log")
    threads: get_threads(config, rule_name)
    params: opt = get_opt(config, rule_name)
    resources: mem_mb = get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    shell: "pyBioTools Fastq Filter {params.opt} -i {input.fastq_dir} -o {output.fastq} --verbose &> {log}"
