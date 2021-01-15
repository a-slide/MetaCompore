# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name = "tombo"

rule_name="tombo_preprocess"
rule tombo_preprocess:
    input:
        fast5_dir=get_fast5,
        fastq=rules.merge_fastq.output.fastq,
        fasta=rules.get_transcriptome.output.fasta
    output:
        fast5_dir=directory(join("results", module_name, rule_name, "{cond}_{rep}"))
    log: join("logs",module_name, rule_name, "{cond}_{rep}.log"),
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/tombo:1.5.1"
    script: f"../scripts/tombo_preprocess.py"
    # shell: """
    #     echo '## multi_to_single_fast5 ##' > {log}
    #     multi_to_single_fast5 -t {threads} -i {input.fast5_dir} -s {output.fast5_dir} &>> {log}
    #     echo '## tombo preprocess annotate_raw_with_fastqs ##' >> {log}
    #     tombo preprocess annotate_raw_with_fastqs --fast5-basedir {output.fast5_dir} --fastq-filenames {input.fastq} --overwrite --processes {threads} &>> {log}
    #     echo '## tombo resquiggle ##' >> {log}
    #     tombo resquiggle --rna --overwrite --processes {threads} {output.fast5_dir} {input.fasta} &>> {log}
    #     """

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
    # shell: """
    #     tombo detect_modifications level_sample_compare \
    #     {params.opt} --processes {threads} \
    #     --fast5-basedirs {input.test_fast5} \
    #     --alternate-fast5-basedirs {input.control_fast5} \
    #     --statistics-file-basename {output.stats_h5[:-12]} \
    #     &> {log}
    #     """

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
