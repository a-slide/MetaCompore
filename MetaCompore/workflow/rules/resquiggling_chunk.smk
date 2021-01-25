# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name="resquiggling"

rule_name="nanopolish_index"
rule nanopolish_index:
    input:
        fast5_dir=get_fast5,
        fastq=rules.merge_fastq.output.fastq,
        seqsum=rules.ont_guppy.output.seqsum
    output:
        index=rules.merge_fastq.output.fastq+".index",
        index_fai=rules.merge_fastq.output.fastq+".index.fai",
        index_gzi=rules.merge_fastq.output.fastq+".index.gzi",
        index_readdb=rules.merge_fastq.output.fastq+".index.readdb"
    log: join("logs", module_name, rule_name, "{cond}_{rep}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/nanopolish:0.13.2"
    shell: "nanopolish index {params.opt} -v -d {input.fast5_dir} -s {input.seqsum} {input.fastq} 2> {log}"

# Chunked processing to speed up nanopolish
try:
    nchunk=int(config["alignment_split"]["n_chunks"])
except:
    nchunk=4
chunk_list=list(range(nchunk))

rule_name="alignment_split"
rule alignment_split:
    input:
        bam=rules.alignmemt_postfilter.output.bam,
        bam_index=rules.alignmemt_postfilter.output.bam_index
    output:
        bam=temp(expand(join("results", module_name, rule_name, "{{cond}}_{{rep}}", "{chunk}.bam"), chunk=chunk_list)),
        bam_index=temp(expand(join("results", module_name, rule_name, "{{cond}}_{{rep}}", "{chunk}.bam.bai"), chunk=chunk_list))
    log: join("logs", module_name, rule_name, "{cond}_{rep}.log")
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/pybiotools:0.2.4"
    shell: "pyBioTools Alignment Split {params.opt} -i {input.bam} -l {output.bam} --verbose &>> {log}"

rule_name="nanopolish_eventalign"
rule nanopolish_eventalign:
    input:
        fastq=rules.merge_fastq.output.fastq,
        index=rules.nanopolish_index.output.index,
        fasta=rules.get_transcriptome.output.fasta,
        bam=join("results", module_name, "alignment_split", "{cond}_{rep}", "{chunk}.bam"),
        bam_index=join("results", module_name, "alignment_split", "{cond}_{rep}", "{chunk}.bam.bai")
    output:
        tsv=temp(join("results", module_name, rule_name, "{cond}_{rep}","{chunk}_data.tsv")),
        summary=temp(join("results", module_name, rule_name, "{cond}_{rep}","{chunk}_summary.tsv")),
    log: join("logs", module_name, rule_name, "{cond}_{rep}", "{chunk}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/nanopolish:0.13.2"
    shell: "nanopolish eventalign {params.opt} -v -t {threads} -r {input.fastq} -b {input.bam} -g {input.fasta} --summary {output.summary} > {output.tsv} 2> {log}"

rule_name="eventalign_merge"
rule eventalign_merge:
    input:
        tsv_list=expand(join("results", module_name, "nanopolish_eventalign","{{cond}}_{{rep}}", "{chunk}_data.tsv"), chunk=chunk_list),
        summary_list=expand(join("results", module_name, "nanopolish_eventalign","{{cond}}_{{rep}}", "{chunk}_summary.tsv"), chunk=chunk_list)
    output:
        tsv=join("results", module_name, rule_name, "{cond}_{rep}_data.tsv"),
        summary=join("results", module_name, rule_name, "{cond}_{rep}_summary.tsv")
    log: join("logs", module_name, rule_name, "{cond}_{rep}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name),
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    script: "../scripts/eventalign_merge.py"
