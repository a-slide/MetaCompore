# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name="epinano"

# Scattergather config
scattergather:
    epinano_bam_split=10


if config["gpu_acceleration"]:
    guppy_container="library://aleg/default/ont_guppy:gpu-3.1.5"
else:
    guppy_container="library://aleg/default/ont_guppy:cpu-3.1.5"

rule_name="ont_guppy_epinano"
rule ont_guppy_epinano:
    input: fast5_dir=get_fast5
    output:
        seqsum=join("results", module_name, rule_name, "{cond}_{rep}","sequencing_summary.txt"),
        fastq_dir=directory(join("results", module_name, rule_name, "{cond}_{rep}"))
    log: join("logs", module_name, rule_name, "{cond}_{rep}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: guppy_container
    shell: "guppy_basecaller {params.opt} -i {input.fast5_dir} -s {output.fastq_dir} &> {log}"

rule_name="merge_fastq_epinano"
rule merge_fastq_epinano:
    input: fastq_dir=rules.ont_guppy_epinano.output.fastq_dir
    output: fastq=join("results", module_name, rule_name, "{cond}_{rep}.fastq")
    log: join("logs",module_name, rule_name, "{cond}_{rep}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/pybiotools:0.2.7"
    shell: "pyBioTools Fastq Filter {params.opt} -i {input.fastq_dir} -o {output.fastq} --verbose &> {log}"

rule_name="minimap2_align_epinano"
rule minimap2_align_epinano:
    input:
        idx=rules.minimap2_index.output.idx,
        fastq=rules.merge_fastq_epinano.output.fastq
    output:
        bam=temp(join("results", module_name, rule_name, "{cond}_{rep}.bam")),
        bam_index=temp(join("results", module_name, rule_name, "{cond}_{rep}.bam.bai"))
    log: join("logs",module_name, rule_name, "{cond}_{rep}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/minimap2:2.17"
    shell: "minimap2 -t {threads} {params.opt} {input.idx} {input.fastq} 2> {log} | \
            samtools view -bh 2>> {log} | \
            samtools sort > {output.bam} 2>> {log} \
            && samtools index {output.bam} 2>> {log}"

rule_name="alignmemt_prefilter_epinano"
rule alignmemt_prefilter_epinano:
    input:
        bam=rules.minimap2_align_epinano.output.bam,
        bam_index=rules.minimap2_align_epinano.output.bam_index
    output:
        bam=temp(join("results", module_name, rule_name, "{cond}_{rep}.bam")),
        bam_index=temp(join("results", module_name, rule_name, "{cond}_{rep}.bam.bai")),
        reads_index=temp(join("results", module_name, rule_name, "{cond}_{rep}.bam.idx.gz"))
    log: join("logs",module_name, rule_name, "{cond}_{rep}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/pybiotools:0.2.7"
    shell: "pyBioTools Alignment Filter {params.opt} -i {input.bam} -o {output.bam} --verbose &> {log}"

rule_name="min_ref_coverage_epinano"
rule min_ref_coverage_epinano:
    input:
        reads_index_list=expand(join("results", module_name, "alignmemt_prefilter_epinano", "{cond}_{rep}.bam.idx.gz"), cond=condition_list, rep=replicates_list)
    output:
        ref_list=join("results", module_name, rule_name, "valid_references_list.txt"),
    log: join("logs",module_name, rule_name, "out.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    script: "../scripts/min_ref_coverage.py"

rule_name="alignmemt_postfilter_epinano"
rule alignmemt_postfilter_epinano:
    input:
        bam=rules.alignmemt_prefilter_epinano.output.bam,
        bam_index=rules.alignmemt_prefilter_epinano.output.bam_index,
        ref_list=rules.min_ref_coverage_epinano.output.ref_list
    output:
        bam=join("results", module_name, rule_name, "{cond}_{rep}.bam"),
        bam_index=join("results", module_name, rule_name, "{cond}_{rep}.bam.bai"),
        reads_index=join("results", module_name, rule_name, "{cond}_{rep}.bam.idx.gz"),
        selected_reads_fn=join("results", module_name, rule_name, "{cond}_{rep}_selected_read_ids.txt")
    log: join("logs",module_name, rule_name, "{cond}_{rep}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/pybiotools:0.2.7"
    shell: "pyBioTools Alignment Filter {params.opt} -i {input.bam} --select_ref_fn {input.ref_list} -o {output.bam} -l {output.selected_reads_fn} --verbose &> {log}"

rule_name="alignmemt_merge_epinano"
rule alignmemt_merge_epinano:
    input:
        bam_list=expand(join("results", module_name, "alignmemt_postfilter_epinano", "{{cond}}_{rep}.bam"), rep=replicates_list),
        bam_index_list=expand(join("results", module_name, "alignmemt_postfilter_epinano", "{{cond}}_{rep}.bam.bai"), rep=replicates_list)
    output:
        bam=join("results", module_name, rule_name, "{cond}.bam"),
        bam_index=join("results", module_name, rule_name, "{cond}.bam.bai"),
    log: join("logs",module_name, rule_name, "{cond}.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/minimap2:2.17"
    shell: "samtools merge {params.opt} -@ {threads} {output.bam} {input.bam_list} &> {log} && samtools index {output.bam} &> {log}"


rule_name="epinano_splitbam"
rule epinano_splitbam:
    input:
        bam=rules.alignmemt_merge_epinano.output.bam
    output:
        temp(scatter.epinano_bam_split(join("results", module_name, rule_name, "{{cond}}_split/{{cond}}.{scatteritem}.bam")))
    container: "library://aleg/default/pybiotools:0.2.7"
    shell: "pyBioTools Alignment Split --index -i {input.bam} -n 10 --output_fn_list {output}"

rule_name="epinano_variants"
rule epinano_variants:
    input:
        bam=join("results", module_name, "epinano_splitbam", "{cond}_split/{cond}.{scatteritem}.bam"),
        fasta=rules.get_transcriptome.output.fasta,
        picard_dict=rules.generate_transcriptome_picard_index.output.picard_dict
    output:
        variants=temp(join("results", module_name, rule_name, "{cond}.{scatteritem}.plus_strand.per_site.5mer.csv"))
    log: join("logs",module_name, rule_name, "{cond}.{scatteritem}.log"),
    threads: get_threads(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/epinano:1.2.0"
    shell:
        """
        Epinano_Variants -R {input.fasta} -b {input.bam} -n {threads} -T t -s /usr/EpiNano/misc/sam2tsv.jar
        # This mv is needed because Epinano_Variants and Slide_Variants.py don't allow changing the output path
        mv $(dirname {input.bam})/{wildcards.cond}.{wildcards.scatteritem}.plus_strand.per.site.csv $(dirname {output.variants})/{wildcards.cond}.{wildcards.scatteritem}.plus_strand
        python /usr/EpiNano/misc/Slide_Variants.py $(dirname {output.variants})/{wildcards.cond}.{wildcards.scatteritem}.plus_strand 5
        """

rule_name="epinano_filter_rrach_variants"
rule epinano_filter_rrach_variants:
    input:
        variants=rules.epinano_variants.output.variants
    output:
        filteredvariants=temp(join("results", module_name, rule_name, "{cond}.{scatteritem}.plus_strand.per_site.5mer.csv.filtered"))
    log: join("logs",module_name, rule_name, "{cond}.{scatteritem}.log"),
    threads: get_threads(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/epinano:1.2.0"
    script: "../scripts/epinano_filter_kmers.py"

rule_name="epinano_gather_variants"
rule epinano_gather_variants:
    input:
        gather.epinano_bam_split(join("results", module_name, "epinano_filter_rrach_variants", "{{cond}}.{scatteritem}.plus_strand.per_site.5mer.csv.filtered"))
    output:
        filteredvariants=join("results", module_name, rule_name, "{cond}.all.plus_strand.per_site.5mer.csv")
    log: join("logs",module_name, rule_name, "{cond}.log"),
    threads: get_threads(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/epinano:1.2.0"
    shell: "cat {input} | awk 'NR==1 || !/^#/' > {output}"

rule_name="epinano_predict"
rule epinano_predict:
    input:
        variants=rules.epinano_gather_variants.output.filteredvariants
    output:
        predictions=temp(join("results", module_name, rule_name, "{cond}.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"))
    log: join("logs",module_name, rule_name, "{cond}.log"),
    threads: get_threads(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/epinano:1.2.0"
    shell: "Epinano_Predict -o $(dirname {output.predictions})/{wildcards.cond} -M /usr/EpiNano/models/rrach.q3.mis3.del3.linear.dump -p {input.variants} -cl 8,13,23"

rule_name="epinano_delta_variants"
rule epinano_delta_variants:
    input:
        control_variants=join("results", module_name, "epinano_gather_variants", "control.all.plus_strand.per_site.5mer.csv"),
        test_variants=join("results", module_name, "epinano_gather_variants", "test.all.plus_strand.per_site.5mer.csv")
    output:
        delta=join("results", module_name, rule_name, "test_control_delta.5mer.csv")
    log: join("logs",module_name, rule_name, "delta_variants.log"),
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/epinano:1.2.0"
    shell: "python /usr/EpiNano/misc/Epinano_make_delta.py {input.test_variants} {input.control_variants} {params.opt[min_cov]} 5 > {output.delta}"

rule_name="epinano_delta_predict"
rule epinano_delta_predict:
    input:
        delta_variants=rules.epinano_delta_variants.output.delta
    output:
        predictions=temp(join("results", module_name, rule_name, "test_control.DeltaMis3.DeltaDel3.DeltaQ3.MODEL.rrach.deltaQ3.deltaMis3.deltaDel3.linear.dump.csv"))
    log: join("logs",module_name, rule_name, "delta_variants.log"),
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/epinano:1.2.0"
    shell: "Epinano_Predict -o $(dirname {output.predictions})/test_control -M /usr/EpiNano/models/rrach.deltaQ3.deltaMis3.deltaDel3.linear.dump -p {input.delta_variants} -cl 7,12,22"



rule_name="epinano_postprocess"
rule epinano_postprocess:
    input:
        join("results", module_name, "epinano_predict", "test.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"),
        join("results", module_name, "epinano_predict", "control.q3.mis3.del3.MODEL.rrach.q3.mis3.del3.linear.dump.csv"),
        join("results", module_name, "epinano_delta_predict", "test_control.DeltaMis3.DeltaDel3.DeltaQ3.MODEL.rrach.deltaQ3.deltaMis3.deltaDel3.linear.dump.csv"),
    output:
        res_tsv=join("results", "final", "{}_results.tsv".format(module_name))
    wildcard_constraints:
        cond="[^.]"
    log: join("logs",module_name, rule_name, "out.log"),
    threads: get_threads(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/epinano:1.2.0"
    shell: "touch {output.res_tsv}"


#python ../../Epinano_Variants.py -R ref.fa -b wt.bam -n 6 -T t -s ../../misc/sam2tsv.jar
#python ../../Epinano_Variants.py -R ref.fa -b ko.bam -n 6 -T t -s ../../misc/sam2tsv.jar
#python ../../Epinano_Predict.py -o SVM_Predict -M /usr/EpiNano/models/rrach.q3.mis3.del3.linear.dump -p wt.plus_strand.per_site.5mer.csv -cl 8,13,23
#
#
#
#echo "generate delta-features"
#python ../../misc/Epinano_make_delta.py wt.plus_strand.per_site.5mer.csv ko.plus_strand.per_site.5mer.csv 5 > wt_ko_delta.5mer.csv
#echo "predict using pretrained SVM models with delta features"
#python ../../Epinano_Predict.py -o SVM_Predict_delta_features -M /usr/EpiNano//models/rrach.deltaQ3.deltaMis3.deltaDel3.linear.dump -p wt_ko_delta.5mer.csv -cl 7,12,22

#python $EPINANO_HOME/Epinano_Predict.py
#	--model q3.mis3.del3.MODEL.linear.model.dump
#	--predict some_sample.per_site.5mer.csv
#	--columns 8,13,23
#	--out_prefix some_sample.modification
#
#
#python $EPINANO_HOME/Epinano_Variants.py -n 6 -R reference.fasta -b sample.reads.bam -s /path/to/sam2tsv/sam2tsv.jar --type t
#
#
