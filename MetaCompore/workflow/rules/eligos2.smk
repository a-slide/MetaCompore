# -*- coding: utf-8 -*-

##### Imports #####

# Std lib
from os.path import join

##### Rules #####
module_name = "eligos2"

rule_name="eligos2_fasta_to_bed"
rule eligos2_fasta_to_bed:
    input:
        fasta = rules.get_transcriptome.output.fasta
    output:
        bed = temp(join("results", module_name, rule_name, "transcriptome_reference.bed"))
    log: join("logs",module_name, rule_name, "out.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/metacompore_python:3.8.6"
    shell: "faidx --transform bed {input.fasta} > {output.bed}"

rule_name="eligos2_pair_diff_mod"
rule eligos2_pair_diff_mod:
    input:
        control_bam=join("results", "alignment", "alignmemt_merge", "control.bam"),
        test_bam=join("results", "alignment", "alignmemt_merge", "test.bam"),
        fasta=rules.get_transcriptome.output.fasta,
        bed=rules.eligos2_fasta_to_bed.output.bed
    output:
        outdir=directory(join("results", module_name, rule_name)),
        stats_tsv=join("results", module_name, rule_name, "test_vs_control_on_transcriptome_reference_combine.txt")
    log: join("logs",module_name, rule_name, "out.log")
    threads: get_threads(config, rule_name)
    params: opt=get_opt(config, rule_name)
    resources: mem_mb=get_mem(config, rule_name)
    container: "library://aleg/default/eligos2:2.0.0"
    shell: "eligos2 pair_diff_mod -t {threads} {params.opt} -cbam {input.control_bam} -tbam {input.test_bam} -ref {input.fasta} -reg {input.bed} -o {output.outdir} &> {log}"

# rule_name="eligos2_postprocess"
# rule eligos2_postprocess:

"""
usage: eligos2 pair_diff_mod [-tbam TEST_BAMS [TEST_BAMS ...]]
                            [-cbam CTRL_BAMS [CTRL_BAMS ...]] [-reg REGION]
                            [-ref REFERENCE] [-m MODEL]
                            [-bcf CDNA_BCF [CDNA_BCF ...]] [-p PREFIX]
                            [-o OUTDIR] [--sub_bam_dir SUB_BAM_DIR]
                            [--max_depth MAX_DEPTH] [--min_depth MIN_DEPTH]
                            [--esb ESB] [--oddR ODDR] [--pval PVAL]
                            [--adjPval ADJPVAL] [-t THREADS] [-h] [--force]

Required Arguments:
  -tbam TEST_BAMS [TEST_BAMS ...], --test_bams TEST_BAMS [TEST_BAMS ...]
                        Test RNA/DNA data in one or more sorted bam files
  -cbam CTRL_BAMS [CTRL_BAMS ...], --ctrl_bams CTRL_BAMS [CTRL_BAMS ...]
                        Control RNA/DNA data in one or more sorted bam files
  -reg REGION, --region REGION
                        BED(BED6/BED12) format of genes or regions for finding
                        modification.
  -ref REFERENCE, --reference REFERENCE
                        FA/FASTA format of the reference genome.

Optional Arguments:
  -m MODEL, --model MODEL
                        Path of rBEM5+2 model in JSON file format
  -bcf CDNA_BCF [CDNA_BCF ...], --cdna_bcf CDNA_BCF [CDNA_BCF ...]
                        Direct cDNA data in sorted bcf file(s).

Output Arguments:
  -p PREFIX, --prefix PREFIX
                        Set the output file prefix (default: None)
  -o OUTDIR, --outdir OUTDIR
                        Path of directory name to store output (default:
                        results)
  --sub_bam_dir SUB_BAM_DIR
                        Temporary directory path for storing subset of bam
                        files (default: tmp)

Criteria Arguments:
  --max_depth MAX_DEPTH
                        Maximum number of reads. default: (10000)
  --min_depth MIN_DEPTH
                        Minimum number of reads. default: (20)
  --esb ESB             Minimum cut-off for ratio of error at specific base
                        (ESB). default: (0.2)
  --oddR ODDR           Minimum cut-off for Odd ratio. default: (2.5)
  --pval PVAL           p-value cut-off. default: (0.05)
  --adjPval ADJPVAL     Adjusted p-value cut-off. default: (1)

Multiprocessing Arguments:
  -t THREADS, --threads THREADS
                        Number of threads (default: 3)

"""
