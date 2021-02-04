# -*- coding: utf-8 -*-

##### Imports #####
from snakemake.shell import shell
import os

##### RUN SCRIPT FUNCTION #####
control_bam_list = snakemake.input.control_bam
test_bam_list = snakemake.input.test_bam
fasta = snakemake.input.fasta
res_bed = snakemake.output.res_bed
log = snakemake.log[0]
threads = snakemake.threads
opt = snakemake.params.opt

control_bam_opt = ""
for i in control_bam_list:
    control_bam_opt += f" -a {i}"

test_bam_opt = ""
for i in test_bam_list:
    test_bam_opt += f" -b {i}"

shell("differr -p {threads} {opt} {control_bam_opt} {test_bam_opt} -r {fasta} -o {res_bed} &> {log}")
