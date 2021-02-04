# -*- coding: utf-8 -*-

##### Imports #####
from snakemake.shell import shell
import os
from tempfile import TemporaryDirectory

##### RUN SCRIPT FUNCTION #####
control_bam = snakemake.input.control_bam
test_bam = snakemake.input.test_bam
fasta = snakemake.input.fasta
bed = snakemake.input.bed
res_tsv = snakemake.output.res_tsv
log = snakemake.log[0]
threads = snakemake.threads
opt = snakemake.params.opt
outdir = os.path.dirname(res_tsv)

with TemporaryDirectory() as tempdir:
    shell("eligos2 pair_diff_mod -t {threads} {opt} -cbam {control_bam} -tbam {test_bam} -ref {fasta} -reg {bed} --sub_bam_dir {tempdir} -o {outdir} &> {log}")
