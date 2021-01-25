# -*- coding: utf-8 -*-

##### Imports #####
from snakemake.shell import shell
import os

##### RUN SCRIPT FUNCTION #####
test_tsv = ",".join(snakemake.input.test_tsv)
control_tsv = ",".join(snakemake.input.control_tsv)
fasta = snakemake.input.fasta
outdir = os.path.dirname(snakemake.output.res_tsv)
log = snakemake.log[0]
threads = snakemake.threads
opt = snakemake.params.opt

shell("nanocompore sampcomp -t {threads} {opt} -w -1 {control_tsv} -2 {test_tsv} --label1 control --label2 test -f {fasta} -o {outdir} &> {log}")
