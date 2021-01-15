# -*- coding: utf-8 -*-

##### Imports #####
from snakemake.shell import shell
import os

##### RUN SCRIPT FUNCTION #####
test_fast5 = snakemake.input.test_fast5
control_fast5 = snakemake.input.control_fast5
stats_h5 = snakemake.output.stats_h5
log = snakemake.log[0]
threads = snakemake.threads
opt = snakemake.params.opt

# remove tombo added sufix
stats_h5_dir = stats_h5[:-12]

shell("tombo detect_modifications level_sample_compare {opt} --processes {threads} --fast5-basedirs {test_fast5} --alternate-fast5-basedirs {control_fast5} --statistics-file-basename {stats_h5_dir} &> {log}")
