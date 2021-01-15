# -*- coding: utf-8 -*-

##### Imports #####
from snakemake.shell import shell

##### RUN SCRIPT FUNCTION #####
output = snakemake.output
input = snakemake.input
log = snakemake.log[0]
threads = snakemake.threads

shell("echo '## multi_to_single_fast5 ##' > {log}")
shell("multi_to_single_fast5 -t {threads} -i {input.fast5_dir} -s {output.fast5_dir} &>> {log}")
shell("echo '## tombo preprocess annotate_raw_with_fastqs ##' >> {log}")
shell("tombo preprocess annotate_raw_with_fastqs --fast5-basedir {output.fast5_dir} --fastq-filenames {input.fastq} --overwrite --processes {threads} &>> {log}")
shell("echo '## tombo resquiggle ##' >> {log}")
shell("tombo resquiggle --rna --overwrite --processes {threads} {output.fast5_dir} {input.fasta} &>> {log}")
# shell("touch {output.done}")
