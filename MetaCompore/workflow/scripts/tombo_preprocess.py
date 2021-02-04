# -*- coding: utf-8 -*-

##### Imports #####
from snakemake.shell import shell
import tempfile

##### RUN SCRIPT FUNCTION #####
input_fast5_dir = snakemake.input.fast5_dir
fastq = snakemake.input.fastq
fasta = snakemake.input.fasta
selected_reads_fn = snakemake.input.selected_reads_fn
output_fast5_dir = snakemake.output.fast5_dir
log = snakemake.log[0]
threads = snakemake.threads

with tempfile.TemporaryDirectory() as temp_dir:

    shell("echo '## fast5_subset ##' > {log}")
    shell("fast5_subset -t {threads} -i {input_fast5_dir} -s {temp_dir} -l {selected_reads_fn} -r &>> {log}")

    shell("echo '## multi_to_single_fast5 ##' >> {log}")
    shell("multi_to_single_fast5 -t {threads} -i {temp_dir} -s {output_fast5_dir} &>> {log}")

shell("echo '## tombo preprocess annotate_raw_with_fastqs ##' >> {log}")
shell("tombo preprocess annotate_raw_with_fastqs --fast5-basedir {output_fast5_dir} --fastq-filenames {fastq} --overwrite --processes {threads} &>> {log}")

shell("echo '## tombo resquiggle ##' >> {log}")
shell("tombo resquiggle --rna --overwrite --processes {threads} {output_fast5_dir} {fasta} &>> {log}")
