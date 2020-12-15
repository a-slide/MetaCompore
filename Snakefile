# -*- coding: utf-8 -*-

##### Imports #####
from os.path import join
from snakemake.logging import logger
from snakemake.utils import min_version
min_version("5.30.0")

include: "workflow/rules/common.smk"

##### load config and sample sheets #####

logger.info("Loading and checking configuration file")
config = config_load_validate(configfile="config.yaml", schema="workflow/schemas/config.schema.yaml")

logger.info("Loading and checking sample file")
samples_df = samples_load_validate(samplefile="samples.tsv", schema="workflow/schemas/samples.schema.yaml")
samples_list=list(samples_df.index)


##### Define all output files depending on config file #####

logger.info("Define conditional target files")
target_files=[]

# Add input target files
target_files.append(join("results", "input", "get_transcriptome", "transcriptome_reference.fa"))
target_files.extend(expand(join("results", "alignment", "minimap2_align", "{sample}.bam"), sample=samples_list))
target_files.extend(expand(join("results", "alignment", "pbt_alignmemt_filter", "{sample}.bam"), sample=samples_list))
target_files.extend(expand(join("results", "resquiggling", "f5c_eventalign", "{sample}.tsv"), sample=samples_list))

# if config.get("nanocompore", None):
#     target_files.append("nanocompore_out_files")
# if config.get("eligos2", None):
#     target_files.append("eligos2_out_files")
# if config.get("xpore", None):
#     target_files.append("xpore_out_files")
# if config.get("tombo", None):
#     target_files.append("tombo_out_files")
# if config.get("mines", None):
#     target_files.append("mines_out_files")
# if config.get("epinano", None):
#     target_files.append("epinano_out_files")
# if config.get("differr_nanopore_DRS", None):
#     target_files.append("differr_nanopore_DRS_out_files")

##### Set main rule #####

rule all:
    input: target_files

##### Snakemake Include #####

include: "workflow/rules/input.smk"
include: "workflow/rules/basecalling.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/resquiggling.smk"
# include: "workflow/rules/qc.smk"
# include: "workflow/rules/quantification.smk"
# include: "workflow/rules/nanocompore.smk"
# include: "workflow/rules/eligos2.smk"
# include: "workflow/rules/xpore.smk"
# include: "workflow/rules/tombo.smk"
# include: "workflow/rules/mines.smk"
# include: "workflow/rules/epinano.smk"
# include: "workflow/rules/differr_nanopore_DRS.smk"
