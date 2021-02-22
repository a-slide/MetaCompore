# -*- coding: utf-8 -*-

##### Imports #####
from os.path import join
from snakemake.logging import logger
from snakemake.utils import min_version
from glob import glob
min_version("5.30.0")

include: "workflow/rules/common.smk"

##### load config and sample sheets #####

logger.info("Loading and checking configuration file")
config = config_load_validate(configfile="config.yaml", schema="workflow/schemas/config.schema.yaml")

logger.info("Loading and checking sample file")
samples_df = samples_load_validate(samplefile="samples.tsv", schema="workflow/schemas/samples.schema.yaml")
replicates_list=list(samples_df["replicate"].unique())
condition_list=list(samples_df["condition"].unique())
logger.info(f"replicates found: {replicates_list}")
logger.info(f"condition found: {condition_list}")

##### Define all output files depending on config file #####

logger.info("Defining target files")
target_files=[]

# Add input target files
target_files.extend(expand(join("results", "alignment", "alignmemt_postfilter", "{cond}_{rep}.bam"), cond=condition_list, rep=replicates_list))
target_files.extend(expand(join("results", "alignment", "alignmemt_merge", "{cond}.bam"), cond=condition_list))
target_files.extend(expand(join("results", "resquiggling", "f5c_eventalign", "{cond}_{rep}_data.tsv"), cond=condition_list, rep=replicates_list))
target_files.extend(expand(join("results", "resquiggling", "f5c_eventalign", "{cond}_{rep}_summary.tsv"), cond=condition_list, rep=replicates_list))

if config.get("quality_control", None):
    logger.info("Defining target files for `quality_control` rules")
    target_files.extend(expand(join("results", "quality_control", "pycoQC", "pycoQC_{cond}_{rep}.json"), cond=condition_list, rep=replicates_list))
    target_files.extend(expand(join("results", "quality_control", "pycoQC", "pycoQC_{cond}_{rep}.html"), cond=condition_list, rep=replicates_list))

if config.get("nanocompore", None):
    logger.info("Defining target files for `nanocompore` rules")
    target_files.append(join("results", "final", "nanocompore_results_GMM_context_0.tsv"))
    target_files.append(join("results", "final", "nanocompore_results_GMM_context_2.tsv"))
    target_files.append(join("results", "final", "nanocompore_results_KS_dwell_context_0.tsv"))
    target_files.append(join("results", "final", "nanocompore_results_KS_dwell_context_2.tsv"))
    target_files.append(join("results", "final", "nanocompore_results_KS_intensity_context_0.tsv"))
    target_files.append(join("results", "final", "nanocompore_results_KS_intensity_context_2.tsv"))

if config.get("tombo", None):
    logger.info("Defining target files for `tombo` rules")
    target_files.append(join("results", "final", "tombo_results.tsv"))

if config.get("differr", None):
    logger.info("Defining target files for `differr` rules")
    target_files.append(join("results", "final", "differr_results.tsv"))

if config.get("eligos2", None):
    logger.info("Defining target files for `eligos2` rules")
    target_files.append(join("results", "final", "eligos2_results.tsv"))

if config.get("mines", None):
    logger.info("Defining target files for `mines` rules")
    target_files.append(join("results", "final", "mines_results.tsv"))

# if config.get("xpore", None):
#     target_files.append("xpore_out_files")

if config.get("epinano", None):
    target_files.append(join("results", "final", "epinano_results.tsv"))

##### Set main rule #####

rule all:
    input: target_files

##### Snakemake Include #####

include: "workflow/rules/input.smk"
include: "workflow/rules/basecalling.smk"
include: "workflow/rules/alignment.smk"
include: "workflow/rules/resquiggling.smk"
include: "workflow/rules/quality_control.smk"
include: "workflow/rules/nanocompore.smk"
include: "workflow/rules/tombo.smk"
include: "workflow/rules/differr.smk"
include: "workflow/rules/eligos2.smk"
include: "workflow/rules/mines.smk"
# include: "workflow/rules/xpore.smk"
include: "workflow/rules/epinano.smk"
