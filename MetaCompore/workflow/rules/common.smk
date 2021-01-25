# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
import sys
from collections import *
import shutil

# Third party lib
import pandas as pd
from snakemake.logging import logger
from snakemake.io import load_configfile
from snakemake.utils import validate

##### Input files validation functions #####

def config_load_validate(configfile, schema):
    config = load_configfile(configfile)
    validate(config, schema=schema)
    logger.debug(config)
    return config

def samples_load_validate(samplefile, schema):
    samples_df = pd.read_csv (samplefile, comment="#", skip_blank_lines=True, sep="\t")
    validate(samples_df, schema=schema)
    samples_df.set_index("sample_id", drop=True, inplace=True)

    # Check that conditions and control and test
    conditions = sorted(list(samples_df["condition"].unique()))
    if conditions != ["control", "test"]:
        logger.error("Metacompore requires exactly 2 conditions: `control` and `test`")
        sys.exit()

    # Check that the conditions have the same replicates
    rep1 = sorted(list(samples_df["replicate"][samples_df["condition"]=="control"]))
    rep2 = sorted(list(samples_df["replicate"][samples_df["condition"]=="test"]))
    if rep1 != rep2:
        logger.error("The 2 condition groups requires the same number of replicates")
        sys.exit()

    logger.debug(samples_df)
    return samples_df

##### Getter functions #####

def get_fast5 (wildcards):
    res = samples_df[(samples_df["condition"]==wildcards.cond)&(samples_df["replicate"]==int(wildcards.rep))]
    return res["fast5_dir"][0]

def get_threads (config, rule_name, default=1):
    try:
        return config[rule_name]["threads"]
    except (KeyError, TypeError):
        logger.error("Could not find value `threads` for rule `{}` in config file".format(rule_name))
        return default

def get_opt (config, rule_name, default=""):
    try:
        return config[rule_name]["opt"]
    except KeyError:
        logger.error("Could not find value `opt` for rule `{}` in config file".format(rule_name))
        return default

def get_mem (config, rule_name, default=1000):
    try:
        return config[rule_name]["mem_mb"]
    except KeyError:
        logger.error("Could not find value `mem_mb` for rule `{}` in config file".format(rule_name))
        return default

##### Helper functions #####

def mkdir (fn, exist_ok=False):
    """ Create directory recursivelly. Raise IO error if path exist or if error at creation """
    try:
        os.makedirs (fn, exist_ok=exist_ok)
    except:
        raise NanocomporeError ("Error creating output folder `{}`".format(fn))

def access_file (fn, **kwargs):
    """ Check if the file is readable """
    return os.path.isfile (fn) and os.access (fn, os.R_OK)
